#include "ConjugateGradientSolver.hpp"

ConjugateGradientSolver::ConjugateGradientSolver(double tol, unsigned int maxIter)
    : tol_(tol), maxIter_(maxIter)
{}

void ConjugateGradientSolver::solve(
    Eigen::VectorXd& x0, //initial guess of solution
    const Eigen::SparseMatrix<double>& A, // Ax = b
    const Eigen::VectorXd& b
){
    Eigen::VectorXd r_k = b - A*x0;

    Eigen::VectorXd p_k = r_k;

    unsigned int k = 0;
    while(k < maxIter_){
        Eigen::VectorXd Ap = A*p_k;
        Eigen::VectorXd alpha = (r_k.transpose()*r_k)/(p_k.transpose()*Ap); 
        double alpha_k = alpha(0);

        x0 += alpha_k*p_k; //update the solution in place

        Eigen::VectorXd r_kp1 = r_k - alpha_k*Ap;

        if(r_kp1.norm() < tol_){
            std::cout << "CG solve complete in " << k << " iterations" << std::endl;
            return;
        }
        else{
            Eigen::VectorXd beta = (r_kp1.transpose()*r_kp1)/(r_k.transpose()*r_k);
            double beta_k = beta(0);
            Eigen::VectorXd p_kp1 = r_kp1 + beta_k*p_k;
            r_k = r_kp1;
            p_k = p_kp1;
            k++;
        }
        if(k == maxIter_){
            std::cout << "CG did not converge in " << maxIter_ << " iterations. " << "Final residual = " << r_k.norm() << "\n";
        }
    }
}


void ConjugateGradientSolver::solve_parallel(
    Eigen::VectorXd& x0,          // initial guess, solution in-place
    const Eigen::SparseMatrix<double>& A_csc,      // KUU in Eigen default CSC format
    const Eigen::VectorXd& b,          // right-hand side (-RU)
    int numThreads
){
    const int n = b.size();

    // ── Convert CSC → CSR once per Newton iteration ───────────────────────────
    // CSC (column-major): SpMV outer loop over columns → scattered writes to
    //   out[row], causing data races when parallelized.
    // CSR (row-major):    SpMV outer loop over rows    → each row writes only
    //   to out[row], fully independent → safe to parallelize with OMP.
    // Cost: O(nnz) conversion, done once here, amortized over all CG iterations.
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(A_csc);

    // Raw pointers into Eigen's CSR storage — no copies, valid while A lives
    const int*    rowPtr = A.outerIndexPtr();  // row start indices, size n+1
    const int*    colIdx = A.innerIndexPtr();  // column indices, size nnz
    const double* vals   = A.valuePtr();       // values, size nnz

    // ── Working vectors ───────────────────────────────────────────────────────
    Eigen::VectorXd r(n), p(n), Ap(n);

    // ── Initial residual: r = b - A*x0 ─────────────────────────────────────────
    // PARALLEL: each row of the SpMV is independent (CSR)
    // PARALLEL: the subtraction b[i] - Ap[i] is per-element
    #pragma omp parallel for schedule(static) num_threads(numThreads)
    for(int row = 0; row < n; row++){
        double sum = 0.0;
        for(int k = rowPtr[row]; k < rowPtr[row+1]; k++)
            sum += vals[k] * x0(colIdx[k]);
        r(row) = b(row) - sum;   // r = b - A*x0, one row at a time
    }

    // p = r (initial search direction)
    p = r;

    // rTr = r^T * r  (will be reused across iterations — no recomputation)
    // PARALLEL: dot product reduction
    double rTr = 0.0;
    #pragma omp parallel for reduction(+:rTr) schedule(static) num_threads(numThreads)
    for(int i = 0; i < n; i++)
        rTr += r(i) * r(i);

    // ── CG iteration loop ──────────────────────────────────────────────────────
    for(unsigned int k = 0; k < maxIter_; k++){

        // ── SpMV: Ap = A * p ───────────────────────────────────────────────────
        // This is the dominant cost (~1.1M FLOPs, called once per iteration).
        // PARALLEL: each row independently computes one element of Ap.
        // No write conflicts — row i writes only to Ap(i).
        #pragma omp parallel for schedule(static) num_threads(numThreads)
        for(int row = 0; row < n; row++){
            double sum = 0.0;
            for(int kk = rowPtr[row]; kk < rowPtr[row+1]; kk++)
                sum += vals[kk] * p(colIdx[kk]);
            Ap(row) = sum;
        }

        // ── pAp = p^T * Ap  (denominator of alpha) ────────────────────────────
        // FIX 1: was "p.T * A * p" which cost 2 SpMVs.
        // Now: 1 SpMV (above) + 1 dot product = correct and cheaper.
        // PARALLEL: dot product reduction
        double pAp = 0.0;
        #pragma omp parallel for reduction(+:pAp) schedule(static) num_threads(numThreads)
        for(int i = 0; i < n; i++)
            pAp += p(i) * Ap(i);

        double alpha = rTr / pAp;

        // ── Update solution: x0 = x0 + alpha * p ────────────────────────────────
        // PARALLEL: each element independent
        #pragma omp parallel for schedule(static) num_threads(numThreads)
        for(int i = 0; i < n; i++)
            x0(i) += alpha * p(i);

        // ── Update residual: r = r - alpha * Ap ───────────────────────────────
        // FIX 2: was "r_kp1 = b - A*x0" which cost a full SpMV every iteration.
        // Now: just an axpy using the Ap we already computed — O(n) not O(nnz).
        // Mathematical identity: r_{k+1} = r_k - alpha * A*p_k
        // PARALLEL: each element independent
        double rTr_new = 0.0;
        #pragma omp parallel for reduction(+:rTr_new) schedule(static) num_threads(numThreads)
        for(int i = 0; i < n; i++){
            r(i) -= alpha * Ap(i);
            rTr_new += r(i) * r(i);   // accumulate new rTr in the same pass
            // Combining the axpy and the dot product in one loop saves
            // a full pass over n elements — better cache utilization.
        }

        // ── Convergence check ──────────────────────────────────────────────────
        if(std::sqrt(rTr_new) < tol_){
            std::cout << "CG solve complete in " << k << " iterations" << std::endl;
            return;
        }

        // ── beta = rTr_new / rTr  (direction update scalar) ───────────────────
        // FIX 3: was recomputing r.T*r and r_old.T*r_old separately.
        // rTr is carried from the previous iteration — no extra dot product.
        double beta = rTr_new / rTr;
        rTr = rTr_new;  // carry forward for next iteration

        // ── Update search direction: p = r + beta * p ─────────────────────────
        // PARALLEL: each element independent
        #pragma omp parallel for schedule(static) num_threads(numThreads)
        for(int i = 0; i < n; i++)
            p(i) = r(i) + beta * p(i);
    }

    std::cout << "CG did not converge in " << maxIter_ << " iterations. " << "Final residual = " << std::sqrt(rTr) << "\n";
}