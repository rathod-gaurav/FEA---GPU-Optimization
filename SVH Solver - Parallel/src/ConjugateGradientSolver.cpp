#include "ConjugateGradientSolver.hpp"

ConjugateGradientSolver::ConjugateGradientSolver(double tol, unsigned int maxIter)
    : tol_(tol), maxIter_(maxIter)
{}

void ConjugateGradientSolver::solve(
    Eigen::VectorXd& x0, //initial guess of solution
    const Eigen::SparseMatrix<double>& A, // Ax = b
    const Eigen::VectorXd& b,
    int numThreads
){
    Eigen::VectorXd r_k = b - A*x0;

    Eigen::VectorXd p_k = r_k;

    unsigned int k = 0;
    while(k < maxIter_){
        Eigen::VectorXd alpha = (r_k.transpose()*r_k)/(p_k.transpose()*A*p_k); 
        double alpha_k = alpha(0);

        x0 += alpha_k*p_k; //update the solution in place

        Eigen::VectorXd r_kp1 = b - A*x0;

        if(r_kp1.norm() < tol_){
            std::cout << "cg solve complete in " << k << " iterations" << std::endl;
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
    }
}