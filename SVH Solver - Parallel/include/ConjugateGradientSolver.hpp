#pragma once //include this only once during compilation

#include <Eigen/Dense>
#include <Eigen/Sparse>

class ConjugateGradientSolver{
    public:
        ConjugateGradientSolver(double tol, unsigned int maxIter);

        void solve(
            Eigen::VectorXd& x0, //initial guess of solution
            Eigen::SparseMatrix<double>& A, // Ax = b
            Eigen::VectorXd& b,
            int numThreads
        );

    private:
        double tol_;
        unsigned int maxIter_;
};