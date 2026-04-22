#pragma once //include this only once during compilation

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

class ConjugateGradientSolver{
    public:
        ConjugateGradientSolver(double tol, unsigned int maxIter);

        void solve(
            Eigen::VectorXd& x0, //initial guess of solution
            const Eigen::SparseMatrix<double>& A, // Ax = b
            const Eigen::VectorXd& b,
            int numThreads
        );

    private:
        double tol_;
        unsigned int maxIter_;
};