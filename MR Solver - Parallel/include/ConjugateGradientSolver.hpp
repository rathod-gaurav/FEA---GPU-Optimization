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
            const Eigen::VectorXd& b
        );

        void solve_parallel(
            Eigen::VectorXd& x0,          // initial guess, solution in-place
            const Eigen::SparseMatrix<double>& A_csc,      // KUU in Eigen default CSC format
            const Eigen::VectorXd& b,          // right-hand side (-RU)
            int numThreads
        );

    private:
        double tol_;
        unsigned int maxIter_;
};