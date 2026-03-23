#pragma once //include this only once during compilation

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "Assembler.hpp"
#include "BoundaryConditions.hpp"

template <unsigned int Nne, unsigned int Nsd>
class NonlinearSolver{
    public:
        NonlinearSolver(double tol, unsigned int maxIncr, unsigned int maxIter);

        void solve(
            Eigen::VectorXd& u, //displacement vector, modified in place
            const Assembler<Nne, Nsd>& assembler, //provides Kglobal, Rglobal
            const BoundaryConditions<Nne>& bcs //provides dirischlet indexes and values
        );

    private:
        double tol_; //tolerance for convergence
        unsigned int maxIncr_; //maximum number of increments (timesteps)
        unsigned int maxIter_; //maximum number of iterations per increment
};

#include "NonLinearSolver.tpp" //include the implementation of the NonlinearSolver class
