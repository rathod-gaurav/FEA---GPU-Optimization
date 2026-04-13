#pragma once //include this only once during compilation

template <unsigned int Nne, unsigned int Nsd>
NonlinearSolver<Nne, Nsd>::NonlinearSolver(double tol, unsigned int maxIncr, unsigned int maxIter)
    : tol_(tol), maxIncr_(maxIncr), maxIter_(maxIter)
{}

template <unsigned int Nne, unsigned int Nsd>
void NonlinearSolver<Nne, Nsd>::solve(
    Eigen::VectorXd& u, //displacement vector, modified in place
    const Assembler<Nne, Nsd>& assembler, //provides Kglobal, Rglobal
    int outerThreads, //number of threads for parallel execution //goes to Assembler.assembleSystem
    const BoundaryConditions<Nne>& bcs, //provides dirischlet indexes and values
    std::function<void(unsigned int, unsigned int, double)> iterCallback
){
    Eigen::VectorXd Rglobal, RU; //global residual vector
    Eigen::SparseMatrix<double> Kglobal, KUU, KUD; //global stiffness matrix
    for(unsigned int incr = 0; incr < maxIncr_; incr++){
        double incrFraction = (incr+1)/static_cast<double>(maxIncr_); //factor to scale dirischlet values for current incr
        bcs.applyToSolution(u, incrFraction); //apply dirischlet boundary conditions to the solution vector for the current incr

        for(unsigned int iter = 0; iter < maxIter_; iter++){
            
            assembler.assembleSystem(u, Kglobal, Rglobal, outerThreads); //assemble the global stiffness matrix and residual vector based on the current solution vector
            
            assembler.partition(Kglobal, Rglobal, bcs, KUU, KUD, RU); //partition the global stiffness matrix and residual vector into submatrices/vectors corresponding to unknown and dirischlet degrees of freedom

            // solve the linear system
            // std::cout << "Initilising solver for incr " << incr+1 << ", iteration " << iter+1 << "\n";

            // Eigen::FullPivLU<Eigen::MatrixXd> solver(KUU);
            Eigen::SparseLU<Eigen::SparseMatrix<double>> linear_solver;
            linear_solver.analyzePattern(KUU);
            linear_solver.factorize(KUU);
            if(linear_solver.info() != Eigen::Success) {
                std::cout << "Decomposition failed for incr " << incr+1 << ", iteration " << iter+1 << "\n";
                return;
            }

            // std::cout << "Initilised solver for incr " << incr+1 << ", iteration " << iter+1 << "\n";

            Eigen::VectorXd duU = linear_solver.solve(-RU); //solve for the incral displacements at the unknown degrees of freedom
            
            // std::cout << "Solved for incr " << incr+1 << ", iteration " << iter+1 << "\n";

            //construct full du vector including known values at dirischlet boundary
            const auto& unknownIndexes = bcs.getUnknownIndexes();
            for(int i = 0 ; i < unknownIndexes.size() ; i++){
                u(unknownIndexes[i]) += duU(i);
            }

            //check residual norm for convergence
            double residualNorm = RU.norm();
            std::cout << "incr: " << incr+1 << ", Iteration: " << iter+1 << "\n";
            std::cout << "Modified residual norm: " << residualNorm << "\n"; //print the norm of the modified residual to monitor convergence of the unknown degrees of freedom
            std::cout << "-----------------------------------" << "\n";
            
            if(iterCallback){
                iterCallback(incr, iter, residualNorm); //call the iteration callback function if provided
            }

            if(residualNorm < tol_){
                std::cout << "Convergence achieved for incr " << incr+1 << " in " << iter+1 << " iterations." << "\n";
                break; 
            }

        }
    }
}