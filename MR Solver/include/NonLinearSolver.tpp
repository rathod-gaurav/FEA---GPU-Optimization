#pragma once //include this only once during compilation

template <unsigned int Nne, unsigned int Nsd>
NonlinearSolver<Nne, Nsd>::NonlinearSolver(double tol, unsigned int maxIncr, unsigned int maxIter)
    : tol_(tol), maxIncr_(maxIncr), maxIter_(maxIter)
{}

template <unsigned int Nne, unsigned int Nsd>
void NonlinearSolver<Nne, Nsd>::numericalTangentCheck(
    const Eigen::VectorXd& u,
    const Assembler<Nne, Nsd>& assembler,
    const BoundaryConditions<Nne>& bcs,
    const Eigen::SparseMatrix<double>& KUU_analytical)
{
    const auto& unknownIndexes = bcs.getUnknownIndexes();
    int nFree = unknownIndexes.size();

    Eigen::VectorXd Rglobal, RU_plus, RU_minus;
    Eigen::SparseMatrix<double> Kglobal, KUU_dummy, KUD_dummy;

    // Get baseline R at current u
    assembler.assembleSystem(u, Kglobal, Rglobal);
    assembler.partition(Kglobal, Rglobal, bcs, KUU_dummy, KUD_dummy, RU_minus);
    Eigen::VectorXd RU_base = RU_minus;

    double eps = 1e-7;  // perturbation size
    Eigen::MatrixXd KUU_numerical = Eigen::MatrixXd::Zero(nFree, nFree);

    for(int j = 0; j < nFree; j++){
        Eigen::VectorXd u_pert = u;
        u_pert(unknownIndexes[j]) += eps;

        assembler.assembleSystem(u_pert, Kglobal, Rglobal);
        assembler.partition(Kglobal, Rglobal, bcs, KUU_dummy, KUD_dummy, RU_plus);

        // Numerical column of K
        KUU_numerical.col(j) = (RU_plus - RU_base) / eps;
    }

    // Compare with analytical
    Eigen::MatrixXd KUU_analytical_dense = KUU_analytical;
    Eigen::MatrixXd diff = KUU_numerical - KUU_analytical_dense;

    double rel_error = diff.norm() / KUU_analytical_dense.norm();

    std::cout << "=== Numerical Tangent Check ===\n";
    std::cout << "||K_numerical||   = " << KUU_numerical.norm() << "\n";
    std::cout << "||K_analytical||  = " << KUU_analytical_dense.norm() << "\n";
    std::cout << "||K_num - K_ana|| = " << diff.norm() << "\n";
    std::cout << "Relative error    = " << rel_error << "\n";

    // Find worst column (which DOF has biggest K error)
    for(int j = 0; j < nFree; j++){
        double col_err = diff.col(j).norm();
        if(col_err > 1e-3 * KUU_analytical_dense.col(j).norm()){
            std::cout << "  Large error in col " << j
                      << " (global DOF " << unknownIndexes[j] << ")"
                      << "  err=" << col_err << "\n";
        }
    }
    std::cout << "================================\n";
}


template <unsigned int Nne, unsigned int Nsd>
void NonlinearSolver<Nne, Nsd>::solve(
    Eigen::VectorXd& u, //displacement vector, modified in place
    const Assembler<Nne, Nsd>& assembler, //provides Kglobal, Rglobal
    const BoundaryConditions<Nne>& bcs, //provides dirischlet indexes and values
    std::function<void(unsigned int, unsigned int, double)> iterCallback
){
    Eigen::VectorXd Rglobal, RU; //global residual vector
    Eigen::SparseMatrix<double> Kglobal, KUU, KUD; //global stiffness matrix

    for(unsigned int incr = 0; incr < maxIncr_; incr++){
        double R0_norm = 0.0; //initial residual norm for convergence monitoring
        double incrFraction = (incr+1)/static_cast<double>(maxIncr_); //factor to scale dirischlet values for current incr

        bcs.applyToSolution(u, incrFraction); //apply dirischlet boundary conditions to the solution vector for the current incr

        for(unsigned int iter = 0; iter < maxIter_; iter++){
            // std::cout << "Assembling system for incr " << incr+1 << ", iteration " << iter+1 << "\n";
            
            assembler.assembleSystem(u, Kglobal, Rglobal); //assemble the global stiffness matrix and residual vector based on the current solution vector
            
            // std::cout << "displacement vector u:\n" << u << "\n";

            assembler.partition(Kglobal, Rglobal, bcs, KUU, KUD, RU); //partition the global stiffness matrix and residual vector into submatrices/vectors corresponding to unknown and dirischlet degrees of freedom

            // if(incr == 0 && iter == 0){
            //     numericalTangentCheck(u, assembler, bcs, KUU);
            // }

            



            // solve the linear system
            // std::cout << "Initilising solver for incr " << incr+1 << ", iteration " << iter+1 << "\n";
            
            // for (int i = 0; i < Kglobal.outerSize(); ++i) {
            //     if (Kglobal.innerVector(i).nonZeros() == 0) {
            //         std::cerr << "Column " << i << " is completely empty!" << std::endl;
            //         return;
            //     }
            // }

            // double max_diag = KUU.diagonal().maxCoeff();
            // double min_diag = KUU.diagonal().minCoeff();
            // std::cout << "Stiffness Ratio: " << max_diag / min_diag << std::endl;

            // Eigen::FullPivLU<Eigen::MatrixXd> solver(KUU);
            Eigen::SparseLU<Eigen::SparseMatrix<double>> linear_solver;
            linear_solver.analyzePattern(KUU);
            linear_solver.factorize(KUU);
            if(linear_solver.info() != Eigen::Success) {
                std::cout << "Decomposition failed for incr " << incr+1 << ", iteration " << iter+1 << "\n";
                std::cout << linear_solver.lastErrorMessage() << "\n";

                return;
            }

            // std::cout << "Initilised solver for incr " << incr+1 << ", iteration " << iter+1 << "\n";
            Eigen::VectorXd duU = linear_solver.solve(-RU); //solve for the incral displacements at the unknown degrees of freedom        
            
            // double R_norm_before = RU.norm();
            // const auto& unknownIndexes = bcs.getUnknownIndexes();

            // double scale_factor = 1.0;
            // bool line_search_success = false;

            // for(int s_iter = 0; s_iter < 10; s_iter++){
                
            //     // Reset u_trial to current u at every iteration
            //     Eigen::VectorXd u_trial = u;  // ← move inside loop
                
            //     // Apply scaled step
            //     for(int i = 0; i < unknownIndexes.size(); i++){
            //         u_trial(unknownIndexes[i]) += scale_factor * duU(i);
            //     }
                
            //     // Evaluate residual at trial point
            //     assembler.assembleSystem(u_trial, Kglobal, Rglobal);
            //     assembler.partition(Kglobal, Rglobal, bcs, KUU, KUD, RU);
            //     double R_norm_trial = RU.norm();
                
            //     std::cout << "  LS iter " << s_iter 
            //             << " alpha=" << scale_factor 
            //             << " ||R||=" << R_norm_trial << "\n";
                
            //     if(R_norm_trial < R_norm_before){
            //         u = u_trial;  // accept
            //         line_search_success = true;
            //         std::cout << "  LS accepted: alpha=" << scale_factor << "\n";
            //         break;
            //     }
                
            //     scale_factor *= 0.5;
            // }

            // if(!line_search_success){
            //     // Accept smallest step rather than stagnating
            //     Eigen::VectorXd u_trial = u;
            //     for(int i = 0; i < unknownIndexes.size(); i++){
            //         u_trial(unknownIndexes[i]) += scale_factor * duU(i);
            //     }
            //     u = u_trial;
            //     std::cout << "  LS failed — accepting smallest step\n";
            // }


            
            // std::cout << "Solved for incr " << incr+1 << ", iteration " << iter+1 << "\n";

            // construct full du vector including known values at dirischlet boundary
            const auto& unknownIndexes = bcs.getUnknownIndexes();
            for(int i = 0 ; i < unknownIndexes.size() ; i++){
                u(unknownIndexes[i]) += duU(i);
            }

            //check residual norm for convergence
            if(iter == 0){
                R0_norm = RU.norm(); //store the initial residual norm for convergence monitoring
            }
            double residualNorm = RU.norm()/R0_norm;
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