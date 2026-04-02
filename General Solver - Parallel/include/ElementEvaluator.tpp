#pragma once

template <unsigned int Nne, unsigned int Nsd>
ElementEvaluator<Nne, Nsd>::ElementEvaluator(
    const Mesh<Nne>& mesh,
    const MaterialModel& material,
    const QuadratureRule& quadRule
) : 
    mesh_(mesh),
    material_(material),
    quadRule_(quadRule)
{}

template <unsigned int Nne, unsigned int Nsd>
Eigen::Matrix3d ElementEvaluator<Nne, Nsd>::computeJacobian(unsigned int e, double xi1, double xi2, double xi3) const{
    Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
        
    for(int A = 0 ; A < Nne ; A++){
        auto [basis_gradient_xi1, basis_gradient_xi2, basis_gradient_xi3] = ShapeFunction::basis_gradient(A, xi1, xi2, xi3);
        unsigned int Aglobal = mesh_.elements[e].node[A];
        J(0,0) += basis_gradient_xi1*mesh_.nodes[Aglobal].x1; //dx1/dxi1
        J(0,1) += basis_gradient_xi2*mesh_.nodes[Aglobal].x1; //dx1/dxi2
        J(0,2) += basis_gradient_xi3*mesh_.nodes[Aglobal].x1; //dx1/dxi3
        J(1,0) += basis_gradient_xi1*mesh_.nodes[Aglobal].x2; //dx2/dxi1
        J(1,1) += basis_gradient_xi2*mesh_.nodes[Aglobal].x2; //dx2/dxi2
        J(1,2) += basis_gradient_xi3*mesh_.nodes[Aglobal].x2; //dx2/dxi3
        J(2,0) += basis_gradient_xi1*mesh_.nodes[Aglobal].x3; //dx3/dxi1
        J(2,1) += basis_gradient_xi2*mesh_.nodes[Aglobal].x3; //dx3/dxi2
        J(2,2) += basis_gradient_xi3*mesh_.nodes[Aglobal].x3; //dx3/dxi3
    }
    return J;
}

template <unsigned int Nne, unsigned int Nsd>
Eigen::Matrix3d ElementEvaluator<Nne, Nsd>::computeGradU(const Eigen::VectorXd& u_e, Eigen::Vector3d (&dN_dx)[Nne]) const {
    Eigen::Matrix3d grad_u = Eigen::Matrix3d::Zero();
    //compute the gradient of the displacement field at the quadrature point using the basis function gradients and the nodal displacements
    for(int A = 0 ; A < Nne ; A++){
        Eigen::Vector3d& dNA_dx = dN_dx[A]; // #optimization - alias, no copy
        grad_u(0,0) += dNA_dx(0) * u_e(A*3 + 0); //du1/dx1
        grad_u(0,1) += dNA_dx(1) * u_e(A*3 + 0); //du1/dx2
        grad_u(0,2) += dNA_dx(2) * u_e(A*3 + 0); //du1/dx3

        grad_u(1,0) += dNA_dx(0) * u_e(A*3 + 1); //du2/dx1
        grad_u(1,1) += dNA_dx(1) * u_e(A*3 + 1); //du2/dx2
        grad_u(1,2) += dNA_dx(2) * u_e(A*3 + 1); //du2/dx3

        grad_u(2,0) += dNA_dx(0) * u_e(A*3 + 2); //du3/dx1
        grad_u(2,1) += dNA_dx(1) * u_e(A*3 + 2); //du3/dx2
        grad_u(2,2) += dNA_dx(2) * u_e(A*3 + 2); //du3/dx3
    }
    return grad_u;
}

template <unsigned int Nne, unsigned int Nsd>
void ElementEvaluator<Nne, Nsd>::computeElement(
    unsigned int e, //element index
    const Eigen::VectorXd& u_e, //element nodal displacements (Nne*3 x 1 vector)
    Eigen::MatrixXd& Klocal, //element stiffness matrix (Nne*3 x Nne*3 matrix)
    Eigen::VectorXd& Rlocal //element internal force vector (Nne*3 x 1 vector)
) const {
    
    const auto& quad_points = quadRule_.points;
    const auto& quad_weights = quadRule_.weights;
    unsigned int quadOrder = quad_points.size(); //number of quadrature points in each direction

    // Rlocal = Eigen::VectorXd::Zero(Nne * Nsd); //local residual vector for the element
    // Klocal = Eigen::MatrixXd::Zero(Nne * Nsd, Nne * Nsd); //local tangent stiffness matrix for the element

    #pragma omp parallel
    {
        Eigen::VectorXd Rlocal_priv = Eigen::VectorXd::Zero(Nne * Nsd); //local residual vector for the element
        Eigen::MatrixXd Klocal_priv = Eigen::MatrixXd::Zero(Nne * Nsd, Nne * Nsd); //local tangent stiffness matrix for the element

        Eigen::Matrix3d S_priv = Eigen::Matrix3d::Zero(); //Second Piola-Kirchhoff stress tensor
        Eigen::Matrix3d P_priv = Eigen::Matrix3d::Zero(); //First Piola-Kirchhoff stress tensor
        Eigen::MatrixXd C_mat_priv = Eigen::MatrixXd::Zero(9,9); //material tangent stiffness matrix in Voigt notation (3x3 block for each pair of nodes)

        //Gaussian quadrature loop
        #pragma omp for collapse(3) schedule(static) //parallelize the quadrature loop with OpenMP, and use reduction to safely accumulate contributions to Rlocal and Klocal from different threads
        for(int I = 0 ; I < quadOrder ; I++){
            for(int J = 0 ; J < quadOrder ; J++){
                for(int K = 0 ; K < quadOrder ; K++){
                    //Get the quadrature point coordinates and weights
                    double xi1 = quad_points[I]; 
                    double xi2 = quad_points[J];
                    double xi3 = quad_points[K];
                    double weight = quad_weights[I] * quad_weights[J] * quad_weights[K];

                    Eigen::Matrix3d Jac = computeJacobian(e, xi1, xi2, xi3); //compute the Jacobian matrix at the quadrature point
                    double JacDet = Jac.determinant(); //compute the determinant of the Jacobian
                    Eigen::Matrix3d JacInv = Jac.inverse(); //compute the inverse of the Jacobian
                    Eigen::Matrix3d JacInvT = JacInv.transpose(); //#optimization - transpose of the inverse Jacobian only once per quadrature point
                    double JacDetWeight = JacDet * weight; // #optimization - combine Jacobian determinant and quadrature weight to reduce redundant multiplications

                    alignas(32) Eigen::Vector3d dN_dx[Nne]; //#optimization - cache the gradient of the basis functions in global coordinates for all nodes at the quadrature point to avoid redundant calculations
                    for(unsigned int B = 0 ; B < Nne ; B++){
                        auto [dN_dxi1, dN_dxi2, dN_dxi3] = ShapeFunction::basis_gradient(B, xi1, xi2, xi3);
                        dN_dx[B] = JacInvT*Eigen::Vector3d(dN_dxi1, dN_dxi2, dN_dxi3); //gradient of the basis function in global coordinates
                    }

                    Eigen::Matrix3d grad_u = computeGradU(u_e, dN_dx); //compute the gradient of the displacement field at the quadrature point #optimization - pass cached dN_dx to avoid redundant calculations in computeGradU
                    Eigen::Matrix3d F = Eigen::Matrix3d::Identity() + grad_u; //deformation gradient

                    material_.compute(F, S_priv, P_priv, C_mat_priv); //compute the stress tensors E,S,P and material tangent stiffness matrix at the quadrature point using the material model
                    
                    for(int B = 0 ; B < Nne ; B++){//Loop to calculate Residual
                        // #optimization - use cached dN_dx to avoid redundant calculations, and use Eigen's fixed size segment
                        Rlocal_priv.segment<Nsd>(B*Nsd) += P_priv * dN_dx[B] * JacDetWeight; //contribution to the local residual vector 
                    }
                    
                    // cout << "Calculated Rlocal for element " << e+1 << "/" << Nel_t << "\r";
                    
                    for(int A = 0 ; A < Nne ; A++){//Loops to calculate tangent matrix
                        Eigen::Vector3d& dNA = dN_dx[A]; // #optimization - alias, no copy
                        
                        for(int B = 0 ; B < Nne ; B++){
                            Eigen::Vector3d& dNB = dN_dx[B]; // #optimization - alias, no copy

                            //Kgeometric
                            double Kgeo_scalar = (dNA.transpose() * S_priv * dNB)(0,0);
                            double KgeoAB =  Kgeo_scalar * JacDetWeight;
                            Klocal_priv.block<3,3>(3*A,3*B).diagonal().array() += KgeoAB; //#optimization - only add to the diagonal entries of the 3x3 block for nodes A,B since KgeoAB is a scalar times the identity matrix, and use Eigen's fixed size block to avoid redundant calculations in block indexing

                            // Correct KmatAB (3x3 block for nodes A,B) #optimization - compute the geometric stiffness contribution to the tangent matrix using matrix operations instead of triple nested loops

                            Eigen::Matrix3d Gmat = Eigen::Matrix3d::Zero();
                            for(int P = 0; P < 3; P++){
                                for(int M = 0; M < 3; M++){
                                    Eigen::Matrix3d C_block = C_mat_priv.block<3,3>(3*P, 3*M); //material tangent stiffness block for nodes P,M
                                    Gmat(P,M) = dNA.dot(C_block * dNB); //contribution to the geometric stiffness from nodes P,M
                                }
                            }

                            Eigen::Matrix3d KmatAB = (F * Gmat * F.transpose())*JacDetWeight; //material stiffness contribution to the tangent matrix for nodes A,B
                            Klocal_priv.block<3,3>(3*A,3*B) += KmatAB; //add the material stiffness contribution to the local tangent stiffness matrix for nodes A,B
                            
                        }
                    }
                    
                    // cout << "Calculated Klocal for element " << e+1 << "/" << Nel_t << "\r";
                    
                }
            }
        }

        #pragma omp critical //safely accumulate contributions from different threads to the local residual and tangent stiffness matrix
        {
            Rlocal += Rlocal_priv;
            Klocal += Klocal_priv;
        }
    }

}