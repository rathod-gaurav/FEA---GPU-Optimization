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
Eigen::Matrix3d ElementEvaluator<Nne, Nsd>::computeGradU(const Eigen::VectorXd& u_e, double xi1, double xi2, double xi3, Eigen::Matrix3d& JacInv) const {
    Eigen::Matrix3d grad_u = Eigen::Matrix3d::Zero();
    //compute the gradient of the displacement field at the quadrature point using the basis function gradients and the nodal displacements
    for(int A = 0 ; A < Nne ; A++){
        auto [dN_dxi1, dN_dxi2, dN_dxi3] = ShapeFunction::basis_gradient(A, xi1, xi2, xi3);
        Eigen::Vector3d dN_dx = JacInv.transpose()*Eigen::Vector3d(dN_dxi1, dN_dxi2, dN_dxi3);
        grad_u(0,0) += dN_dx[0] * u_e(A*3 + 0); //du1/dx1
        grad_u(0,1) += dN_dx[1] * u_e(A*3 + 0); //du1/dx2
        grad_u(0,2) += dN_dx[2] * u_e(A*3 + 0); //du1/dx3

        grad_u(1,0) += dN_dx[0] * u_e(A*3 + 1); //du2/dx1
        grad_u(1,1) += dN_dx[1] * u_e(A*3 + 1); //du2/dx2
        grad_u(1,2) += dN_dx[2] * u_e(A*3 + 1); //du2/dx3

        grad_u(2,0) += dN_dx[0] * u_e(A*3 + 2); //du3/dx1
        grad_u(2,1) += dN_dx[1] * u_e(A*3 + 2); //du3/dx2
        grad_u(2,2) += dN_dx[2] * u_e(A*3 + 2); //du3/dx3
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
    Rlocal = Eigen::VectorXd::Zero(Nne * Nsd); //local residual vector for the element
    Klocal = Eigen::MatrixXd::Zero(Nne * Nsd, Nne * Nsd); //local tangent stiffness matrix for the element

    const auto& quad_points = quadRule_.points;
    const auto& quad_weights = quadRule_.weights;
    unsigned int quadOrder = quad_points.size(); //number of quadrature points in each direction

    Eigen::Matrix3d P = Eigen::Matrix3d::Zero(); //First Piola-Kirchhoff stress tensor
    Tensor4D C_mat(3,3,3,3); //material tangent stiffness matrix in Voigt notation (3x3 block for each pair of nodes)

    //Gaussian quadrature loop
    
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

                Eigen::Matrix3d grad_u = computeGradU(u_e, xi1, xi2, xi3, JacInv); //compute the gradient of the displacement field at the quadrature point
                Eigen::Matrix3d F = Eigen::Matrix3d::Identity() + grad_u; //deformation gradient

                if(F.determinant() <= 0){
                    std::cerr << "Warning: Negative F determinant at element " << e << " quadrature point (" << xi1 << ", " << xi2 << ", " << xi3 << ")." << std::endl;
                    throw std::runtime_error("Negative F determinant encountered. Simulation cannot proceed.");
                }

                material_.compute(F, P, C_mat); //compute the stress tensors E,S,P and material tangent stiffness matrix at the quadrature point using the material model
                
                // Eigen::Matrix3d b = F*F.transpose(); //left Cauchy-Green deformation tensor for debugging
                std::cout << "---------------------------------------" << std::endl;
                std::cout << "Element: " << e << ", Quadrature Point: (" << xi1 << ", " << xi2 << ", " << xi3 << ")\n";
                std::cout << "Jacobian:\n" << Jac << "\n";
                std::cout << "Deformation Gradient:\n" << F << "\n";
                std::cout << "JJ:\n" << F.determinant() << "\n";
                // std::cout << "b:'\n" << b << "\n";
                // std::cout << "dev(b):\n" << b - (1.0/3.0)*b.trace()*Eigen::Matrix3d::Identity() << "\n";
                std::cout << "sigma:\n" << (1/F.determinant())*P*F.transpose() << "\n";
                std::cout << "---------------------------------------" << std::endl;
                
                for(int B = 0 ; B < Nne ; B++){//Loop to calculate Residual
                    auto [dN_dxi1, dN_dxi2, dN_dxi3] = ShapeFunction::basis_gradient(B, xi1, xi2, xi3);
                    Eigen::Vector3d dN_dx = JacInv.transpose()*Eigen::Vector3d(dN_dxi1, dN_dxi2, dN_dxi3); //gradient of the basis function in global coordinates
                    Rlocal.segment(B*Nsd, Nsd) += P * dN_dx * weight * JacDet; //contribution to the local residual vector
                }
                
                // cout << "Calculated Rlocal for element " << e+1 << "/" << Nel_t << "\r";
                 
                for(int A = 0 ; A < Nne ; A++){//Loops to calculate tangent matrix
                    auto [dNA_dxi1, dNA_dxi2, dNA_dxi3] = ShapeFunction::basis_gradient(A, xi1, xi2, xi3);
                    Eigen::Vector3d dNA_dx = JacInv.transpose()*Eigen::Vector3d(dNA_dxi1, dNA_dxi2, dNA_dxi3);
                    
                     
                    for(int B = 0 ; B < Nne ; B++){

                        auto [dNB_dxi1, dNB_dxi2, dNB_dxi3] = ShapeFunction::basis_gradient(B, xi1, xi2, xi3);
                        Eigen::Vector3d dNB_dx = JacInv.transpose()*Eigen::Vector3d(dNB_dxi1, dNB_dxi2, dNB_dxi3);
                        
                        Eigen::Matrix3d K_AB = Eigen::Matrix3d::Zero(); //tangent stiffness matrix contribution from node A to node B

                        for(int i = 0 ; i < 3 ; i++){
                            for(int J = 0 ; J < 3 ; J++){
                                for(int j = 0 ; j < 3 ; j++){
                                    for(int KK = 0 ; KK < 3 ; KK++){
                                        K_AB(i,j) += C_mat(i,J,j,KK)*dNA_dx(J)*dNB_dx(KK)*JacDet*weight;
                                    }
                                }
                            }
                        }

                        Klocal.block<3,3>(3*A,3*B) += K_AB;
                    }
                }
                
                // cout << "Calculated Klocal for element " << e+1 << "/" << Nel_t << "\r";
                
            }
        }
    }


}