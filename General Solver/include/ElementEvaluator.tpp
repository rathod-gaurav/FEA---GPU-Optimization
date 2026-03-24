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

    Eigen::Matrix3d S = Eigen::Matrix3d::Zero(); //Second Piola-Kirchhoff stress tensor
    Eigen::Matrix3d P = Eigen::Matrix3d::Zero(); //First Piola-Kirchhoff stress tensor
    Eigen::MatrixXd C_mat = Eigen::MatrixXd::Zero(9,9); //material tangent stiffness matrix in Voigt notation (3x3 block for each pair of nodes)

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

                material_.compute(F, S, P, C_mat); //compute the stress tensors E,S,P and material tangent stiffness matrix at the quadrature point using the material model

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

                        //Kgeometric
                        double Kgeo_scalar = (dNA_dx.transpose() * S * dNB_dx)(0,0);
                        Eigen::Matrix3d KgeoAB =  Kgeo_scalar * JacDet * weight * Eigen::Matrix3d::Identity();
                        Klocal.block<3,3>(3*A,3*B) += KgeoAB;


                        // Correct KmatAB (3x3 block for nodes A,B)
                        Eigen::Matrix3d KmatAB = Eigen::Matrix3d::Zero();
                        for(int i = 0; i < 3; i++){
                            for(int j = 0; j < 3; j++){
                                double val = 0.0;
                                for(int P = 0; P < 3; P++){
                                    for(int Q = 0; Q < 3; Q++){
                                        for(int M = 0; M < 3; M++){
                                            for(int N = 0; N < 3; N++){
                                                double C_PQMN = C_mat(3*P + Q, 3*M + N); //material tangent stiffness in Voigt notation
                                                val += F(i,P)*C_PQMN*F(j,M)*dNA_dx(Q)*dNB_dx(N);
                                            }
                                        }
                                    }
                                }
                                KmatAB(i,j) = val;
                            }
                        }
                        KmatAB *= JacDet * weight;
                        Klocal.block<3,3>(3*A,3*B) += KmatAB;
                    }
                }
                
                // cout << "Calculated Klocal for element " << e+1 << "/" << Nel_t << "\r";
                
            }
        }
    }


}