#include "StVenantKirchhoff.hpp"

StVenantKirchhoff::StVenantKirchhoff( //assign the function parameters to the class member variables using an initializer list
    double lambda,
    double mu
): lambda_(lambda), mu_(mu)
{}

void StVenantKirchhoff::compute(
    const Eigen::Matrix3d& F, //deformation gradient
    Eigen::Matrix3d& S, //second Piola-Kirchhoff stress tensor
    Eigen::Matrix3d& P //first Piola-Kirchhoff stress tensor
) const {

    Eigen::Matrix3d E = 0.5 * (F.transpose() * F - Eigen::Matrix3d::Identity()); //Green-Lagrange strain
    S = 2*mu_*E + lambda_*E.trace()*Eigen::Matrix3d::Identity(); //second Piola-Kirchhoff stress using St. Venant-Kirchhoff model
    P = F * S; //first Piola-Kirchhoff stress
}

void StVenantKirchhoff::computeCmat(
    Eigen::MatrixXd& C_mat //because material tangent matrix is constant for St. Venant Kirchhoff, we separate it out from the compute function so that it doesn't get calculated on every element
) const {
    // C_mat = Eigen::MatrixXd::Zero(9,9); //material tangent stiffness matrix (4th order elasticity tensor in Voigt notation)
    for(int P = 0; P < 3; P++){
        for(int Q = 0; Q < 3; Q++){
            for(int M = 0; M < 3; M++){
                for(int N = 0; N < 3; N++){
                    double C_PQMN = lambda_*(P==Q ? 1:0)*(M==N ? 1:0) + mu_*((P==M ? 1:0)*(Q==N ? 1:0) + (P==N ? 1:0)*(Q==M ? 1:0));
                    C_mat(P*3+Q, M*3+N) = C_PQMN;
                }
            }
        }
    }
}
