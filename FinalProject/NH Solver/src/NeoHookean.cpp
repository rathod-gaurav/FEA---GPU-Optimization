#include "NeoHookean.hpp"

NeoHookean::NeoHookean(
    double C10,
    double D1
): C10_(C10), D1_(D1)
{}

void NeoHookean::compute(
    const Eigen::Matrix3d& F, //deformation gradient
    Eigen::Matrix3d& P, //first Piola-Kirchhoff stress tensor
    Tensor4D& C_mat
) const {
    double JJ = F.determinant();
    Eigen::Matrix3d C = F.transpose() * F; //right Cauchy-Green deformation tensor
    double I1 = C.trace(); //first invariant of C

    Eigen::Matrix3d FinvT = F.inverse().transpose();

    double vol_multiplier = (2.0/D1_)*JJ;
    double iso_multiplier = C10_*std::pow(JJ, -2.0/3.0);

    // P = vol_multiplier*(JJ-1)*FinvT;
    P = vol_multiplier*(JJ-1)*FinvT + 2*iso_multiplier*(F - I1*(1.0/3.0)*FinvT); //total first Piola-Kirchhoff stress

    for(int i = 0 ; i < 3 ; i++){
        for(int J = 0 ; J < 3 ; J++){
            double vol_term1 = vol_multiplier*(2*JJ - 1)*FinvT(i,J);
            double iso_term1 = iso_multiplier*(-4.0/3.0)*F(i,J); 
            double iso_term2 = iso_multiplier*(-4.0/3.0)*FinvT(i,J);
            double iso_term3 = iso_multiplier*(4.0/9.0)*I1*FinvT(i,J);
            for(int j = 0 ; j < 3 ; j++){
                double vol_term2 = vol_multiplier*(JJ-1)*FinvT(j,J);
                double iso_term4 = iso_multiplier*(2.0/3.0)*I1*FinvT(j,J);
                double iso_term5 = iso_multiplier*2*(i==j ? 1:0);
                for(int K = 0 ; K < 3 ; K++){
                    C_mat(i,J,j,K) = vol_term1*FinvT(j,K) - vol_term2*FinvT(i,K) + iso_term1*FinvT(j,K) + iso_term2*F(j,K) + iso_term3*FinvT(j,K) + iso_term4*FinvT(i,K) + iso_term5*(J==K ? 1:0);
                }
            }
        }
    }
}