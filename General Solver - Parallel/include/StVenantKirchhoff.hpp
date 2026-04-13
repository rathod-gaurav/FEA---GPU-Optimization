#pragma once //include this only once during compilation

#include "MaterialModel.hpp"

class StVenantKirchhoff : public MaterialModel{ //derived from MaterialModel class
    public: 
        StVenantKirchhoff(double lambda, double mu); //default constructor - takes Lame' parameters

        void compute(
            const Eigen::Matrix3d& F, //deformation gradient
            Eigen::Matrix3d& S, //second Piola-Kirchhoff stress tensor
            Eigen::Matrix3d& P //first Piola-Kirchhoff stress tensor
        ) const override; //override the pure virtual function from MaterialModel class

        void computeCmat(Eigen::MatrixXd& C_mat) const override; //material tangent stiffness matrix (4th order elasticity tensor in Voigt notation)

    private:
        double lambda_; //Lame' parameter lambda
        double mu_; //Lame' parameter mu
};
