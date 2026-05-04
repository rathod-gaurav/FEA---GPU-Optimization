#pragma once

#include "MaterialModel.hpp"

class MooneyRivlin : public MaterialModel{ //derived from MaterialModel class
    public:
        MooneyRivlin(double C10, double C01, double D1); //default constructor - takes material parameters C10 and D1

        void compute(
            const Eigen::Matrix3d& F, //deformation gradient
            Eigen::Matrix3d& P, //first Piola-Kirchhoff stress tensor
            Tensor4D& C_mat //material tangent stiffness matrix (4th order elasticity tensor in Voigt notation)
        ) const override; //override the pure virtual function from MaterialModel class

    private:
        double C10_; //material parameter C10
        double C01_; //material parameter C01
        double D1_; //material parameter D1
};