#pragma once //include this only once during compilation

#include <Eigen/Dense>

//Any class that inherits from this MaterialModel class must provide these methods

struct Tensor4D{ //struct to store 4th order tensor in contiguous array
    int d1, d2, d3, d4; //dimensions of each axis in the 4th order tensor
    std::vector<double> data; //contiguous array to store the tensor data

    //initializer list to initialize the dimensions and allocate memory for the data vector
    Tensor4D(int n1, int n2, int n3, int n4): d1(n1), d2(n2), d3(n3), d4(n4), data(n1*n2*n3*n4, 0.0){};
    
    //Tensor4D(int n1, int n2, int n3, int n4)
    ///This defines what you need to provide to create a tensor. You’re telling the compiler: "To make a 4D tensor, I need four integers representing the size of each dimension."

    //The Initializer List (The part after the :)
    // d1(n1): This sets our internal variable d1 to whatever value you passed in for n1.
    // data(n1 * n2 * n3 * n4): This is the most important part. It initializes the std::vector and tells it exactly how much memory to reserve immediately.

    //function to access tensor elements using 4 indices
    double& operator()(int i, int j, int k, int l){
        return data[i*d2*d3*d4 + j*d3*d4 + k*d4 + l]; //calculate the index in the data vector based on the 4D indices
    }
};

class MaterialModel{
    public:
        virtual void compute(
            const Eigen::Matrix3d& F, //deformation gradient
            Eigen::Matrix3d& P, //first Piola-Kirchhoff stress tensor
            Tensor4D& C_mat //material tangent stiffness matrix (4th order elasticity tensor in Voigt notation)
        ) const = 0; // = 0 means this is a pure virtual function. which meahc that MaterialModel itself cannot be instantiated, but any derived class from this must inplement this compute method

        virtual ~MaterialModel() = default; //virtual destructor to ensure proper cleanup of derived classes
        //this is required wherever polymorphism is used, i.e., when we have a base class pointer pointing to a derived class object. This ensures that when we delete the base class pointer, the destructor of the derived class is also called, preventing memory leaks and ensuring proper resource management.
};