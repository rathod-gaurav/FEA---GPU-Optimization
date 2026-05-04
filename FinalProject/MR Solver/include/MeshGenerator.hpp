#pragma once //include this file only once during compilation
#include "Mesh.hpp"

template <unsigned int Nne>
class MeshGenerator{
    public:
        MeshGenerator( //default constructor for MeshGenerator class
            //domain parameters
            double x1_ll, double x1_ul,
            double x2_ll, double x2_ul,
            double x3_ll, double x3_ul,
            //mesh parameters
            unsigned int Nel_x1, unsigned int Nel_x2, unsigned int Nel_x3
        );

        Mesh<Nne> buildMesh() const; //function to build the mesh and return a Mesh object | this is a const member function and does not modify the state of the MeshGenerator object

    private: //these variables are private and can only be accessed within the MeshGenerator class
        //domain parameters
        double x1_ll_, x1_ul_; //the trailing underscore says - this belongs to the class and is not a local variable in the function | this is a common naming convention for class member variables
        double x2_ll_, x2_ul_;
        double x3_ll_, x3_ul_;
        //mesh parameters
        unsigned int Nel_x1_, Nel_x2_, Nel_x3_;
};

#include "MeshGenerator.tpp" //include the implementation of the member functions of the MeshGenerator class