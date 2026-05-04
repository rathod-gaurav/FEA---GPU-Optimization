#pragma once //include this file only once during compilation

#include <tuple>

class ShapeFunction{
    public:
        //when "static" is used - it means that the function belongs to the class itself, not an instance. This means you can call the function without creating an object of the class. For example, you can call ShapeFunction::xi_at_node(node) directly without needing to instantiate a ShapeFunction object.
        static std::tuple<double,double,double> xi_at_node(unsigned int node); //function to return xi1, xi2, and xi3 for given node A

        static double basis_function(unsigned int node, double xi1, double xi2, double xi3); //function to calculate basis function value for given node A and xi1, xi2, xi3

        static std::tuple<double,double,double> basis_gradient(unsigned int node, double xi1, double xi2, double xi3); //function to calculate basis function gradient with respect to xi1, xi2, and xi3 for given node A and xi1, xi2, xi3
};