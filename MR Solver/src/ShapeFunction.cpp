#include "ShapeFunction.hpp"
#include <stdexcept> //for std::invalid_argument exception

std::tuple<double,double,double> ShapeFunction::xi_at_node(unsigned int node){ //function to return xi1, xi2, and xi3 for given node A
        double xi1, xi2, xi3;
        switch(node){
            case 0:
                xi1 = -1.0;
                xi2 = -1.0;
                xi3 = -1.0;
                break;
            case 1:
                xi1 = 1.0;
                xi2 = -1.0;
                xi3 = -1.0;
                break;
            case 2:
                xi1 = 1.0;
                xi2 = 1.0;
                xi3 = -1.0;
                break;
            case 3:
                xi1 = -1.0;
                xi2 = 1.0;
                xi3 = -1.0;
                break;
            case 4:
                xi1 = -1.0;
                xi2 = -1.0;
                xi3 = 1.0;
                break;
            case 5:
                xi1 = 1.0;
                xi2 = -1.0;
                xi3 = 1.0;
                break;
            case 6:
                xi1 = 1.0;
                xi2 = 1.0;
                xi3 = 1.0;
                break;
            case 7:
                xi1 = -1.0;
                xi2 = 1.0;
                xi3 = 1.0;
                break;
            default:
                throw std::invalid_argument("xi_at_node mapping not implemented for this local node number");
        }
        return {xi1, xi2, xi3};
};

double ShapeFunction::basis_function(unsigned int node, double xi1, double xi2, double xi3){
        auto [xi1_node , xi2_node , xi3_node] = xi_at_node(node);
        double value = 0.125*(1 + xi1*xi1_node)*(1 + xi2*xi2_node)*(1 + xi3*xi3_node);
        return value;
};

std::tuple<double,double,double> ShapeFunction::basis_gradient(unsigned int node, double xi1, double xi2, double xi3){
    auto [xi1_node,xi2_node,xi3_node] = xi_at_node(node);
    double basis_gradient_xi1, basis_gradient_xi2, basis_gradient_xi3;
    basis_gradient_xi1 = 0.125*xi1_node*(1 + xi2*xi2_node)*(1 + xi3*xi3_node);
    basis_gradient_xi2 = 0.125*xi2_node*(1 + xi1*xi1_node)*(1 + xi3*xi3_node);
    basis_gradient_xi3 = 0.125*xi3_node*(1 + xi1*xi1_node)*(1 + xi2*xi2_node);
    return {basis_gradient_xi1, basis_gradient_xi2, basis_gradient_xi3};
}

