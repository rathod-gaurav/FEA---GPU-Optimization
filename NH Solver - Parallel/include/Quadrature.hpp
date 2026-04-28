#pragma once //include this only once during compilation

#include <vector>

struct QuadratureRule {
    std::vector<double> points;
    std::vector<double> weights;
};

class Quadrature{
    public:
        static QuadratureRule gauss_legendre(unsigned int n);
};
