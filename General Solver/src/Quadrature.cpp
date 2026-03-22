#include "Quadrature.hpp"
#include <stdexcept> // for std::invalid_argument

QuadratureRule Quadrature::gauss_legendre(unsigned int n) {
    QuadratureRule rule;

    switch(n) {
        case 1:
            rule.points  = { 0.0 };
            rule.weights = { 2.0 };
            break;

        case 2:
            rule.points  = { -0.5773502691896257,  0.5773502691896257 };
            rule.weights = {  1.0,                 1.0 };
            break;

        case 3:
            rule.points  = { -0.7745966692414834, 0.0, 0.7745966692414834 };
            rule.weights = {  0.5555555555555556, 0.8888888888888888, 0.5555555555555556 };
            break;

        case 4:
            rule.points  = { -0.8611363115940526, -0.3399810435848563,
                              0.3399810435848563,  0.8611363115940526 };
            rule.weights = {  0.3478548451374539,  0.6521451548625461,
                              0.6521451548625461,  0.3478548451374539 };
            break;

        case 5:
            rule.points  = { -0.9061798459386640, -0.5384693101056831,
                              0.0,
                              0.5384693101056831,  0.9061798459386640 };
            rule.weights = {  0.2369268850561891,  0.4786286704993665,
                              0.5688888888888889,  0.4786286704993665,
                              0.2369268850561891 };
            break;

        case 6:
            rule.points  = { -0.9324695142031521, -0.6612093864662645,
                             -0.2386191860831969,  0.2386191860831969,
                              0.6612093864662645,  0.9324695142031521 };
            rule.weights = {  0.1713244923791704,  0.3607615730481386,
                              0.4679139345726910,  0.4679139345726910,
                              0.3607615730481386,  0.1713244923791704 };
            break;

        default:
            throw std::invalid_argument("Gauss-Legendre quadrature not implemented for this n");
    }

    return rule;
}