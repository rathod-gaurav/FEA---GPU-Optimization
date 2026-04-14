#pragma once //include this only once during compilation

#include <Eigen/Dense>
#include "Mesh.hpp"
#include "ShapeFunction.hpp"
#include "Quadrature.hpp"
#include "MaterialModel.hpp"

template <unsigned int Nne, unsigned int Nsd>
class ElementEvaluator{
    public:
        ElementEvaluator( //default constructor
            const Mesh<Nne>& mesh,
            const MaterialModel& material,
            const QuadratureRule& quadRule
        );

        void computeElement(
            unsigned int e, //element index
            const Eigen::VectorXd& u_e, //element nodal displacements (Nne*3 x 1 vector)
            Eigen::MatrixXd& Klocal, //element stiffness matrix (Nne*3 x Nne*3 matrix)
            Eigen::VectorXd& Rlocal //element internal force vector (Nne*3 x 1 vector)
        ) const;
    
    private:
        Eigen::Matrix3d computeJacobian(unsigned int e, double xi1, double xi2, double xi3) const; //function to compute the Jacobian matrix for the element at given quadrature point (xi1, xi2, xi3)

        Eigen::Matrix3d computeGradU(const Eigen::VectorXd& u_e, double xi1, double xi2, double xi3, Eigen::Matrix3d& JacInv) const; //function to compute the gradient of the displacement field at the quadrature point using the basis function gradients and the nodal displacements

        const Mesh<Nne>& mesh_; //reference to the mesh object
        const MaterialModel& material_; //reference to the material model object
        const QuadratureRule& quadRule_; //reference to the quadrature rule object
};

#include "ElementEvaluator.tpp" //include the implementation of the template class