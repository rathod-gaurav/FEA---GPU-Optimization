#pragma once //include this only once during compilation

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Mesh.hpp"
#include <ElementEvaluator.hpp>
#include <BoundaryConditions.hpp>
#include <unordered_set>

template <unsigned int Nne, unsigned int Nsd>
class Assembler{
    public: 
        Assembler(const Mesh<Nne>& mesh, const ElementEvaluator<Nne,Nsd>& elem_evaluator);

        void assembleSystem(
            const Eigen::VectorXd& u, //global nodal displacement vector (Nnodes*Nsd x 1 vector)
            Eigen::SparseMatrix<double>& Kglobal, //global stiffness matrix (Nnodes*Nsd x Nnodes*Nsd sparse matrix)
            Eigen::VectorXd& Rglobal //global internal force vector (Nnodes*Nsd x 1 vector)
        ) const;

        void partition(
            const Eigen::SparseMatrix<double>& Kglobal, //global stiffness matrix (Nnodes*Nsd x Nnodes*Nsd sparse matrix)
            Eigen::VectorXd& Rglobal, //global internal force vector (Nnodes*Nsd x 1 vector)
            const BoundaryConditions<Nne>& bcs, //boundary conditions object containing the indexes of the dirischlet DOFs

            Eigen::SparseMatrix<double>& KUU, //extract the submatrix of K corresponding to the unknown degrees of freedom
            Eigen::SparseMatrix<double>& KUD, //extract the submatrix of K corresponding to the coupling between unknown and dirischlet degrees of freedom
            Eigen::VectorXd& RU //extract the subvector of R corresponding to the unknown degrees of freedom
        ) const;
    
    private:
        Eigen::SparseMatrix<double> extractSparseSubmatrix(
            const Eigen::SparseMatrix<double>& K,
            const std::vector<unsigned int>& rows,
            const std::vector<unsigned int>& cols) const; //function to extract a sparse submatrix from the global stiffness matrix given row and column indexes
    
        const Mesh<Nne>& mesh_; //reference to the mesh object
        const ElementEvaluator<Nne,Nsd>& elem_evaluator_; //reference to the element evaluator object


};

#include "Assembler.tpp" //include the implementation of the Assembler class