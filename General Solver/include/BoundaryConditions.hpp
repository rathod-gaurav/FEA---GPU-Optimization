#pragma once // include this only once during compilation

#include "Mesh.hpp"
#include <map>
#include <vector>
#include <Eigen/Dense>

template <unsigned int Nne>
class BoundaryConditions{
    public:
        BoundaryConditions(const Mesh<Nne>& mesh, unsigned int Nsd); //default constructor for BoundaryConditions class

        //function to store dirischlet boundary conditions map
        void addDirischlet(unsigned int node, int dof, double value); //function to add a dirischlet boundary condition for a given LOCAL node, degree of freedom, and value

        void buildBCs(); //function to build boundary conditions | must be called after all addDirichlet calls have been made to finalize the boundary conditions

        void applyToSolution(Eigen::VectorXd& u, double incrementFraction) const; //function to apply the boundary conditions to the solution vector u based on the current increment fraction (used for incremental loading)

        //Query methods
        const std::vector<unsigned int>& getUnknownIndexes() const {return unknownIndexes_;} //function to return the indexes of the unknown degrees of freedom | this is a const member function and does not modify the state of the object
        const std::vector<unsigned int>& getDirischletIndexes() const {return dirischletIndexes_;} //function to return the indexes of the dirischlet degrees of freedom | this is a const member function and does not modify the state of the object
        bool isDirischlet(unsigned int globalDOF) const {return isDirischlet_[globalDOF];} //function to check if a given global degree of freedom is subject to dirischlet boundary conditions | this is a const member function and does not modify the state of the object

        void printSummary() const; //function to print a summary of the boundary conditions | this is a const member function and does not modify the state of the object

    private:
        //these variables are private and can only be accessed within the BoundaryConditions class
        const Mesh<Nne>& mesh_; //the trailing underscore says - this belongs to the class and is not a local variable in the function | this is a common naming convention for class member variables
        unsigned int Nsd_; //number of spatial dimensions | this is also a member variable of the class
        unsigned int totalDOFs_; //total degrees of freedom in the system | this is calculated based on the number of nodes and spatial dimensions and is also a member variable of the class

        std::map<unsigned int, double> dirischletVals_; //map to store dirischlet values with key as global node index and value as the dirischlet value for that node
        std::vector<bool> isDirischlet_; //vector to indicate which degrees of freedom are subject to dirischlet boundary conditions
        std::vector<unsigned int> dirischletIndexes_; //indexes corresponding to dirischlet nodes and degrees of freedom
        std::vector<unsigned int> unknownIndexes_; //indexes corresponding to unknown nodes and degrees of freedom

};

#include "BoundaryConditions.tpp" //include the implementation of the template class in a separate .tpp file