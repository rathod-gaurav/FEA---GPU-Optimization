#pragma once //include this file only once during compilation
#include <vector>
#include <string>

struct Node{
    double x1, x2, x3; //global coordinates of the node
};

template <unsigned int Nne>
struct Element{
    unsigned int node[Nne]; //IDs of the nodes that form the element
};


template <unsigned int Nne>
class Mesh{
    public:
        //Mesh class holds the nodes and elements lists
        std::vector<Node> nodes;
        std::vector<Element<Nne>> elements;
        
        //functions to return number of nodes and number of elements in the mesh
        unsigned int Nnodes() const { return nodes.size(); } //const here means this function does not modify the state of the object | read-only access to the mesh object
        unsigned int Nelements() const { return elements.size(); }

        void writeToFiles(const std::string& dir) const; //function to write the mesh into files | this is also a const member function and does not modify the state of the mesh object
};

#include "Mesh.tpp" //include the implementation of the member functions of the Mesh class