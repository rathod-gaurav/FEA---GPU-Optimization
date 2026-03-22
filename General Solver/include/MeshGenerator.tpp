#include "MeshGenerator.hpp"

template <unsigned int Nne>
MeshGenerator<Nne>::MeshGenerator( //assign the function parameters to the class member variables using an initializer list
    //domain parameters
    double x1_ll, double x1_ul,
    double x2_ll, double x2_ul,
    double x3_ll, double x3_ul,
    //mesh parameters
    unsigned int Nel_x1, unsigned int Nel_x2, unsigned int Nel_x3
):
    x1_ll_(x1_ll), x1_ul_(x1_ul),
    x2_ll_(x2_ll), x2_ul_(x2_ul),
    x3_ll_(x3_ll), x3_ul_(x3_ul),
    Nel_x1_(Nel_x1), Nel_x2_(Nel_x2), Nel_x3_(Nel_x3)
{}

template <unsigned int Nne>
Mesh<Nne> MeshGenerator<Nne>::buildMesh() const{
    Mesh<Nne> mesh; //create an empty mesh object

    unsigned int Nnodes_x1 = Nel_x1_ + 1; //number of nodes in x1 direction
    unsigned int Nnodes_x2 = Nel_x2_ + 1; //number of nodes in x2 direction
    unsigned int Nnodes_x3 = Nel_x3_ + 1; //number of nodes in x3 direction
    unsigned int Nt = Nnodes_x1 * Nnodes_x2 * Nnodes_x3; //total number of nodes

    double dx1 = (x1_ul_ - x1_ll_) / Nel_x1_; //spacing between nodes in x1 direction
    double dx2 = (x2_ul_ - x2_ll_) / Nel_x2_; //spacing between nodes in x2 direction
    double dx3 = (x3_ul_ - x3_ll_) / Nel_x3_; //spacing between nodes in x3 direction

    //Build the nodes list of the mesh
    mesh.nodes.reserve(Nt);
    for(unsigned int k = 0 ; k < Nnodes_x3 ; k++){
        for(unsigned int j = 0 ; j < Nnodes_x2 ; j++){
            for(unsigned int i = 0 ; i < Nnodes_x1 ; i++){
                Node n;
                n.x1 = x1_ll_ + i*dx1;
                n.x2 = x2_ll_ + j*dx2;
                n.x3 = x3_ll_ + k*dx3;
                mesh.nodes.push_back(n);
            }
        }
    }

    //variables required for element connectivity
    unsigned int Nel_t = Nel_x1_ * Nel_x2_ * Nel_x3_; //total number of elements

    //Build the elements list of the mesh
    mesh.elements.reserve(Nel_t);
    for(unsigned int k = 0 ; k < Nel_x3_ ; k++){
        for(unsigned int j = 0 ; j < Nel_x2_ ; j++){
            for(unsigned int i = 0 ; i < Nel_x1_ ; i++){
                Element<Nne> elem;
                unsigned int base = i 
                     + j * Nnodes_x1 
                     + k * (Nnodes_x1 * Nnodes_x2);

                unsigned int n0 = base;
                unsigned int n1 = base + 1;
                unsigned int n3 = base + Nnodes_x1;
                unsigned int n2 = n3 + 1;

                unsigned int n4 = base + Nnodes_x1 * Nnodes_x2;
                unsigned int n5 = n4 + 1;
                unsigned int n7 = n4 + Nnodes_x1;
                unsigned int n6 = n7 + 1;

                elem.node[0] = n0;
                elem.node[1] = n1;
                elem.node[2] = n2;
                elem.node[3] = n3;
                elem.node[4] = n4;
                elem.node[5] = n5;
                elem.node[6] = n6;
                elem.node[7] = n7;

                mesh.elements.push_back(elem);
            }
        }
    }

    return mesh; //return the built mesh object

}