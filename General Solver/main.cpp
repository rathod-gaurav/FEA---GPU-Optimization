#include <iostream>
#include "MeshGenerator.hpp"

int main(){
    //Problem parameters
    // constexpr unsigned int Nsd = 3; //3D problem
    constexpr unsigned int Nne = 8; //hexahedral elements | 8 nodes per element

    //Domain parameters
    double x1_ll = 0.0, x1_ul = 0.1; //lower and upper limits in x1 direction
    double x2_ll = 0.0, x2_ul = 0.03; //lower and upper limits in x2 direction
    double x3_ll = 0.0, x3_ul = 0.03; //lower and upper limits in x3 direction
    
    //Mesh parameters
    unsigned int Nel_x1 = 80; //number of elements in x1 direction
    unsigned int Nel_x2 = 24; //number of elements in x2 direction
    unsigned int Nel_x3 = 24; //number of elements in x3 direction

    //Generate the mesh using the MeshGenerator class
    MeshGenerator<Nne> meshGen(x1_ll, x1_ul, x2_ll, x2_ul, x3_ll, x3_ul, Nel_x1, Nel_x2, Nel_x3);
    Mesh<Nne> mesh = meshGen.buildMesh();

    //Write the mesh into files ###### see if this can be made an asynchronous operation in the future to speed up the code
    mesh.writeToFiles("mesh"); //writes points.txt and elems.txt in the "mesh" directory

    std::cout << "Mesh built: " << mesh.Nnodes() << " nodes, " << mesh.Nelements() << " elements" << std::endl;
}