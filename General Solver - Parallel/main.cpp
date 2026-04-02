#include <iostream>
#include "MeshGenerator.hpp"
#include "Quadrature.hpp"
#include "BoundaryConditions.hpp"
#include "ElementEvaluator.hpp"
#include "StVenantKirchhoff.hpp"
#include "Assembler.hpp"
#include "NonLinearSolver.hpp"
#include "OutputWriter.hpp"

#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration;

int main(){
    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;
    duration<double, std::milli> duration_msec;

    std::cout << "Starting the timer..." << std::endl;
    start = high_resolution_clock::now();

    //Problem parameters
    constexpr unsigned int Nsd = 3; //3D problem
    constexpr unsigned int Nne = 8; //hexahedral elements | 8 nodes per element

    //Quadrature order
    unsigned int quadOrder = 3; //number of quadrature points in each direction for

    //Problem parameters
    double lambda = 6*1e10; //first Lamé parameter
    double mu = 2*1e10; //second Lamé parameter (shear modulus)

    //Solver parameters
    double tol = 1e-6; //tolerance for convergence of the nonlinear solver
    unsigned int maxIncr = 10; //maximum number of increments (timesteps)
    unsigned int maxIter = 20; //maximum number of iterations per increment

    //Domain parameters
    double x1_ll = 0.0, x1_ul = 0.1; //lower and upper limits in x1 direction
    double x2_ll = 0.0, x2_ul = 0.03; //lower and upper limits in x2 direction
    double x3_ll = 0.0, x3_ul = 0.03; //lower and upper limits in x3 direction
    
    //Mesh parameters
    unsigned int Nel_x1 = 10; //number of elements in x1 direction
    unsigned int Nel_x2 = 3; //number of elements in x2 direction
    unsigned int Nel_x3 = 3; //number of elements in x3 direction

    //Generate the mesh using the MeshGenerator class
    MeshGenerator<Nne> meshGen(x1_ll, x1_ul, x2_ll, x2_ul, x3_ll, x3_ul, Nel_x1, Nel_x2, Nel_x3);
    Mesh<Nne> mesh = meshGen.buildMesh();

    //Write the mesh into files ###### see if this can be made an asynchronous operation in the future to speed up the code
    mesh.writeToFiles("mesh"); //writes points.txt and elems.txt in the "mesh" directory

    std::cout << "Mesh built: " << mesh.Nnodes() << " nodes, " << mesh.Nelements() << " elements" << std::endl;
    std::cout << "--------------------" << std::endl;

    //Boundary conditions
    BoundaryConditions<Nne> bcs(mesh, Nsd);
    for(unsigned int i = 0 ; i < mesh.Nnodes() ; i++){
        if(mesh.nodes[i].x1 == x1_ll){ //if the node is on the left face of the domain
            bcs.addDirischlet(i, 0, 0.0); //apply dirischlet boundary condition u1 = 0 at this node
            bcs.addDirischlet(i, 1, 0.0); //apply dirischlet boundary condition u1 = 0 at this node
            bcs.addDirischlet(i, 2, 0.0); //apply dirischlet boundary condition u1 = 0 at this node
        }
        if(mesh.nodes[i].x1 == x1_ul){ //if the node is on the right face of the domain
            bcs.addDirischlet(i, 0, 0.01); //apply dirischlet boundary condition u1 = 0.001 at this node
        }
    }
    bcs.buildBCs(); //finalize the boundary conditions
    bcs.printSummary(); //print a summary of the boundary conditions
    std::cout << "--------------------" << std::endl;


    //Problem physics stack
    QuadratureRule              quadRule = Quadrature::gauss_legendre(quadOrder); //get the quadrature points and weights for the specified quadrature order
    StVenantKirchhoff           material(lambda, mu); //create an instance of the St. Venant-Kirchhoff material model with the specified Lamé parameters
    ElementEvaluator<Nne, Nsd>  elemEval(mesh, material, quadRule); //create an instance of element evaluator with the mesh, material model, and quadrature rule
    Assembler<Nne, Nsd>         assembler(mesh, elemEval); //create an instance of the assembler with the mesh and element evaluator
    OutputWriter<Nne>           writer("solutions"); //create an instance of the output writer to write results to the "output" directory

    NonlinearSolver<Nne, Nsd>   solver(tol, maxIncr, maxIter); //create an instance of the nonlinear solver with a tolerance of 1e-6, maximum 10 increments, and maximum 20 iterations per increment

    Eigen::VectorXd u = Eigen::VectorXd::Zero(mesh.Nnodes()*Nsd); //initialize the global displacement vector to zero

    std::cout << "Starting nonlinear solve..." << std::endl;
    std::cout << "--------------------" << std::endl;
    solver.solve(u, assembler, bcs,
        [&](unsigned int incr, unsigned int iter, double residualNorm){
            writer.sendResidual(incr, iter, residualNorm); //send the residual norm to the output writer for visualization
            
            if(iter == 0){ //write the solution at the first iteration of each increment
                writer.writeVTU(mesh, u, incr); //write the current solution vector and mesh information to files for visualization
            }
        }
    ); //solve the nonlinear system to get the nodal displacements
    std::cout << "--------------------" << std::endl;
    std::cout << "Nonlinear solve completed." << std::endl;

    end = high_resolution_clock::now();
    duration_msec = std::chrono::duration_cast<duration<double, std::milli>>(end-start);

    std::cout << "Total time taken (in ms): " << duration_msec.count() << std::endl;
}