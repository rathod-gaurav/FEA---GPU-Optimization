// nonlinear FE Problem for finite strain mechanics using St. Venant-Kirchhoff strain energy density function
// only Dirischlet boundary conditions specified at left and right faces. No Neumann BC, No body force.

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <Eigen/Dense>
#include <map>
#include <chrono>

using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration;

struct Node {
    double x1, x2, x3; //global coordinates of the node
};

template <unsigned int Nne>
struct Element {
    unsigned int node[Nne]; //IDs of the nodes that form the element
};

std::tuple<double,double,double> xi_at_node(unsigned int node){ //function to return xi1, xi2, and xi3 for given node A
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

double basis_function(unsigned int node, double xi1, double xi2, double xi3){
        auto [xi1_node , xi2_node , xi3_node] = xi_at_node(node);
        double value = 0.125*(1 + xi1*xi1_node)*(1 + xi2*xi2_node)*(1 + xi3*xi3_node);
        return value;
};

std::tuple<double,double,double> basis_gradient(unsigned int node, double xi1, double xi2, double xi3){
    auto [xi1_node,xi2_node,xi3_node] = xi_at_node(node);
    double basis_gradient_xi1, basis_gradient_xi2, basis_gradient_xi3;
    basis_gradient_xi1 = 0.125*xi1_node*(1 + xi2*xi2_node)*(1 + xi3*xi3_node);
    basis_gradient_xi2 = 0.125*xi2_node*(1 + xi1*xi1_node)*(1 + xi3*xi3_node);
    basis_gradient_xi3 = 0.125*xi3_node*(1 + xi1*xi1_node)*(1 + xi2*xi2_node);
    return {basis_gradient_xi1, basis_gradient_xi2, basis_gradient_xi3};
}

struct QuadratureRule {
    std::vector<double> points;
    std::vector<double> weights;
};

QuadratureRule gauss_legendre(unsigned int n) {
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

template <unsigned int Nne>
Eigen::MatrixXd compute_grad_u(Eigen::VectorXd u_e, double xi1, double xi2, double xi3, Eigen::MatrixXd JacInv){
    Eigen::MatrixXd grad_u = Eigen::MatrixXd::Zero(3,3);
    //compute the gradient of the displacement field at the quadrature point using the basis function gradients and the nodal displacements
    for(int A = 0 ; A < Nne ; A++){
        auto [dN_dxi1, dN_dxi2, dN_dxi3] = basis_gradient(A, xi1, xi2, xi3);
        Eigen::VectorXd dN_dx = JacInv.transpose()*Eigen::Vector3d(dN_dxi1, dN_dxi2, dN_dxi3);
        grad_u(0,0) += dN_dx[0] * u_e(A*3 + 0); //du1/dx1
        grad_u(0,1) += dN_dx[1] * u_e(A*3 + 0); //du1/dx2
        grad_u(0,2) += dN_dx[2] * u_e(A*3 + 0); //du1/dx3

        grad_u(1,0) += dN_dx[0] * u_e(A*3 + 1); //du2/dx1
        grad_u(1,1) += dN_dx[1] * u_e(A*3 + 1); //du2/dx2
        grad_u(1,2) += dN_dx[2] * u_e(A*3 + 1); //du2/dx3

        grad_u(2,0) += dN_dx[0] * u_e(A*3 + 2); //du3/dx1
        grad_u(2,1) += dN_dx[1] * u_e(A*3 + 2); //du3/dx2
        grad_u(2,2) += dN_dx[2] * u_e(A*3 + 2); //du3/dx3
    }
    return grad_u;
}

Eigen::MatrixXd extractSubmatrix(const Eigen::MatrixXd& OriginalMatrix , const vector<int> rows , const vector<int> cols){//extract submatrix from original matrix given row and column indexes
    Eigen::MatrixXd subMatrix(rows.size(), cols.size());

    for(int i = 0 ; i < rows.size() ; i++){
        for(int j = 0 ; j < cols.size() ; j++){
            subMatrix(i,j) = OriginalMatrix(rows[i],cols[j]);
        }
    }
    return subMatrix;
}

#include <sstream>
#include <curl/curl.h>
void sendResidual(int incr, int iter, double residual)
{
    CURL *curl = curl_easy_init();

    if(curl)
    {
        std::stringstream json;

        json << "{";
        json << "\"increment\": " << incr << ",";
        json << "\"iteration\": " << iter << ",";
        json << "\"residual\": " << residual;
        json << "}";

        std::string jsonStr = json.str();

        struct curl_slist *headers = NULL;
        headers = curl_slist_append(headers, "Content-Type: application/json");

        curl_easy_setopt(curl, CURLOPT_URL, "http://127.0.0.1:8000/residual");
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, jsonStr.c_str());
        curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

        curl_easy_perform(curl);

        curl_easy_cleanup(curl);
    }
}

template <unsigned int Nne>
void write_vtu(
    const std::string& filename,
    const std::vector<Node>& nodes,
    const std::vector<Element<Nne>>& elements,
    const Eigen::VectorXd& displacement   // size = 3*Nnodes
)
{
    std::ofstream file(filename);

    int Nnodes = nodes.size();
    int Nelems = elements.size();

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "<UnstructuredGrid>\n";
    file << "<Piece NumberOfPoints=\"" << Nnodes 
         << "\" NumberOfCells=\"" << Nelems << "\">\n";

    // ---------------------
    // POINTS
    // ---------------------
    file << "<Points>\n";
    file << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (int i = 0; i < Nnodes; ++i)
        file << nodes[i].x1 << " "
             << nodes[i].x2 << " "
             << nodes[i].x3 << "\n";

    file << "</DataArray>\n</Points>\n";

    // ---------------------
    // CELLS (Hex = VTK type 12)
    // ---------------------
    file << "<Cells>\n";

    // connectivity
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int e = 0; e < Nelems; ++e)
        for (int A = 0; A < 8; ++A)
            file << elements[e].node[A] << " ";
    file << "\n</DataArray>\n";

    // offsets
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int e = 0; e < Nelems; ++e)
        file << (e+1)*8 << " ";
    file << "\n</DataArray>\n";

    // types (12 = hexahedron)
    file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int e = 0; e < Nelems; ++e)
        file << 12 << " ";
    file << "\n</DataArray>\n";

    file << "</Cells>\n";

    // ---------------------
    // POINT DATA (Displacement)
    // ---------------------
    file << "<PointData Vectors=\"Displacement\">\n";
    file << "<DataArray type=\"Float32\" Name=\"Displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (int i = 0; i < Nnodes; ++i)
    {
        file << displacement(3*i)   << " "
             << displacement(3*i+1) << " "
             << displacement(3*i+2) << "\n";
    }

    file << "</DataArray>\n</PointData>\n";

    file << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

void write_pvd(
    const std::string& filename,
    const std::vector<std::string>& vtu_files,
    const std::vector<double>& times)
{
    std::ofstream file(filename);

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "<Collection>\n";

    for (size_t i = 0; i < vtu_files.size(); ++i)
    {
        file << "<DataSet timestep=\"" << times[i]
             << "\" group=\"\" part=\"0\" file=\""
             << vtu_files[i] << "\"/>\n";
    }

    file << "</Collection>\n</VTKFile>\n";
}

int main(){
    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;
    duration<double, std::milli> duration_msec;

    cout << "Starting the timer..." << endl;
    start = high_resolution_clock::now();

    unsigned int Nsd = 3; //number of spatial dimensions - 3D problem
    constexpr int Nne = 8; //number of nodes per element - 8-node hexahedral element
    unsigned int quadRule = 3; //number of quadrature points in each direction for numerical integration
    double epsilon = 1e-6; //Newton Raphson solver tolerance

    //Setup Newton-Raphson increment and iteration parameters
    unsigned int maxIncrement = 10; //maximum number of load increments
    unsigned int maxIter = 100; //maximum number of iterations per increment

    //problem variables
    double lambda = 6*1e10; //first Lamé parameter
    double mu = 2*1e10; //second Lamé parameter (shear modulus)

    //domain
    double x1_ll = 0.0;
    double x1_ul = 0.1;
    double x2_ll = 0.0;
    double x2_ul = 0.03;
    double x3_ll = 0.0;
    double x3_ul = 0.03;

    //Mesh
    unsigned int Nel_x1 = 40; //number of elements in x1 direction
    unsigned int Nel_x2 = 12; //number of elements in x2 direction
    unsigned int Nel_x3 = 12; //number of elements in x3 direction

    unsigned int Nnodes_x1 = Nel_x1 + 1; //number of nodes in x1 direction
    unsigned int Nnodes_x2 = Nel_x2 + 1; //number of nodes in x2 direction
    unsigned int Nnodes_x3 = Nel_x3 + 1; //number of nodes in x3 direction

    double dx1 = (x1_ul - x1_ll) / Nel_x1; //element size in x1 direction
    double dx2 = (x2_ul - x2_ll) / Nel_x2; //element size in x2 direction
    double dx3 = (x3_ul - x3_ll) / Nel_x3; //element size in x3 direction

    unsigned int Nel_t = Nel_x1 * Nel_x2 * Nel_x3; //total number of elements
    unsigned int Nt = Nnodes_x1 * Nnodes_x2 * Nnodes_x3; //total number of nodes

    vector<Node> nodes;
    nodes.reserve(Nt);
    for(unsigned int k = 0 ; k < Nnodes_x3 ; k++){
        for(unsigned int j = 0 ; j < Nnodes_x2 ; j++){
            for(unsigned int i = 0 ; i < Nnodes_x1 ; i++){
                Node n;
                n.x1 = x1_ll + i*dx1;
                n.x2 = x2_ll + j*dx2;
                n.x3 = x3_ll + k*dx3;
                nodes.push_back(n);
            }
        }
    }

    //Local-Global node number mapping for every element
    using Element3D = Element<Nne>;
    vector<Element3D> elements;
    elements.reserve(Nel_t);
    for(unsigned int k = 0 ; k < Nel_x3 ; k++){
        for(unsigned int j = 0 ; j < Nel_x2 ; j++){
            for(unsigned int i = 0 ; i < Nel_x1 ; i++){
                Element3D elem;
                // int n0 = i + j*Nnodes_x1 + k*(Nnodes_x1*Nnodes_x2);
                // int n1 = n0 + 1;
                // int n2 = n1 + (Nnodes_x1*Nnodes_x2);
                // int n3 = n2 - 1;
                // int n4 = i + (j+1)*Nnodes_x1 + k*(Nnodes_x1*Nnodes_x2);
                // int n5 = n4 + 1;
                // int n6 = n5 + (Nnodes_x1*Nnodes_x2);
                // int n7 = n6 - 1;

                int base = i 
                     + j * Nnodes_x1 
                     + k * (Nnodes_x1 * Nnodes_x2);

                int n0 = base;
                int n1 = base + 1;
                int n3 = base + Nnodes_x1;
                int n2 = n3 + 1;

                int n4 = base + Nnodes_x1 * Nnodes_x2;
                int n5 = n4 + 1;
                int n7 = n4 + Nnodes_x1;
                int n6 = n7 + 1;

                elem.node[0] = n0;
                elem.node[1] = n1;
                elem.node[2] = n2;
                elem.node[3] = n3;
                elem.node[4] = n4;
                elem.node[5] = n5;
                elem.node[6] = n6;
                elem.node[7] = n7;

                elements.push_back(elem);
            }
        }
    }

    //store the mesh into points and hexa files
    std::ofstream points_file("mesh/points.txt");
    for(auto& node : nodes){
        points_file << node.x1 << " " << node.x2  << " " << node.x3 << "\n";
    }

    std::ofstream hexas_file("mesh/hexas.txt");
    for(auto& elem : elements){
        hexas_file << elem.node[0] << " " << elem.node[1] << " " << elem.node[2] << " " << elem.node[3] << " " << elem.node[4] << " " << elem.node[5] << " " << elem.node[6] << " " << elem.node[7] << "\n";
    }

    // Initialize the solution vector (displacements at each node)
    Eigen::VectorXd u = Eigen::VectorXd::Zero(Nt * Nsd); //displacement vector initialized to zero
    Eigen::VectorXd du = Eigen::VectorXd::Zero(Nt * Nsd); //incremental displacement vector initialized to zero


    auto calculate_Jacobian_3D = [Nsd, elements, nodes](int e, double xi1, double xi2, double xi3){//function to calculate jacobian
        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(Nsd,Nsd);
        
        for(int A = 0 ; A < Nne ; A++){
            auto [basis_gradient_xi1, basis_gradient_xi2, basis_gradient_xi3] = basis_gradient(A, xi1, xi2, xi3);
            int Aglobal = elements[e].node[A];
            J(0,0) += basis_gradient_xi1*nodes[Aglobal].x1; //dx1/dxi1
            J(0,1) += basis_gradient_xi2*nodes[Aglobal].x1; //dx1/dxi2
            J(0,2) += basis_gradient_xi3*nodes[Aglobal].x1; //dx1/dxi3
            J(1,0) += basis_gradient_xi1*nodes[Aglobal].x2; //dx2/dxi1
            J(1,1) += basis_gradient_xi2*nodes[Aglobal].x2; //dx2/dxi2
            J(1,2) += basis_gradient_xi3*nodes[Aglobal].x2; //dx2/dxi3
            J(2,0) += basis_gradient_xi1*nodes[Aglobal].x3; //dx3/dxi1
            J(2,1) += basis_gradient_xi2*nodes[Aglobal].x3; //dx3/dxi2
            J(2,2) += basis_gradient_xi3*nodes[Aglobal].x3; //dx3/dxi3
        }
        return J;
    };

    //Quadrature points
    QuadratureRule q = gauss_legendre(quadRule);
    std::vector<double> points(quadRule), weights(quadRule);
    points = q.points;
    weights = q.weights;
    Eigen::VectorXd quad_points = Eigen::Map<Eigen::VectorXd>(points.data(), points.size());
    Eigen::VectorXd quad_weights = Eigen::Map<Eigen::VectorXd>(weights.data(), weights.size());

    //Boundary Conditions
    //global nodeLocations where dirischlet boundary conditions are specified
    std::map<int, std::vector<int>> nodeLocationsD_map;
    for(int i = 0 ; i < Nt ; i++){
        if(nodes[i].x1 == x1_ll){
            nodeLocationsD_map[i].push_back(0); // 0 => X1 displacement specified on this node
            nodeLocationsD_map[i].push_back(1); // 1 => X2 displacement specified on this node
            nodeLocationsD_map[i].push_back(2); // 2 => X3 displacement specified on this node
        }
        if(nodes[i].x1 == x1_ul){
            nodeLocationsD_map[i].push_back(0); // 0 => X1 displacement specified on this node
        }
    }
    vector<bool> isDirischlet(Nt,false);
    for(const auto& [key,vec] : nodeLocationsD_map){
        isDirischlet[key] = true;
    }

    //indexes to remove from solution array corresponding to dirischlet nodes and degrees of freedom
    vector<int> dirischletIndexes;
    for(int i = 0 ; i < Nt ; i++){
        if(isDirischlet[i]){
            for(int dof : nodeLocationsD_map[i]){
                dirischletIndexes.push_back(i*Nsd + dof);
            }
        }
    }
    vector<bool> isDirischletIndex(Nt*Nsd,false);
    for(int index : dirischletIndexes){
        isDirischletIndex[index] = true;
    }
    vector<int> unknownIndexes;
    for(int i = 0 ; i < Nt*Nsd ; i++){
        if(!isDirischletIndex[i]){
            unknownIndexes.push_back(i);
        }
    }

    cout << "Number of unknowns: " << unknownIndexes.size() << "\n";
    cout << "Number of dirischlet conditions: " << dirischletIndexes.size() << "\n";
    cout << "Total number of DOFs: " << Nt*Nsd << "\n";

    //given values of displacement field at dirischlet boundary
    Eigen::VectorXd dirischletVal(dirischletIndexes.size());
    for(int i = 0 ; i < dirischletIndexes.size() ; i++){
        int nodeD = dirischletIndexes[i]/3;
        int dof = dirischletIndexes[i]%3;
        if(nodes[nodeD].x1 == x1_ll){
            dirischletVal(i) = 0.0; //displacements at x1 = 0 are fixed to zero
        }
        else if(nodes[nodeD].x1 == x1_ul){
            if(dof == 0){
                dirischletVal(i) = 0.01; //X1 displacements at x1 = 10 are fixed to 0.01
            }
        }
    }
    Eigen::VectorXd dirischletValDot = Eigen::VectorXd::Zero(dirischletIndexes.size()); //time derivative of dirischlet values for dynamic problems, initialized to zero for static problem

    //solution files and increments (refered to here as timesteps)
    std::vector<std::string> solution_files;
    std::vector<double> solution_timesteps;

    cout << "Starting Newton-Raphson Iterations..." << "\n";
    for(unsigned int increment = 0; increment < maxIncrement; increment++){
        //Apply dirischlet BCs
        for(int i = 0 ; i < dirischletIndexes.size() ; i++){
            u(dirischletIndexes[i]) = ((increment+1)/static_cast<double>(maxIncrement))*dirischletVal(i);
        }
        // cout << "Applied Dirischlet BCs for increment " << increment+1 << "\n";

        for(unsigned int iter = 0; iter < maxIter; iter++){
            Eigen::VectorXd Rglobal = Eigen::VectorXd::Zero(Nt * Nsd); //residual vector initialized to zero
            Eigen::MatrixXd Kglobal = Eigen::MatrixXd::Zero(Nt * Nsd, Nt * Nsd); //tangent stiffness matrix initialized to zero
            cout << "Initialized Kglobal and Rglobal to zero for increment " << increment+1 << ", iteration " << iter+1 << "\n";
            //Loop over elements to compute element-level contributions to R and K
            for(unsigned int e = 0; e < Nel_t; e++){
                //Get the nodes of the current element`
                Eigen::VectorXd Rlocal = Eigen::VectorXd::Zero(Nne * Nsd); //local residual vector for the element
                Eigen::MatrixXd Klocal = Eigen::MatrixXd::Zero(Nne * Nsd, Nne * Nsd); //local tangent stiffness matrix for the element
                
                //Element nodal displacements
                Eigen::VectorXd u_e = Eigen::VectorXd::Zero(Nne * Nsd); //displacement vector for the current element
                for(unsigned int i = 0; i < Nne; i++){
                    unsigned int global_node_id = elements[e].node[i];
                    u_e.segment(i*Nsd, Nsd) = u.segment(global_node_id*Nsd, Nsd); //extract the displacements for the nodes of the current element 
                }

                //Gaussian quadrature loop
                for(int I = 0 ; I < quadRule ; I++){
                    for(int J = 0 ; J < quadRule ; J++){
                        for(int K = 0 ; K < quadRule ; K++){
                            //Get the quadrature point coordinates and weights
                            double xi1 = quad_points[I]; 
                            double xi2 = quad_points[J];
                            double xi3 = quad_points[K];
                            double weight = quad_weights[I] * quad_weights[J] * quad_weights[K];
                            Eigen::MatrixXd Jac = calculate_Jacobian_3D(e, xi1, xi2, xi3); //compute the Jacobian matrix at the quadrature point
                            double JacDet = Jac.determinant(); //compute the determinant of the Jacobian
                            Eigen::MatrixXd JacInv = Jac.inverse(); //compute the inverse of the Jacobian

                            Eigen::MatrixXd grad_u = compute_grad_u<Nne>(u_e, xi1, xi2, xi3, JacInv); //compute the gradient of the displacement field at the quadrature point
                            //Compute the deformation gradient F, Green-Lagrange strain E, and the second Piola-Kirchhoff stress S at the quadrature point
                            Eigen::MatrixXd F = Eigen::MatrixXd::Identity(3,3) + grad_u; //deformation gradient
                            Eigen::MatrixXd E = 0.5 * (F.transpose() * F - Eigen::MatrixXd::Identity(3,3)); //Green-Lagrange strain
                            Eigen::MatrixXd S = 2*mu*E + lambda*E.trace()*Eigen::MatrixXd::Identity(3,3); //second Piola-Kirchhoff stress using St. Venant-Kirchhoff model
                            Eigen::MatrixXd P = F * S; //first Piola-Kirchhoff stress

                            for(int B = 0 ; B < Nne ; B++){//Loop to calculate Residual
                                auto [dN_dxi1, dN_dxi2, dN_dxi3] = basis_gradient(B, xi1, xi2, xi3);
                                Eigen::VectorXd dN_dx = JacInv.transpose()*Eigen::Vector3d(dN_dxi1, dN_dxi2, dN_dxi3); //gradient of the basis function in global coordinates
                                Rlocal.segment(B*Nsd, Nsd) += P * dN_dx * weight * JacDet; //contribution to the local residual vector
                            }
                            
                            cout << "Calculated Rlocal for element " << e+1 << "/" << Nel_t << "\r";
                            
                            for(int A = 0 ; A < Nne ; A++){//Loops to calculate tangent matrix
                                for(int B = 0 ; B < Nne ; B++){

                                    auto [dNA_dxi1, dNA_dxi2, dNA_dxi3] = basis_gradient(A, xi1, xi2, xi3);
                                    auto [dNB_dxi1, dNB_dxi2, dNB_dxi3] = basis_gradient(B, xi1, xi2, xi3);

                                    Eigen::VectorXd dNA_dx = JacInv.transpose()*Eigen::Vector3d(dNA_dxi1, dNA_dxi2, dNA_dxi3);
                                    Eigen::VectorXd dNB_dx = JacInv.transpose()*Eigen::Vector3d(dNB_dxi1, dNB_dxi2, dNB_dxi3);

                                    //Kgeometric
                                    double Kgeo_scalar = (dNA_dx.transpose() * S * dNB_dx)(0,0);
                                    Eigen::MatrixXd KgeoAB =  Kgeo_scalar * JacDet * weight * Eigen::MatrixXd::Identity(3,3);
                                    Klocal.block<3,3>(3*A,3*B) += KgeoAB;

                                    //Kmaterial
                                    // Eigen::MatrixXd FA = F.transpose()*dNA_dx;
                                    // Eigen::MatrixXd FB = F.transpose()*dNB_dx;

                                    // Eigen::MatrixXd KmatAB = ((lambda + mu)*FA*FB.transpose() + mu*(F * dNA_dx.cwiseProduct(dNB_dx).asDiagonal() * F.transpose())) * JacDet * weight;
                                    // Klocal.block<3,3>(3*A,3*B) += KmatAB;


                                    // Correct KmatAB (3x3 block for nodes A,B)
                                    Eigen::MatrixXd KmatAB = Eigen::MatrixXd::Zero(3,3);
                                    for(int i = 0; i < 3; i++){
                                        for(int j = 0; j < 3; j++){
                                            double val = 0.0;
                                            for(int P = 0; P < 3; P++){
                                                for(int Q = 0; Q < 3; Q++){
                                                    for(int M = 0; M < 3; M++){
                                                        for(int N = 0; N < 3; N++){
                                                            double C_PQMN = lambda*(P==Q ? 1:0)*(M==N ? 1:0) + mu*((P==M ? 1:0)*(Q==N ? 1:0) + (P==N ? 1:0)*(Q==M ? 1:0));
                                                            val += F(i,P)*C_PQMN*F(j,M)*dNA_dx(Q)*dNB_dx(N);
                                                        }
                                                    }
                                                }
                                            }
                                            KmatAB(i,j) = val;
                                        }
                                    }
                                    KmatAB *= JacDet * weight;
                                    Klocal.block<3,3>(3*A,3*B) += KmatAB;

                                }
                            }
                            
                            cout << "Calculated Klocal for element " << e+1 << "/" << Nel_t << "\r";
                            
                        }
                    }
                }

                

                //Assemble Rlocal and Klocal into Rglobal and Kglobal
                for(int A = 0; A < Nne; A++){
                    int Aglobal = elements[e].node[A];
                    for(int B = 0; B < Nne ; B++)
                    {
                        int Bglobal = elements[e].node[B];
                        Kglobal.block<3,3>(3*Aglobal,3*Bglobal) += Klocal.block<3,3>(3*A,3*B);
                    }
                    Rglobal.segment(3*Aglobal,3) += Rlocal.segment(3*A,3);
                }

                cout << "Assembled element " << e+1 << "/" << Nel_t << "\r";
            }

            Eigen::MatrixXd KUU = extractSubmatrix(Kglobal, unknownIndexes, unknownIndexes); //extract the submatrix of K corresponding to the unknown degrees of freedom
            Eigen::MatrixXd KUD = extractSubmatrix(Kglobal, unknownIndexes, dirischletIndexes); //extract the submatrix of K corresponding to the coupling between unknown and dirischlet degrees of freedom
            Eigen::VectorXd RU(unknownIndexes.size()); //extract the subvector of R corresponding to the unknown degrees of freedom
            for(int i = 0; i < unknownIndexes.size() ; i++){
                RU(i) = Rglobal(unknownIndexes[i]);
            }
            Eigen::VectorXd R(RU.size()); //final residual after applying dirischlet boundary conditions
            R = RU; //modify the residual to account for the known displacements at the dirischlet nodes

            cout << "Initilised solver for increment " << increment+1 << ", iteration " << iter+1 << "\n";

            Eigen::LDLT<Eigen::MatrixXd> solver(KUU);
            Eigen::VectorXd duU = solver.solve(-R); //solve for the incremental displacements at the unknown degrees of freedom
            
            //construct full du vector including known values at dirischlet boundary
            for(int i = 0 ; i < unknownIndexes.size() ; i++){
                u(unknownIndexes[i]) += duU(i);
            }

            //Update the solution vector with the incremental displacements
            // u += du;
            cout << "Increment: " << increment+1 << ", Iteration: " << iter+1 << "\n";
            cout << "Modified residual norm: " << R.norm() << "\n"; //print the norm of the modified residual to monitor convergence of the unknown degrees of freedom
            cout << "-----------------------------------" << "\n";
            sendResidual(increment, iter, R.norm());
            if(R.norm() < epsilon){
                cout << "Convergence achieved for increment " << increment+1 << " in " << iter+1 << " iterations." << "\n";
                break; 
            }
        }

        std::string filename = "solutions/solution_" + std::to_string(increment+1) + ".vtu";
        write_vtu(filename, nodes, elements, u);
        solution_files.push_back(filename);
        solution_timesteps.push_back(increment+1);
    }

    write_pvd("solutions/final_solution.pvd", solution_files, solution_timesteps);

    end = high_resolution_clock::now();
    duration_msec = std::chrono::duration_cast<duration<double, std::milli>>(end-start);

    cout << "Total time taken (in ms): " << duration_msec.count() << endl;

}   