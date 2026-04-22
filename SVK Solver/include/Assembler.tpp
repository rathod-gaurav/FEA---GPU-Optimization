#pragma once //include this only once during compilation

template <unsigned int Nne, unsigned int Nsd>
Assembler<Nne,Nsd>::Assembler(
    const Mesh<Nne>& mesh, const ElementEvaluator<Nne,Nsd>& elem_evaluator
) : mesh_(mesh), elem_evaluator_(elem_evaluator)
{}

template <unsigned int Nne, unsigned int Nsd>
Eigen::SparseMatrix<double> Assembler<Nne,Nsd>::extractSparseSubmatrix(
    const Eigen::SparseMatrix<double>& K,
    const std::vector<unsigned int>& rows,
    const std::vector<unsigned int>& cols) const {
        
    // Build a lookup set for fast membership testing
    std::unordered_set<unsigned int> rowSet(rows.begin(), rows.end());
    std::unordered_set<unsigned int> colSet(cols.begin(), cols.end());

    // Build index remapping: global index → local index in submatrix
    std::unordered_map<unsigned int,unsigned int> rowMap, colMap;
    for(unsigned int i = 0; i < rows.size(); i++) rowMap[rows[i]] = i;
    for(unsigned int j = 0; j < cols.size(); j++) colMap[cols[j]] = j;

    std::vector<Eigen::Triplet<double>> triplets;

    // Iterate over non-zeros of K
    for(int col = 0; col < K.outerSize(); col++){
        if(colSet.count(col) == 0) continue; // skip columns not in subset
        for(Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it){
            if(rowSet.count(it.row()) == 0) continue; // skip rows not in subset
            triplets.emplace_back(rowMap[it.row()], colMap[col], it.value());
        }
    }

    Eigen::SparseMatrix<double> sub(rows.size(), cols.size());
    sub.setFromTriplets(triplets.begin(), triplets.end());
    return sub;
}

template <unsigned int Nne, unsigned int Nsd>
void Assembler<Nne,Nsd>::assembleSystem(
            const Eigen::VectorXd& u, //global nodal displacement vector (Nnodes*Nsd x 1 vector)
            Eigen::SparseMatrix<double>& Kglobal, //global stiffness matrix (Nnodes*Nsd x Nnodes*Nsd sparse matrix)
            Eigen::VectorXd& Rglobal //global internal force vector (Nnodes*Nsd x 1 vector)
) const {
    unsigned int Nt = mesh_.Nnodes(); //total number of nodes in the mesh
    unsigned int Nel_t = mesh_.Nelements(); //total number of elements in the mesh

    Rglobal = Eigen::VectorXd::Zero(Nt * Nsd); //residual vector initialized to zero
    Kglobal = Eigen::SparseMatrix<double>(Nt * Nsd, Nt * Nsd); //sparse version of the tangent stiffness matrix for solving linear systems
    std::vector<Eigen::Triplet<double>> Kglobal_triplets; //triplet format for constructing the sparse tangent stiffness matrix
    Kglobal_triplets.reserve(Nel_t*Nne*Nne*9); //reserve space for triplets to avoid dynamic resizing during assembly

    Eigen::MatrixXd Klocal = Eigen::MatrixXd::Zero(Nne * Nsd, Nne * Nsd); //local tangent stiffness matrix for the element
    Eigen::VectorXd Rlocal = Eigen::VectorXd::Zero(Nne * Nsd); //local residual vector for the element

    //Loop over elements and assemble the global stiffness matrix and residual vector
    for(unsigned int e = 0 ; e < Nel_t ; e++){
        //Element nodal displacements
        Eigen::VectorXd u_e = Eigen::VectorXd::Zero(Nne * Nsd); //displacement vector for the current element
        for(unsigned int i = 0; i < Nne; i++){
            unsigned int global_node_id = mesh_.elements[e].node[i];
            u_e.segment(i*Nsd, Nsd) = u.segment(global_node_id*Nsd, Nsd); //extract the displacements for the nodes of the current element 
        }

        elem_evaluator_.computeElement(
            e, //element index
            u_e, //element nodal displacements (Nne*3 x 1 vector)
            Klocal, //element stiffness matrix (Nne*3 x Nne*3 matrix)
            Rlocal //element internal force vector (Nne*3 x 1 vector)
        );

        //Assemble Rlocal and Klocal into Rglobal and Kglobal
        for(int A = 0; A < Nne; A++){
            int Aglobal = mesh_.elements[e].node[A];
            for(int B = 0; B < Nne ; B++)
            {
                int Bglobal = mesh_.elements[e].node[B];
                for(int i = 0 ; i < 3 ; i++){
                    for(int j = 0 ; j < 3 ; j++){
                        Kglobal_triplets.emplace_back(3*Aglobal + i, 3*Bglobal + j, Klocal(3*A + i, 3*B + j));
                    }
                }
            }
            Rglobal.segment(3*Aglobal,3) += Rlocal.segment(3*A,3);
            // cout << "Assembled element " << e+1 << "/" << Nel_t << "\r";
        }
    }
    Kglobal.setFromTriplets(Kglobal_triplets.begin(), Kglobal_triplets.end()); //construct the sparse global tangent stiffness matrix from the triplets
    Kglobal.makeCompressed(); //compress the sparse matrix for efficient arithmetic and solving
}

template <unsigned int Nne, unsigned int Nsd>
void Assembler<Nne,Nsd>::partition(
    const Eigen::SparseMatrix<double>& Kglobal, //global stiffness matrix (Nnodes*Nsd x Nnodes*Nsd sparse matrix)
    Eigen::VectorXd& Rglobal, //global internal force vector (Nnodes*Nsd x 1 vector)
    const BoundaryConditions<Nne>& bcs, //boundary conditions object containing the indexes of the dirischlet DOFs

    Eigen::SparseMatrix<double>& KUU, //extract the submatrix of K corresponding to the unknown degrees of freedom
    Eigen::SparseMatrix<double>& KUD, //extract the submatrix of K corresponding to the coupling between unknown and dirischlet degrees of freedom
    Eigen::VectorXd& RU //extract the subvector of R corresponding to the unknown degrees of freedom
) const {
    const auto& dirischletIndexes = bcs.getDirischletIndexes(); //get the indexes of the dirischlet degrees of freedom from the boundary conditions object
    const auto& unknownIndexes = bcs.getUnknownIndexes(); //get the indexes of the unknown degrees of freedom from the boundary conditions object

    KUU = extractSparseSubmatrix(Kglobal, unknownIndexes, unknownIndexes); //extract the KUU submatrix corresponding to the unknown degrees of freedom
    KUD = extractSparseSubmatrix(Kglobal, unknownIndexes, dirischletIndexes); //extract the KUD submatrix corresponding to the coupling between unknown and dirischlet degrees of freedom
    
    RU.resize(unknownIndexes.size()); //extract the subvector of R corresponding to the unknown degrees of freedom
    for(int i = 0; i < unknownIndexes.size() ; i++){
        RU(i) = Rglobal(unknownIndexes[i]);
    }
}