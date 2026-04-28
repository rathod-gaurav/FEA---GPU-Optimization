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
    const std::vector<unsigned int>& cols,
    int numThreads  // ← add this parameter
) const {

    // Build lookup structures — constructed serially once, then READ-ONLY
    // unordered_set/map are safe for concurrent reads from multiple threads
    std::unordered_set<unsigned int> rowSet(rows.begin(), rows.end());
    std::unordered_set<unsigned int> colSet(cols.begin(), cols.end());

    std::unordered_map<unsigned int,unsigned int> rowMap, colMap;
    for (unsigned int i = 0; i < rows.size(); i++) rowMap[rows[i]] = i;
    for (unsigned int j = 0; j < cols.size(); j++) colMap[cols[j]] = j;

    // One triplet vector per thread — same pattern as assembleSystem
    // Avoids all contention: each thread writes only to its own vector
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(numThreads);
    for (auto& v : thread_triplets)
        v.reserve(K.nonZeros() / numThreads + 1);

    // Outer loop: over columns of K (Eigen CSC — column is the natural outer index)
    // Each column is independent: colMap[col] is unique per column, no two
    // threads write to the same output column => no race condition
    #pragma omp parallel for schedule(dynamic, 64) num_threads(numThreads)
    for (int col = 0; col < K.outerSize(); col++) {
        if (colSet.count(col) == 0) continue; // not a column we want — skip

        int tid = omp_get_thread_num();
        int local_col = colMap.at(col); // .at() is read-only, thread-safe

        for (Eigen::SparseMatrix<double>::InnerIterator it(K, col); it; ++it) {
            if (rowSet.count(it.row()) == 0) continue;
            thread_triplets[tid].emplace_back(
                rowMap.at(it.row()),  // local row index
                local_col,            // local col index
                it.value()
            );
        }
    }

    // Serial merge — fires numThreads times, negligible cost
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(K.nonZeros());
    for (int t = 0; t < numThreads; t++) {
        triplets.insert(triplets.end(),
            std::make_move_iterator(thread_triplets[t].begin()),
            std::make_move_iterator(thread_triplets[t].end()));
    }

    Eigen::SparseMatrix<double> sub(rows.size(), cols.size());
    sub.setFromTriplets(triplets.begin(), triplets.end());
    return sub;
}

template <unsigned int Nne, unsigned int Nsd>
void Assembler<Nne,Nsd>::assembleSystem(
            const Eigen::VectorXd& u, //global nodal displacement vector (Nnodes*Nsd x 1 vector)
            Eigen::SparseMatrix<double>& Kglobal, //global stiffness matrix (Nnodes*Nsd x Nnodes*Nsd sparse matrix)
            Eigen::VectorXd& Rglobal, //global internal force vector (Nnodes*Nsd x 1 vector)
            int outerThreads
) const {
    unsigned int Nt = mesh_.Nnodes(); //total number of nodes in the mesh
    unsigned int Nel_t = mesh_.Nelements(); //total number of elements in the mesh

    Rglobal = Eigen::VectorXd::Zero(Nt * Nsd); //residual vector initialized to zero
    Kglobal = Eigen::SparseMatrix<double>(Nt * Nsd, Nt * Nsd); //sparse version of the tangent stiffness matrix for solving linear systems
    std::vector<Eigen::Triplet<double>> Kglobal_triplets; //triplet format for constructing the sparse tangent stiffness matrix
    Kglobal_triplets.reserve(Nel_t*Nne*Nne*9); //reserve space for triplets to avoid dynamic resizing during assembly

    const int numThreads = outerThreads; //set the outer number of threads for parallel execution

    //thread local triplet lists #optimization - use thread local storage for triplet lists to avoid contention during parallel assembly of the global stiffness matrix
    std::vector<std::vector<Eigen::Triplet<double>>> thread_local_triplets(numThreads); //thread local triplet lists
    for(auto& v : thread_local_triplets){
        v.reserve((Nel_t*Nne*Nne*9)/(numThreads+1)); //reserve space for triplets in each thread's local list to avoid dynamic resizing during assembly
    }

    //thread local residual #optimization - use thread local storage for residual vectors to avoid contention during parallel assembly of the global residual vector
    std::vector<Eigen::VectorXd> thread_local_R(numThreads); //thread local residual vectors
    #pragma omp parallel num_threads(numThreads)
    {
        int tid = omp_get_thread_num(); //get the thread ID for indexing into thread local storage
        thread_local_R[tid] = Eigen::VectorXd::Zero(Nt * Nsd); //initialize the thread local residual vector to zero
    }

    //first touch initialization of mesh data structures #optimization - parallelize the first-touch initialization of mesh data structures to ensure they are allocated on the correct NUMA nodes for better memory access performance during assembly
    // mesh_.elements and u are read by all threads during assembly.
    // This loop touches each element's data from the thread that will later process it (schedule(static) gives the same distribution as the assembly loop below). The OS migrates pages toward first-touch locality.
    
    #pragma omp parallel for schedule(static) num_threads(numThreads) //parallelize the first-touch initialization of mesh data structures with OpenMP, and use static scheduling to ensure the same distribution of iterations as the assembly loop below for better NUMA locality
    for(unsigned int e = 0 ; e < Nel_t ; e++){
        volatile unsigned int touch = mesh_.elements[e].node[0]; //touch the first node index of the element to initialize the corresponding page in memory on the thread that will later process this element during assembly
        (void)touch; //suppress unused variable warning
    }

    const int chunk = Nel_t/(numThreads*8); //chunk size for static scheduling of the assembly loop, ensuring each thread processes a contiguous block of elements for better cache performance and NUMA locality
    
    //Loop over elements and assemble the global stiffness matrix and residual vector
    #pragma omp parallel for schedule(dynamic, chunk) num_threads(numThreads) //#optimization - parallelize the assembly loop with OpenMP, and use dynamic scheduling with a chunk size to balance load while maintaining some locality
    for(unsigned int e = 0 ; e < Nel_t ; e++){
        int tid = omp_get_thread_num(); //get the thread ID for indexing into thread local storage

        //Klocal and Rlocal are thread local now
        Eigen::MatrixXd Klocal = Eigen::MatrixXd::Zero(Nne * Nsd, Nne * Nsd); //local tangent stiffness matrix for the element
        Eigen::VectorXd Rlocal = Eigen::VectorXd::Zero(Nne * Nsd); //local residual vector for the element

        //Element nodal displacements // because u is read only, it is safe for all threads to read from it without synchronization. Each thread will read the displacements for the nodes of the elements it processes, and there is no contention since they are only reading.
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

        //Assemble Rlocal and Klocal into thread local storage first to avoid contention, and then safely accumulate into the global residual and triplet list with a critical section
        for(int A = 0; A < Nne; A++){
            int Aglobal = mesh_.elements[e].node[A];
            for(int B = 0; B < Nne ; B++)
            {
                int Bglobal = mesh_.elements[e].node[B];
                for(int i = 0 ; i < 3 ; i++){
                    for(int j = 0 ; j < 3 ; j++){
                        thread_local_triplets[tid].emplace_back(3*Aglobal + i, 3*Bglobal + j, Klocal(3*A + i, 3*B + j));
                    }
                }
            }
            thread_local_R[tid].segment(3*Aglobal,3) += Rlocal.segment(3*A,3);
            // cout << "Assembled element " << e+1 << "/" << Nel_t << "\r";
        }
    }

    for(int t = 0 ; t < numThreads ; t++){
        //accumulate thread local triplets into global triplet list
        Kglobal_triplets.insert(Kglobal_triplets.end(), std::make_move_iterator(thread_local_triplets[t].begin()), std::make_move_iterator(thread_local_triplets[t].end())); //move the triplets from the thread local list into the global list to avoid copying, and then clear the thread local list to free memory
        //accumulate thread local residuals into global residual vector
        Rglobal += thread_local_R[t];
        thread_local_triplets[t].clear(); //clear the thread local triplet list to free memory
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
    Eigen::VectorXd& RU, //extract the subvector of R corresponding to the unknown degrees of freedom
    int numThreads
) const {
    const auto& dirischletIndexes = bcs.getDirischletIndexes(); //get the indexes of the dirischlet degrees of freedom from the boundary conditions object
    const auto& unknownIndexes = bcs.getUnknownIndexes(); //get the indexes of the unknown degrees of freedom from the boundary conditions object

    KUU = extractSparseSubmatrix(Kglobal, unknownIndexes, unknownIndexes, numThreads); //extract the KUU submatrix corresponding to the unknown degrees of freedom
    KUD = extractSparseSubmatrix(Kglobal, unknownIndexes, dirischletIndexes, numThreads); //extract the KUD submatrix corresponding to the coupling between unknown and dirischlet degrees of freedom
    
    RU.resize(unknownIndexes.size()); //extract the subvector of R corresponding to the unknown degrees of freedom
    for(int i = 0; i < unknownIndexes.size() ; i++){
        RU(i) = Rglobal(unknownIndexes[i]);
    }
}