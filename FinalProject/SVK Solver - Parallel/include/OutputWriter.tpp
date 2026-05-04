#pragma once //include this only once during compilation


template<unsigned int Nne>
OutputWriter<Nne>::OutputWriter(const std::string& outputDir)
: outputDir_(outputDir) //initialize the output directory member variable
{
    // Create the output directory if it doesn't exist
    std::filesystem::create_directories(outputDir_);
}

template<unsigned int Nne>
std::string OutputWriter<Nne>::writeVTU(
    const Mesh<Nne>& mesh,
    const Eigen::VectorXd& u,
    unsigned int incr
){
    std::string filename = outputDir_ + "/solution_"
                         + std::to_string(incr + 1) + ".vtu";

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "OutputWriter: could not open " << filename << "\n";
        return filename;
    }

    unsigned int Nn = mesh.Nnodes();
    unsigned int Ne = mesh.Nelements();

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "<UnstructuredGrid>\n"
         << "<Piece NumberOfPoints=\"" << Nn << "\" NumberOfCells=\"" << Ne << "\">\n";

    // --- Points ---
    file << "<Points>\n"
         << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < Nn; i++)
        file << mesh.nodes[i].x1 << " "
             << mesh.nodes[i].x2 << " "
             << mesh.nodes[i].x3 << "\n";
    file << "</DataArray>\n</Points>\n";

    // --- Cells ---
    file << "<Cells>\n"
         << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (unsigned int e = 0; e < Ne; e++) {
        for (unsigned int A = 0; A < Nne; A++)
            file << mesh.elements[e].node[A] << " ";
        file << "\n";
    }
    file << "</DataArray>\n"
         << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (unsigned int e = 0; e < Ne; e++) file << (e+1)*8 << " ";
    file << "\n</DataArray>\n"
         << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (unsigned int e = 0; e < Ne; e++) file << "12 "; // 12 = VTK hexahedron
    file << "\n</DataArray>\n</Cells>\n";

    // --- Point data: displacement ---
    file << "<PointData Vectors=\"Displacement\">\n"
         << "<DataArray type=\"Float32\" Name=\"Displacement\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (unsigned int i = 0; i < Nn; i++)
        file << u(3*i) << " " << u(3*i+1) << " " << u(3*i+2) << "\n";
    file << "</DataArray>\n</PointData>\n";

    file << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";

    vtuFiles_.push_back(filename);
    timestamps_.push_back(incr + 1);
    return filename;
}

template<unsigned int Nne>
void OutputWriter<Nne>::writePVD(
    const std::string& filename
) const {
    std::string path = outputDir_ + "/" + filename;
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "OutputWriter: could not open " << path << "\n";
        return;
    }

    file << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "<Collection>\n";
    for (size_t i = 0; i < vtuFiles_.size(); i++)
        file << "  <DataSet timestep=\"" << timestamps_[i]
             << "\" group=\"\" part=\"0\" file=\"" << vtuFiles_[i] << "\"/>\n";
    file << "</Collection>\n</VTKFile>\n";

    std::cout << "PVD written: " << path << "\n";
}

template<unsigned int Nne>
void OutputWriter<Nne>::sendResidual(
    unsigned int incr,
    unsigned int iter,
    double residualNorm
) const {

    CURL* curl = curl_easy_init();
    if (!curl) return;

    std::ostringstream json;
    json << "{\"increment\":" << incr
        << ",\"iteration\":" << iter
        << ",\"residual\":"  << residualNorm << "}";
    std::string body = json.str();

    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Content-Type: application/json");

    curl_easy_setopt(curl, CURLOPT_URL,       "http://127.0.0.1:8000/residual");
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT,    2L); // don't block the solver

    curl_easy_perform(curl);
    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
}