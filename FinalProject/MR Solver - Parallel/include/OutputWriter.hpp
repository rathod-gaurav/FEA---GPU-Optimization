#pragma once //include this only once during compilation

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Mesh.hpp"
#include <fstream>
#include <filesystem>

#include <sstream>
#include <curl/curl.h>

template <unsigned int Nne>
class OutputWriter{
    public:
        explicit OutputWriter(const std::string& outputDir); //constructor that takes the output directory as an argument

        std::string writeVTU(
            const Mesh<Nne>& mesh,
            const Eigen::VectorXd& u,
            unsigned int incr
        );

        void writePVD(
            const std::string& filename = "final_solution.pvd"
        ) const;

        void sendResidual(
            unsigned int incr,
            unsigned int iter,
            double residualNorm
        ) const;

    private:
        std::string outputDir_; //directory where the output files will be written
        std::vector<std::string> vtuFiles_; //list of VTU files generated during the simulation
        std::vector<double> timestamps_;
};

#include "OutputWriter.tpp" //include the implementation of the OutputWriter class