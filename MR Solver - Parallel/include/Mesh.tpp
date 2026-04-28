#pragma once //include this file only once during compilation

//Member function definitions for Mesh class template
#include <fstream>
#include <string>
#include <filesystem>

template <unsigned int Nne>
void Mesh<Nne>::writeToFiles(const std::string& dir) const {
    std::filesystem::create_directories(dir); //creates dir folder if it does not exist

    std::ofstream points_file(dir + "/points.txt");
    for(const auto& node : nodes){
        points_file << node.x1 << " " << node.x2  << " " << node.x3 << "\n";
    }

    std::ofstream hexas_file(dir + "/elems.txt"); //triangulation file
    for(const auto& elem : elements){ 
        for(unsigned int i = 0 ; i < Nne ; i++){
            hexas_file << elem.node[i] << " ";
        }
        hexas_file << "\n";
    }
}