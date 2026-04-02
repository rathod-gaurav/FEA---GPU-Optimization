cmake_minimum_required(VERSION 3.10)
project(MyProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include headers
include_directories("include")
include_directories("/home/grathod/Downloads/eigen-5.0.0/eigen-5.0.0") #include Eigen headers

# Gather all source files
file(GLOB SOURCES "src/*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -funroll-loops -std=c++17 -fopenmp") #add OpenMP flag for parallelization

# Add executable with main.cpp + all src files
add_executable(main main.cpp ${SOURCES})

# Link curl
target_link_libraries(main PRIVATE curl)