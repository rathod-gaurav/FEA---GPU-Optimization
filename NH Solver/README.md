cmake -B build : Once — only when you first set up the project, or if you change CMakeLists.txt
cmake --build build : Every time you change source code — recompiles only the changed files
./build/main : Every time you want to run