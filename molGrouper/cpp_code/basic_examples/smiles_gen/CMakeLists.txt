cmake_minimum_required(VERSION 3.10)
project(smiles_gen)



set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/raid6/homes/kierannp/bin/rdkit/Code/cmake/Modules")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find RDKit package
find_package(RDKit REQUIRED)

# Add executable
add_executable(smiles_gen src/main.cpp)

# Include RDKit directories
include_directories(/raid6/homes/kierannp/bin/rdkit/Code /raid6/homes/kierannp/bin/boost /raid6/homes/kierannp/bin/cairo)

# Link RDKit libraries
link_directories(/raid6/homes/kierannp/bin/rdkit/lib)
target_link_libraries(smiles_gen
    /raid6/homes/kierannp/bin/rdkit/lib/libRDKitFileParsers.dylib
    /raid6/homes/kierannp/bin/rdkit/lib/libRDKitSmilesParse.dylib
    /raid6/homes/kierannp/bin/rdkit/lib/libRDKitGraphMol.dylib
    /raid6/homes/kierannp/bin/rdkit/lib/libRDKitRDGeneral.dylib)
