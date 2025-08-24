cmake_minimum_required(VERSION 3.15)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(pybind11 REQUIRED CONFIG)

# Get conda prefix
if(DEFINED ENV{CONDA_PREFIX})
    set(CONDA_PREFIX $ENV{CONDA_PREFIX})
else()
    message(FATAL_ERROR "CONDA_PREFIX environment variable not set")
endif()

# Platform-specific settings
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES "arm64")
    set(PLATFORM_COMPILE_FLAGS "-Xpreprocessor" "-fopenmp")
elseif(UNIX AND NOT APPLE)
    set(PLATFORM_COMPILE_FLAGS "-D_Thread_local=thread_local")
endif()

# Include directories
set(INCLUDE_DIRS
    ${CONDA_PREFIX}/include
    ${CONDA_PREFIX}/include/cairo
    ${CONDA_PREFIX}/include/boost
    ${CONDA_PREFIX}/include/rdkit
    ${CONDA_PREFIX}/include/omp
    ${CONDA_PREFIX}/include/nauty
    ${CONDA_PREFIX}/include/libpq
    ${CMAKE_CURRENT_SOURCE_DIR}/Grouper
)

# Library directories
set(LIBRARY_DIRS
    ${CONDA_PREFIX}/lib
    ${CONDA_PREFIX}/lib/cairo
)

# Libraries to link
set(LIBRARIES
    RDKitFileParsers
    RDKitSmilesParse
    RDKitGraphMol
    RDKitRDGeneral
    omp
    nauty
    pq
)

# Source files
set(SOURCES
    Grouper/dataStructures.cpp
    Grouper/binding.cpp
    Grouper/processColoredGraphs.cpp
    Grouper/generate.cpp
    Grouper/fragmentation.cpp
)

# Create pybind11 module
pybind11_add_module(_Grouper ${SOURCES})

# Set target properties
target_include_directories(_Grouper PRIVATE ${INCLUDE_DIRS})
target_link_directories(_Grouper PRIVATE ${LIBRARY_DIRS})
target_link_libraries(_Grouper PRIVATE ${LIBRARIES})

# Compiler flags
target_compile_options(_Grouper PRIVATE 
    -fPIC
    ${PLATFORM_COMPILE_FLAGS}
)

# Install the module
install(TARGETS _Grouper DESTINATION Grouper)