cmake_minimum_required(VERSION 3.18)

project(lineProcessor LANGUAGES CXX CUDA)

# Debug
set(CMAKE_BUILD_TYPE Debug)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify the CUDA standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIE")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fPIE")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Packages
find_package(RDKit REQUIRED)
# find_package(CUDAToolkit REQUIRED)
find_package(CUDA REQUIRED)
find_package(Boost COMPONENTS system REQUIRED)

set(SOURCES
    src/main.cpp
    src/kernel.cu
)

add_executable(lineProcessor ${SOURCES})

# Include the CUDA directories
include_directories(${CUDA_INCLUDE_DIRS})

set_target_properties(lineProcessor PROPERTIES CUDA_SEPARABLE_COMPILATION ON)


target_link_libraries(lineProcessor 
    RDKit::SmilesParse 
    RDKit::MolHash 
    RDKit::GraphMol 
    RDKit::RDGeneral 
    RDKit::FileParsers 
    ${CUDA_LIBRARIES}
)

