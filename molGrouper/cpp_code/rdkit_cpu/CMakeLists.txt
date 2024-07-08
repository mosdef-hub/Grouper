cmake_minimum_required(VERSION 3.18)

project(lineProcessor LANGUAGES CXX)

# Debug
# set(CMAKE_BUILD_TYPE Debug)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# # Specify the CUDA standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIE")
# set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Packages
find_package(RDKit REQUIRED)
find_package(Boost COMPONENTS system REQUIRED)
find_package(OpenMP REQUIRED)

set(SOURCES
    src/main.cpp
)

add_executable(lineProcessor ${SOURCES})

# Add OpenMP flags to the compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

find_library(RT_LIBRARY rt)

target_link_libraries(lineProcessor 
    RDKit::SmilesParse 
    RDKit::MolHash 
    RDKit::GraphMol 
    RDKit::RDGeneral 
    RDKit::FileParsers 
    ${OpenMP_CXX_FLAGS} ${OpenMP_CXX_LIBRARIES}
    rt
)

