# cmake_minimum_required(VERSION 3.18)

# project(lineProcessor LANGUAGES CXX)

# # Debug
# # set(CMAKE_BUILD_TYPE Debug)

# # Specify the C++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

# # # Specify the CUDA standard
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIE")
# # set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# # Packages
# find_package(RDKit REQUIRED)
# find_package(Boost COMPONENTS system REQUIRED)
# find_package(OpenMP REQUIRED)

# set(SOURCES
#     src/main.cpp
# )

# add_executable(lineProcessor ${SOURCES})

# # Add OpenMP flags to the compiler flags
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# find_library(RT_LIBRARY rt)

# target_link_libraries(lineProcessor 
#     RDKit::SmilesParse 
#     RDKit::MolHash 
#     RDKit::GraphMol 
#     RDKit::RDGeneral 
#     RDKit::FileParsers 
#     ${OpenMP_CXX_FLAGS} ${OpenMP_CXX_LIBRARIES}
#     rt
# )

cmake_minimum_required(VERSION 3.10)
project(lineProcessor)

set(CMAKE_BUILD_TYPE Debug)

# Add debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")

# Ensure RDBASE is set
if(NOT DEFINED ENV{RDBASE})
    message(FATAL_ERROR "RDBASE environment variable is not set")
endif()
set(RDBASE $ENV{RDBASE})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${RDBASE}/Code/cmake/Modules")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find RDKit package
find_package(RDKit REQUIRED)
# Find OpenMP package
set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/Cellar/libomp/18.1.8/include")
set(OpenMP_C_LIB_NAMES "omp")

# Add executable
add_executable(lineProcessor src/main.cpp)

# Include directories
include_directories(${RDBASE}/Code /opt/homebrew/Cellar/cairo/1.18.0/include /opt/homebrew/Cellar/boost/1.85.0/include /opt/homebrew/Cellar/libomp/18.1.8/include)

# Link RDKit libraries
link_directories(${RDBASE}/lib)
# target_compile_options(lineProcessor PRIVATE ${OpenMP_CXX_FLAGS})
target_link_libraries(lineProcessor
    ${RDBASE}/lib/libRDKitFileParsers.dylib
    ${RDBASE}/lib/libRDKitSmilesParse.dylib
    ${RDBASE}/lib/libRDKitGraphMol.dylib
    ${RDBASE}/lib/libRDKitRDGeneral.dylib
    /opt/homebrew/Cellar/libomp/18.1.8/lib/libomp.dylib)