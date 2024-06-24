cmake_minimum_required(VERSION 3.10)
project(lineProcessor)

# set_source_files_properties(main.cpp PROPERTIES COMPILE_FLAGS "-g -O0")
# set_source_files_properties(GroupGraph.cpp PROPERTIES COMPILE_FLAGS "-g -O0")
# set_source_files_properties(process_colored_graphs.cpp PROPERTIES COMPILE_FLAGS "-g -O0")

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

# Add executable
add_executable(lineProcessor src/main.cpp)

# Include RDKit directories
include_directories(${RDBASE}/Code /opt/homebrew/Cellar/cairo/1.18.0/include /opt/homebrew/Cellar/boost/1.85.0/include)

# Link RDKit libraries
link_directories(${RDBASE}/lib)
target_link_libraries(lineProcessor
    ${RDBASE}/lib/libRDKitFileParsers.dylib
    ${RDBASE}/lib/libRDKitSmilesParse.dylib
    ${RDBASE}/lib/libRDKitGraphMol.dylib
    ${RDBASE}/lib/libRDKitRDGeneral.dylib)
