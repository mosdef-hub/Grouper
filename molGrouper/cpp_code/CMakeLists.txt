cmake_minimum_required(VERSION 3.12)
project(groupGraphBindings)

# Find Python
find_package(Python REQUIRED COMPONENTS Development)

# Specify the Python interpreter path
set(PYTHON_EXECUTABLE /Users/kieran/miniconda3/envs/molGrouper/bin/python3)

# Find Pybind11
find_package(pybind11 REQUIRED)

# Add the Pybind11 module
pybind11_add_module(groupGraphBindings groupGraphBindings.cpp)

# Link against the GroupGraph library (if it's defined in a separate source file)
# target_link_libraries(group_graph_bindings PRIVATE GroupGraph)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
