#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "groupGraph.cpp" 

namespace py = pybind11;

PYBIND11_MODULE(group_graph_bindings, m) {
    // Bind the GroupGraph class
    py::class_<GroupGraph>(m, "GroupGraph")
        .def(py::init<const std::unordered_map<std::string, std::vector<GroupGraph::PortType>>&>())
        .def("addNode", &GroupGraph::addNode)
        .def("addEdge", &GroupGraph::addEdge)
        .def("n_free_ports", &GroupGraph::n_free_ports)
        .def("make_undirected", &GroupGraph::make_undirected)
        .def("printGraph", &GroupGraph::printGraph);
}
