#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "dataStructures.hpp" 
#include "processColoredGraphs.hpp"
#include "generate.hpp"  

namespace py = pybind11;

PYBIND11_MODULE(_molGrouper, m) {
    m.doc() = "molGrouper bindings for Python";

    // Bind the GroupGraph::Node class
    py::class_<GroupGraph::Node>(m, "Node")
        .def(py::init<>())
        .def(py::init<int, const std::string&, const std::string&, const std::vector<int>&, const std::vector<int>&>())
        .def_readwrite("id", &GroupGraph::Node::id)
        .def_readwrite("type", &GroupGraph::Node::ntype)
        .def_readwrite("smiles", &GroupGraph::Node::smiles)
        .def_readwrite("ports", &GroupGraph::Node::ports)
        .def_readwrite("hubs", &GroupGraph::Node::hubs)
        .def("__eq__", &GroupGraph::Node::operator==)
        .def("__hash__", [](const GroupGraph::Node& node) {
            return std::hash<GroupGraph::Node>{}(node);
        });;

    // Bind the GroupGraph class
    py::class_<GroupGraph>(m, "GroupGraph")
        .def(py::init<>())
        .def_readwrite("nodes", &GroupGraph::nodes)
        .def_readwrite("edges", &GroupGraph::edges)
        .def_readwrite("node_types", &GroupGraph::nodetypes)
        .def("add_node", &GroupGraph::addNode, 
             py::arg("type") = "", 
             py::arg("smiles") = "", 
             py::arg("ports") = std::vector<int>{}, 
             py::arg("hubs") = std::vector<int>{})
        .def("add_edge", &GroupGraph::addEdge, 
             py::arg("from") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0}, 
             py::arg("to") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0},
             py::arg("verbose") = false)
        .def("n_nodes", &GroupGraph::numNodes)
        .def("n_free_ports", &GroupGraph::n_free_ports)
        .def("__str__", &GroupGraph::printGraph)
        .def("to_smiles", &GroupGraph::toSmiles, "Convert GroupGraph to SMILES")
        .def("to_vector", &GroupGraph::toVector, "Convert GroupGraph to group vector")
        .def("__eq__", &GroupGraph::operator==);

    m.def("process_nauty_output", &process_nauty_output, 
        py::arg("line"), 
        py::arg("node_defs"), 
        py::arg("verbose") = false);

    m.def("exhaustive_generate", &exhaustiveGenerate, 
        py::arg("n_nodes"), 
        py::arg("node_defs"), 
        py::arg("input_file_path"), 
        py::arg("num_procs") = 32, 
        py::arg("verbose") = false);
}
