#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "groupGraph.cpp" 
#include "processMultig.cpp"

namespace py = pybind11;

PYBIND11_MODULE(groupGraphBindings, m) {
    m.doc() = "GroupGraph bindings";
    // Bind the GroupGraph class
    py::class_<GroupGraph::Node>(m, "Node")
        .def(py::init<>())
        .def_readwrite("type", &GroupGraph::Node::type)
        .def_readwrite("ports", &GroupGraph::Node::ports);
    py::class_<GroupGraph>(m, "GroupGraph")
        .def(py::init<const std::unordered_map<std::string, std::vector<GroupGraph::PortType>>&>())
        .def_readwrite("nodes", &GroupGraph::nodes)
        .def_readwrite("edges", &GroupGraph::edges)
        .def_readwrite("node_types", &GroupGraph::nodeTypes)
        .def("add_node", &GroupGraph::addNode)
        .def("add_edge", &GroupGraph::addEdge)
        .def("n_nodes", &GroupGraph::numNodes)
        .def("n_free_ports", &GroupGraph::n_free_ports)
        .def("print", &GroupGraph::printGraph);

    // Bind the MultigConverter class
    py::class_<MultigConverter>(m, "MultigConverter")
        .def(py::init<>())
        .def_static("multi_to_pair", &MultigConverter::multi_to_pair, "Convert multi to pair")
        .def_static("parse_multig_file", &MultigConverter::parse_multig_file, "return group graphs from multig text file")
        .def_static("multig_line_to_graph", &MultigConverter::multig_line_to_graph, "Convert multig output to GroupGraph",
            py::arg("line"), py::arg("int_to_node_type"), py::arg("node_int_to_port"), py::arg("node_types"), py::arg("verbose") = false);
}
