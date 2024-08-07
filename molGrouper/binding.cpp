#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "dataStructures.h" 
#include "process_colored_graphs.cpp"

namespace py = pybind11;

PYBIND11_MODULE(molGrouper, m) {
    m.doc() = "molGrouper bindings for Python";
    // Bind the GroupGraph class
    py::class_<GroupGraph::Node>(m, "Node")
        .def(py::init<>())
        .def_readwrite("type", &GroupGraph::Node::ntype)
        .def_readwrite("smiles", &GroupGraph::Node::smiles)
        .def_readwrite("ports", &GroupGraph::Node::ports)
        .def_readwrite("hubs", &GroupGraph::Node::hubs);
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
             py::arg("from") = std::tuple<GroupGraph::NodeIDType,GroupGraph::PortType> {0,0}, 
             py::arg("to") = std::tuple<GroupGraph::NodeIDType,GroupGraph::PortType> {0,0},
             py::arg("verbose") = false)
        .def("n_nodes", &GroupGraph::numNodes)
        .def("n_free_ports", &GroupGraph::n_free_ports)
        .def("print", &GroupGraph::printGraph);
    // m.def("to_smiles", &GroupGraph::toSmiles, "Convert GroupGraph to SMILES",
    //     py::arg("node_type_to_smiles"), py::arg("node_type_port_to_index"));
    m.def("process_nauty_output", &process_nauty_output, 
        py::arg("line"), py::arg("node_types"), py::arg("node_int_to_port"), py::arg("verbose") = false);
}