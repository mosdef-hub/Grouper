#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "dataStructures.hpp"
#include "processColoredGraphs.hpp"
#include "generate.hpp"
#include "autUtils.hpp"
#include "fragmentation.hpp"
#include <nauty/nauty.h>
#include <pybind11/stl.h>

#include <iostream>



namespace py = pybind11;


// Function to convert std::set<GroupGraph> to py::set
py::set convert_unordered_set(const std::unordered_set<GroupGraph>& cpp_set) {
    std::cout << "Converting unordered set" << std::endl;
    py::set py_set;
    for (const auto& item : cpp_set) {
        py_set.add(item);
    }
    return py_set;
}

PYBIND11_MODULE(_Grouper, m) {
    m.doc() = "Grouper bindings for Python";
    py::class_<GroupGraph::Node>(m, "Node")
        .def(py::init<>())
        .def(py::init<const std::string&, const std::string&, const std::vector<int>&>())
        .def_readwrite("type", &GroupGraph::Node::ntype)
        .def_readwrite("smarts", &GroupGraph::Node::smarts)
        .def_readwrite("ports", &GroupGraph::Node::ports)
        .def_readwrite("hubs", &GroupGraph::Node::hubs)
        .def("__eq__", &GroupGraph::Node::operator==)
        .def("__hash__", [](const GroupGraph::Node& node) {
            return std::hash<GroupGraph::Node>{}(node);
        });;
    py::class_<GroupGraph>(m, "GroupGraph")
        .def(py::init<>())
        .def_readwrite("nodes", &GroupGraph::nodes)
        .def_readwrite("edges", &GroupGraph::edges)
        .def_readwrite("node_types", &GroupGraph::nodetypes)
        .def("add_node", &GroupGraph::addNode,
             py::arg("type") = "",
             py::arg("smarts") = "",
             py::arg("hubs") = std::vector<int>{})
        .def("add_edge", &GroupGraph::addEdge,
             py::arg("src") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0},
             py::arg("dst") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0},
             py::arg("verbose") = false)
        .def("n_free_ports", &GroupGraph::n_free_ports)
        .def("compute_edge_orbits", &GroupGraph::computeEdgeOrbits)
        .def("__repr__", &GroupGraph::printGraph)
        .def("to_smiles", &GroupGraph::toSmiles, "Convert GroupGraph to SMILES")
        .def("to_vector", &GroupGraph::toVector, "Convert GroupGraph to group vector")
        .def("to_atom_graph", &GroupGraph::toAtomicGraph, "Convert GroupGraph to AtomGraph")
        .def("to_json", &GroupGraph::serialize, "Turn GroupGraph in JSON")
        .def("__hash__", [](const GroupGraph& g) {
            return std::hash<GroupGraph>{}(g);  // Using your defined hash function
        })
        .def("__eq__", &GroupGraph::operator==);
    py::class_<AtomGraph>(m, "AtomGraph")
        .def(py::init<>())
        .def_readwrite("nodes", &AtomGraph::nodes)
        .def_readwrite("edges", &AtomGraph::edges)
        .def("add_node", &AtomGraph::addNode,
             py::arg("type"),
             py::arg("valency") = 0)
        .def("add_edge", &AtomGraph::addEdge,
             py::arg("src"),
             py::arg("dst"),
             py::arg("order") = 1)
        .def("from_smiles", &AtomGraph::fromSmiles)
        .def("substructure_search", &AtomGraph::substructureSearch)
        .def("free_valency", &AtomGraph::getFreeValency)
        .def("__str__", &AtomGraph::printGraph)
        .def("__eq__", &AtomGraph::operator==);
    m.def("process_nauty_output", &process_nauty_output,
        py::arg("line"),
        py::arg("node_defs"),
        py::arg("graph_basis"),
        py::arg("positive_constraints"),
        py::arg("negative_constraints"),
        py::arg("verbose") = false,
        py::arg("g"),
        py::arg("lab"),
        py::arg("ptn"),
        py::arg("orbits"),
        py::arg("options"),
        py::arg("stats"));
    m.def("exhaustive_generate", [](int n_nodes,
                                    const std::unordered_set<GroupGraph::Node>& node_defs,
                                    const std::string& nauty_path,
                                    const std::string& input_file_path,
                                    int num_procs,
                                    const std::unordered_map<std::string, int>& positive_constraints,
                                    const std::unordered_set<std::string>& negative_constraints,
                                    const std::string& config_path,
                                    bool verbose) {
        std::unordered_set<GroupGraph> result = exhaustiveGenerate(n_nodes, node_defs, nauty_path, input_file_path, num_procs, positive_constraints, negative_constraints, config_path, verbose);
        return convert_unordered_set(result);
    },
        py::arg("n_nodes"),
        py::arg("node_defs"),
        py::arg("nauty_path"),
        py::arg("input_file_path") = "",
        py::arg("num_procs") = -1,
        py::arg("positive_constraints") = std::unordered_map<std::string, int>{},
        py::arg("negative_constraints") = std::unordered_set<std::string>{},
        py::arg("config_path") = "",
        py::arg("verbose") = false);
    m.def("fragment", &fragment,
        py::arg("smiles"),
        py::arg("node_defs"));
}
