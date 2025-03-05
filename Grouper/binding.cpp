#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "dataStructures.hpp"
#include "processColoredGraphs.hpp"
#include "generate.hpp"
#include "fragmentation.hpp"
#include <nauty/nauty.h>
#include <pybind11/stl.h>

#include <iostream>



namespace py = pybind11;


// Function to convert std::set<GroupGraph> to py::set
py::set convert_unordered_set(const std::unordered_set<GroupGraph>& cpp_set) {
    py::set py_set;
    for (const auto& item : cpp_set) {
        py_set.add(item);
    }
    return py_set;
}

PYBIND11_MODULE(_Grouper, m) {
    m.doc() = "Grouper bindings for Python";
    py::class_<GroupGraph::Group>(m, "Group")
        .def(py::init<>())
        .def(py::init<const std::string&, const std::string&, const std::vector<int>&, bool>(),
            py::arg("ntype"), py::arg("pattern"), py::arg("hubs"), py::arg("is_smarts") = false)
        .def_readwrite("type", &GroupGraph::Group::ntype)
        .def_readwrite("pattern", &GroupGraph::Group::pattern)
        .def_readwrite("ports", &GroupGraph::Group::ports)
        .def_readwrite("hubs", &GroupGraph::Group::hubs)
        .def_readwrite("is_smarts", &GroupGraph::Group::isSmarts)
        .def("compute_hub_orbits", &GroupGraph::Group::hubOrbits)
        .def("possible_attachments", &GroupGraph::Group::getPossibleAttachments)
        .def("__eq__", &GroupGraph::Group::operator==)
        .def("__ne__", &GroupGraph::Group::operator!=)
        .def("__str__", &GroupGraph::Group::toString)
        .def("__repr__", &GroupGraph::Group::toString)
        .def("__hash__", [](const GroupGraph::Group& node) {
            return std::hash<GroupGraph::Group>{}(node);
        });;
    py::class_<GroupGraph>(m, "GroupGraph")
        .def(py::init<>())
        .def_readwrite("nodes", &GroupGraph::nodes)
        .def_readwrite("edges", &GroupGraph::edges)
        .def_readwrite("node_types", &GroupGraph::nodetypes)
        .def("add_node", &GroupGraph::addNode,
             py::arg("type") = "",
             py::arg("pattern") = "",
             py::arg("hubs") = std::vector<int>{},
             py::arg("isSmarts") = false
        )
        .def("add_edge", &GroupGraph::addEdge,
             py::arg("src") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0},
             py::arg("dst") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0},
             py::arg("order") = 1,
             py::arg("verbose") = false
        )
        .def("n_free_ports", &GroupGraph::numFreePorts)
        .def("compute_orbits",
            static_cast<std::pair<std::vector<int>, std::vector<int>> (GroupGraph::*)(
                const std::vector<std::pair<int, int>>&,
                const std::vector<int>&
            ) const>(&GroupGraph::computeOrbits),
            py::arg("edge_list"),
            py::arg("node_colors")
        )
        .def("__repr__", &GroupGraph::printGraph)
        .def("to_smiles", &GroupGraph::toSmiles, "Convert GroupGraph to SMILES")
        .def("to_vector", &GroupGraph::toVector, "Convert GroupGraph to group vector")
        .def("to_atom_graph", &GroupGraph::toAtomicGraph, "Convert GroupGraph to AtomGraph")
        .def("to_json", &GroupGraph::serialize, "Turn GroupGraph in JSON")
        .def("to_canonical", &GroupGraph::canonize, "Canonicalize GroupGraph")
        .def("from_json", &GroupGraph::deserialize, "Load GroupGraph from JSON")
        .def("__hash__", [](const GroupGraph& g) {
            return std::hash<GroupGraph>{}(g);  // Using your defined hash function
        })
        .def("__eq__", &GroupGraph::operator==)
        .def(py::pickle(
            [](const GroupGraph& g) {
                return g.serialize();
            },
            [](const std::string& json_state) {
                GroupGraph g;
                g.deserialize(json_state);
                return g;
            }
        ));
    py::class_<AtomGraph::Atom>(m, "Atom")
        .def(py::init<>())
        .def(py::init<std::string&>())
        .def(py::init<const std::string&, int>())
        .def_readwrite("type", &AtomGraph::Atom::ntype)
        .def_readwrite("valency", &AtomGraph::Atom::valency)
        .def("__eq__", &AtomGraph::Atom::operator==)
        .def("__str__", &AtomGraph::Atom::toString)
        .def("__repr__", &AtomGraph::Atom::toString)
        .def("__hash__", [](const AtomGraph::Atom& node) {
            return std::hash<AtomGraph::Atom>{}(node);
        });
    py::class_<AtomGraph>(m, "AtomGraph")
        .def(py::init<>())
        .def_readwrite("nodes", &AtomGraph::nodes)
        .def_readwrite("edges", &AtomGraph::edges)
        .def("add_node", &AtomGraph::addNode,
             py::arg("type"),
             py::arg("valency") = -1)
        .def("add_edge", &AtomGraph::addEdge,
             py::arg("src"),
             py::arg("dst"),
             py::arg("order") = 1)
        .def("from_smiles", &AtomGraph::fromSmiles)
        .def("from_smarts", &AtomGraph::fromSmarts)
        .def("substructure_search", &AtomGraph::substructureSearch)
        .def("to_canonical", &AtomGraph::canonize)
        .def("free_valency", &AtomGraph::getFreeValency)
        .def("__str__", &AtomGraph::printGraph)
        .def("__eq__", &AtomGraph::operator==)
        .def("__hash__", [](const AtomGraph& g) {
            return std::hash<AtomGraph>{}(g);  // Using your defined hash function
        }
    );
    m.def("exhaustive_generate", [](int n_nodes,
                                    const std::unordered_set<GroupGraph::Group>& node_defs,
                                    const std::string& vcolg_output_file,
                                    int num_procs,
                                    const std::unordered_map<std::string, int>& positive_constraints,
                                    const std::unordered_set<std::string>& negative_constraints,
                                    const std::string& config_path) {
        std::unordered_set<GroupGraph> result = exhaustiveGenerate(n_nodes, node_defs, num_procs, vcolg_output_file, positive_constraints, negative_constraints, config_path);
        return convert_unordered_set(result);
    },
        py::arg("n_nodes"),
        py::arg("node_defs"),
        py::arg("num_procs") = -1,
        py::arg("vcolg_output_file") = "",
        py::arg("positive_constraints") = std::unordered_map<std::string, int>{},
        py::arg("negative_constraints") = std::unordered_set<std::string>{},
        py::arg("config_path") = ""
    );
    m.def("random_generate", [](int n_nodes,
                                const std::unordered_set<GroupGraph::Group>& node_defs,
                                int num_graphs,
                                int num_procs,
                                const std::unordered_map<std::string, int>& positive_constraints,
                                const std::unordered_set<std::string>& negative_constraints) {
        std::unordered_set<GroupGraph> result = randomGenerate(n_nodes, node_defs, num_graphs, num_procs, positive_constraints, negative_constraints);
        return convert_unordered_set(result);
    },
        py::arg("n_nodes"),
        py::arg("node_defs"),
        py::arg("num_graphs") = 100,
        py::arg("num_procs") = -1,
        py::arg("positive_constraints") = std::unordered_map<std::string, int>{},
        py::arg("negative_constraints") = std::unordered_set<std::string>{}
    );
    m.def("fragment", &fragment,
        py::arg("smiles"),
        py::arg("node_defs")
    );
}
