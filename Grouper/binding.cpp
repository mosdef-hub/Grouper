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
    // exceptions
    py::register_exception<GrouperParseException>(m, "GrouperParseException");
    py::register_exception<GrouperFragmentationError>(m, "GrouperFragmentationError");
    py::register_exception<GrouperNotYetImplementedException>(m, "GrouperNotYetImplementedException");
    py::class_<GroupGraph::Group>(m, "Group", R"doc(Represents a chemical group with a SMARTS pattern and connection ports.

.. code-block:: python

    from Grouper import Group

    # Create a methyl group
    methyl = Group("methyl", "[CH3]", [0])
)doc")
        .def(py::init<>(), "Default constructor.")
        .def(py::init<const std::string&, const std::string&, const std::vector<int>&, std::string>(),
            "Constructor with parameters.",
            py::arg("ntype"), py::arg("pattern"), py::arg("hubs")=py::list(), py::arg("pattern_type") = "SMILES")
        .def_readwrite("type", &GroupGraph::Group::ntype, "The type of the group.")
        .def_readwrite("pattern", &GroupGraph::Group::pattern, "The SMARTS pattern of the group.")
        .def_readwrite("ports", &GroupGraph::Group::ports, "The connection ports of the group.")
        .def_readwrite("hubs", &GroupGraph::Group::hubs, "The hub atoms of the group.")
        .def_readwrite("pattern_type", &GroupGraph::Group::patternType, "The type of the pattern (e.g., SMILES, SMARTS).")
        .def("possible_attachments", &GroupGraph::Group::getPossibleAttachments, "Get the number of possible attachments for the group.")
        .def("__eq__", &GroupGraph::Group::operator==)
        .def("__ne__", &GroupGraph::Group::operator!=)
        .def("__str__", &GroupGraph::Group::toString)
        .def("__repr__", &GroupGraph::Group::toString)
        .def("__hash__", [](const GroupGraph::Group& node) {
            return std::hash<GroupGraph::Group>{}(node);
        });;
    py::class_<GroupGraph>(m, "GroupGraph", R"doc(Represents a graph of chemical groups.

.. code-block:: python

    from Grouper import GroupGraph

    # Create an empty graph
    graph = GroupGraph()
)doc")
        .def(py::init<>(), "Default constructor.")
        .def_readwrite("nodes", &GroupGraph::nodes, "The nodes of the graph.")
        .def_readwrite("edges", &GroupGraph::edges, "The edges of the graph.")
        .def_readwrite("node_types", &GroupGraph::nodetypes, "The types of the nodes in the graph.")
        .def_readwrite("is_coarse_grained", &GroupGraph::isCoarseGrained, "Whether the graph is coarse-grained.")
        .def("add_node", static_cast<void (GroupGraph::*)(std::string, std::string, std::vector<GroupGraph::NodeIDType>, std::string)>(&GroupGraph::addNode),
            R"doc(Adds a node to the graph.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node("methyl", "[CH3]", [0])
)doc",
            py::arg("type") = "",
            py::arg("pattern") = "",
            py::arg("hubs") = std::vector<GroupGraph::NodeIDType>{},
            py::arg("pattern_type") = "SMILES"
        )
        .def("add_node", static_cast<void (GroupGraph::*)(GroupGraph::Group)>(&GroupGraph::addNode),
            "Adds a node to the graph from a Group object.",
            py::arg("group")
        )
        .def("add_edge", &GroupGraph::addEdge,
            R"doc(Adds an edge to the graph.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node("methyl", "[CH3]", [0])
    graph.add_node("hydroxyl", "[OH]", [0])
    graph.add_edge((0, 0), (1, 0))
)doc",
            py::arg("src") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0},
            py::arg("dst") = std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType>{0, 0},
            py::arg("order") = 1.,
            py::arg("strict") = true
        )
        .def("n_free_ports", &GroupGraph::numFreePorts, "Get the number of free ports in the graph.")
        .def("is_port_free", &GroupGraph::isPortFree, "Check if a port is free.", py::arg("node_id"), py::arg("port"))
        .def("clear_edges", &GroupGraph::clearEdges, "Remove all edges from the graph.")
        .def("compute_orbits",
            static_cast<std::tuple<std::vector<int>, std::vector<int>, std::vector<std::vector<int>>> (GroupGraph::*)(
                const std::vector<std::pair<int, int>>&,
                const std::vector<int>&
            ) const>(&GroupGraph::computeOrbits),
            R"doc(Compute the orbits of the graph.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node("ester", "C(O)=O", [0, 1])
    graph.add_node("ester", "C(O)=O", [0, 1])
    graph.add_node("ester", "C(O)=O", [0, 1])
    graph.add_node("not-ester", "C", [0, 0])

    edge_list = [ (0, 2),(0, 3),(1, 2),(1, 3) ]
    node_colors = [0,0,0,1]
    node_orbits, edge_orbits, _ = graph.compute_orbits(edge_list, node_colors)
)doc",
            py::arg("edge_list"),
            py::arg("node_colors")
        )
        .def("__repr__", &GroupGraph::printGraph)
        .def("to_smiles", &GroupGraph::toSmiles, R"doc(Convert GroupGraph to SMILES string.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node("methyl", "[CH3]", [0])
    graph.add_node("hydroxyl", "[OH]", [0])
    graph.add_edge((0, 0), (1, 0))
    # print(graph.to_smiles()) # Output: CO
)doc")
        .def("to_vector", &GroupGraph::toVector, R"doc(Convert GroupGraph to group vector.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node("methyl", "[CH3]", [0])
    graph.add_node("hydroxyl", "[OH]", [0])
    graph.add_edge((0, 0), (1, 0))
    # print(graph.to_vector()) # Output: {'methyl': 1, 'hydroxyl': 1}
)doc")
        .def("to_atom_graph", &GroupGraph::toAtomicGraph, R"doc(Convert GroupGraph to AtomGraph.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node(type="carbon", pattern= "C", hubs = [0,0,0,0])
    graph.add_node(type="nitrogen", pattern = "N", hubs = [0,0,0])
    graph.add_edge((0, 0), (1, 0))
    atom_graph = graph.to_atom_graph()
)doc")
        .def("to_json", &GroupGraph::serialize, R"doc(Serialize the GroupGraph to a JSON string.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node("type1", "C", [0, 0])
    graph.add_node("type1", "C", [0, 0])
    graph.add_edge((0, 0), (1, 0))
    json_str = graph.to_json()
)doc")
        // .def("to_canonical", &GroupGraph::canonize, "Canonicalize GroupGraph")
        .def("from_json", &GroupGraph::deserialize, R"doc(Deserialize a GroupGraph from a JSON string.

.. code-block:: python

    graph = GroupGraph()
    graph.add_node("type1", "C", [0, 0])
    graph.add_node("type1", "C", [0, 0])
    graph.add_edge((0, 0), (1, 0))
    json_str = graph.to_json()
    graph2 = GroupGraph()
    graph2.from_json(json_str)
)doc")
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
    py::class_<AtomGraph::Atom>(m, "Atom", R"doc(Represents an atom in an AtomGraph.

.. code-block:: python

    from Grouper import Atom

    # Create a carbon atom
    carbon = Atom("C", 4)
)doc")
        .def(py::init<>(), "Default constructor.")
        .def(py::init<std::string&>(), "Constructor with atom type.")
        .def(py::init<const std::string&, int>(), "Constructor with atom type and valency.")
        .def_readwrite("type", &AtomGraph::Atom::ntype, "The type of the atom (e.g., 'C', 'O').")
        .def_readwrite("valency", &AtomGraph::Atom::valency, "The valency of the atom.")
        .def("__eq__", &AtomGraph::Atom::operator==)
        .def("__str__", &AtomGraph::Atom::toString)
        .def("__repr__", &AtomGraph::Atom::toString)
        .def("__hash__", [](const AtomGraph::Atom& node) {
            return std::hash<AtomGraph::Atom>{}(node);
        });
    py::class_<AtomGraph>(m, "AtomGraph", R"doc(Represents a graph of atoms.

.. code-block:: python

    from Grouper import AtomGraph

    # Create an empty atom graph
    atom_graph = AtomGraph()
)doc")
        .def(py::init<>(), "Default constructor.")
        .def_readwrite("nodes", &AtomGraph::nodes, "The atoms in the graph.")
        .def_readwrite("edges", &AtomGraph::edges, "The bonds in the graph.")
        .def("add_node", static_cast<void (AtomGraph::*)(const std::string&, int)>(&AtomGraph::addNode),
            "Adds an atom to the graph.",
            py::arg("type"),
            py::arg("valency") = -1
        )
        .def("add_node", static_cast<void (AtomGraph::*)(AtomGraph::Atom)>(&AtomGraph::addNode),
            "Adds an atom to the graph from an Atom object.",
            py::arg("atom")
        )
        .def("add_edge", &AtomGraph::addEdge,
            "Adds a bond to the graph.",
            py::arg("src"),
            py::arg("dst"),
            py::arg("order") = 1.,
            py::arg("validate") = true)
        .def("from_smiles", &AtomGraph::fromSmiles, R"doc(Creates an AtomGraph from a SMILES string.

.. code-block:: python

    from Grouper import AtomGraph

    atom_graph = AtomGraph()
    atom_graph.from_smiles("CCO")
)doc")
        .def("from_smarts", &AtomGraph::fromSmarts, "Creates an AtomGraph from a SMARTS string.")
        .def("substructure_search", &AtomGraph::substructureSearch, R"doc(Performs a substructure search on the graph.

.. code-block:: python

    from Grouper import AtomGraph

    graph = AtomGraph()
    graph.from_smiles("CCO")
    sub = AtomGraph()
    sub.from_smiles("CO")
    matches = graph.substructure_search(sub)
)doc")
        .def("free_valency", &AtomGraph::getFreeValency, "Get the free valency of the graph.")
        .def("__str__", &AtomGraph::printGraph)
        .def("__eq__", &AtomGraph::operator==)
        .def("__hash__", [](const AtomGraph& g) {
            return std::hash<AtomGraph>{}(g);  // Using your defined hash function
        }
    );
    m.def("exhaustive_generate", [](int n_nodes,
                                    const std::unordered_set<GroupGraph::Group>& node_defs,
                                    int num_procs,
                                    const std::string& vcolg_output_file,
                                    const std::unordered_map<std::string, int>& positive_constraints,
                                    const std::unordered_set<std::string>& negative_constraints,
                                    const std::string& config_path) {
        auto result = exhaustiveGenerate(n_nodes, node_defs, num_procs, vcolg_output_file, positive_constraints, negative_constraints, config_path);
        return convert_unordered_set(result);
    },
    R"doc(Exhaustively generate all possible GroupGraphs with a given number of nodes and a set of allowed groups.

    .. code-block:: python

        from Grouper import Group, exhaustive_generate

        node_defs = [
            {"type": "t2", "smarts": "N", "hubs": [0, 0, 0]},
            {"type": "Methyl", "smarts": "C", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "C(=O)O", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "O", "hubs": [0, 0]},
        ]
        node_defs = set(
            Group(n["type"], n["smarts"], n["hubs"], pattern_type="SMILES")
            for n in node_defs
        )

        generated_graphs = exhaustive_generate(2, node_defs)
    )doc",
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
        auto result = randomGenerate(n_nodes, node_defs, num_graphs, num_procs, positive_constraints, negative_constraints);
        return convert_unordered_set(result);
    },
    R"doc(Randomly generate a specified number of GroupGraphs.

    .. code-block:: python

        from Grouper import Group, random_generate

        node_defs = [
            {"type": "t2", "smarts": "N", "hubs": [0, 0, 0]},
            {"type": "Methyl", "smarts": "C", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "C(=O)(O)", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "O", "hubs": [0, 0]},
        ]
        node_defs = set(
            Group(n["type"], n["smarts"], n["hubs"], pattern_type="SMILES")
            for n in node_defs
        )

        generated_graphs = random_generate(2, node_defs, 5)
    )doc",
        py::arg("n_nodes"),
        py::arg("node_defs"),
        py::arg("num_graphs") = 100,
        py::arg("num_procs") = -1,
        py::arg("positive_constraints") = std::unordered_map<std::string, int>{},
        py::arg("negative_constraints") = std::unordered_set<std::string>{}
    );
    m.def("exhaustive_fragment", &fragment,
        R"(Fragment a molecule into a set of GroupGraphs based on a given set of group definitions.

        .. code-block:: python

            from Grouper import Group, exhaustive_fragment

            node_defs = set()
            node_defs.add(Group("hydroxyl", "O", [0]))
            node_defs.add(Group("alkene", "C=C", [0, 0, 1, 1]))

            smiles = "C=CO"

            fragmented_graphs = exhaustive_fragment(smiles, node_defs)
        )",
        py::arg("smiles"),
        py::arg("node_defs")
    );
}
