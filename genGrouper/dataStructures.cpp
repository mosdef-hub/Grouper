#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cstring>

#include "dataStructures.hpp"
#include "autUtils.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>
#include <RDGeneral/types.h>

#include <nauty/nauty.h>
#include <nauty/naututil.h>




GroupGraph::GroupGraph()
    : nodes(), edges(), nodetypes() {}

GroupGraph::GroupGraph(const GroupGraph& other)
    : nodes(other.nodes), edges(other.edges), nodetypes(other.nodetypes) {}

GroupGraph& GroupGraph::operator=(const GroupGraph& other) {
    if (this != &other) {
        nodes = other.nodes;
        edges = other.edges;
        nodetypes = other.nodetypes;
    }
    return *this;
}

// Implementation of Node equality operator
bool GroupGraph::Node::operator==(const Node& other) const {
    return id == other.id &&
           ntype == other.ntype &&
           smiles == other.smiles &&
           ports == other.ports &&
           hubs == other.hubs;
}

bool GroupGraph::operator==(const GroupGraph& other) const {
    // Check if both graphs have the same number of nodes and edges
    if (nodes.size() != other.nodes.size() || edges.size() != other.edges.size()) {
        return false;
    }

    // Create adjacency lists for comparison using unordered_map and unordered_set
    std::unordered_map<NodeIDType, std::unordered_set<NodeIDType>> adjacency_list_1;
    std::unordered_map<NodeIDType, std::unordered_set<NodeIDType>> adjacency_list_2;

    // Populate adjacency lists for the current graph
    for (const auto& edge : edges) {
        NodeIDType node1 = std::get<0>(edge);
        NodeIDType node2 = std::get<2>(edge);
        adjacency_list_1[node1].insert(node2);
        adjacency_list_1[node2].insert(node1);
    }

    // Populate adjacency lists for the other graph
    for (const auto& edge : other.edges) {
        NodeIDType node1 = std::get<0>(edge);
        NodeIDType node2 = std::get<2>(edge);
        adjacency_list_2[node1].insert(node2);
        adjacency_list_2[node2].insert(node1);
    }

    // Create a vector of node IDs for the current graph and the other graph
    std::vector<NodeIDType> nodes_1;
    std::vector<NodeIDType> nodes_2;
    
    for (const auto& node : nodes) {
        nodes_1.push_back(node.first);
    }

    for (const auto& node : other.nodes) {
        nodes_2.push_back(node.first);
    }

    // Sort nodes to ensure consistent permutation checks
    std::sort(nodes_1.begin(), nodes_1.end());
    std::sort(nodes_2.begin(), nodes_2.end());

    // Try all permutations of node mappings
    do {
        std::unordered_map<NodeIDType, NodeIDType> node_mapping;

        for (size_t i = 0; i < nodes_1.size(); ++i) {
            node_mapping[nodes_1[i]] = nodes_2[i];
        }

        bool is_match = true;
        
        for (const auto& entry : adjacency_list_1) {
            NodeIDType node_id_1 = entry.first;
            const auto& neighbors_1 = entry.second;

            NodeIDType node_id_2 = node_mapping[node_id_1];
            const auto& neighbors_2 = adjacency_list_2[node_id_2];

            std::unordered_set<NodeIDType> mapped_neighbors_1;
            for (NodeIDType neighbor : neighbors_1) {
                mapped_neighbors_1.insert(node_mapping[neighbor]);
            }

            if (mapped_neighbors_1 != neighbors_2) {
                is_match = false;
                break;
            }
        }

        if (is_match) {
            return true;
        }

    } while (std::next_permutation(nodes_2.begin(), nodes_2.end()));

    return false;
}


// Non-member function to compare tuples (used for sorting edges)
inline bool operator<(const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType>& lhs,
                      const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType>& rhs) {
    // Use the built-in tuple comparison operator
    return std::tie(lhs) < std::tie(rhs);
}

void GroupGraph::addNode( 
    std::string ntype = "", 
    std::string smiles = "", 
    std::vector<PortType> ports = {}, 
    std::vector<NodeIDType> hubs = {}
) {
    // We assume that either you provide a type or a smiles string, if you provide smiles but not the type we will use it as the type
    int id;

    if (ntype == "" && smiles == "") {
        throw std::invalid_argument("Either smiles or type must be provided");
    }
    else if (ntype != "" && smiles == "") {
        if (nodetypes.find(ntype) != nodetypes.end()) {
            if (ports.size() == 0) {
                ports = nodetypes[ntype];
            }
            if (nodetypes[ntype] != ports) {
                throw std::invalid_argument("Node type already exists with different ports");
            }
            id = nodes.size();
        } 
       else {
            if (ports.size() == 0) {
                throw std::invalid_argument("Ports and hubs must be provided for new node type");
            };
            id = nodes.size();
            nodetypes[ntype] = ports;
        } 
    }
    else if (ntype == "" && smiles != "") {
        if (nodetypes.find(smiles) != nodetypes.end()) {
            if (nodetypes[smiles] != ports) {
                throw std::invalid_argument("Smiles already exists with different ports");
            }
            if (nodetypes[ntype] != ports) {
                throw std::invalid_argument("Smiles already exists with different ports");
            }
        }
        else {
            if (ports.size() == 0) {
                throw std::invalid_argument("Ports and hubs must be provided for new node type");
            };
            ntype = smiles.c_str();
            nodetypes[ntype] = ports;
        }
        
    }
    else {
        if (nodetypes.find(ntype) != nodetypes.end()) {
            if (ports.size() == 0) {
                ports = nodetypes[ntype];
            }
            if (nodetypes[ntype] != ports) {
                throw std::invalid_argument("Node type already exists with different ports");
            }
        }
        else {
            nodetypes[ntype] = ports;
        }
        
    }
    
    id = nodes.size();
    nodes[id] = Node(id, ntype, smiles, ports, hubs);
    
}

bool GroupGraph::addEdge(std::tuple<NodeIDType,PortType> fromNodePort, std::tuple<NodeIDType,PortType>toNodePort, bool verbose) {
    NodeIDType from = std::get<0>(fromNodePort);
    PortType fromPort = std::get<1>(fromNodePort);
    NodeIDType to = std::get<0>(toNodePort);
    PortType toPort = std::get<1>(toNodePort);

    if (n_free_ports(from) <= 0) {
        throw std::invalid_argument("Source node doesn't have enough ports!");
    }
    if (n_free_ports(to) <= 0) {
        throw std::invalid_argument("Destination node doesn't have enough ports!");
    }
    const std::tuple<NodeIDType, PortType, NodeIDType, PortType> edge = std::make_tuple(from, fromPort, to, toPort);
    if (std::find(edges.begin(), edges.end(), edge) != edges.end()) {
            throw std::invalid_argument("Edge already exists");
    }
    for (const auto& existingEdge : edges) {
        if (std::get<0>(existingEdge) == from && std::get<1>(existingEdge) == fromPort) {
            throw std::invalid_argument("Source port already in use");
        }
        if (std::get<2>(existingEdge) == to && std::get<3>(existingEdge) == toPort) {
            throw std::invalid_argument("Destination port already in use");
        }
    }
    if (nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) {
        throw std::invalid_argument("Node does not exist");
    }
    if (std::find(nodes[from].ports.begin(), nodes[from].ports.end(), fromPort) == nodes[from].ports.end()) {
        throw std::invalid_argument("Port does not exist");
    }
    edges.push_back(std::make_tuple(from, fromPort, to, toPort));
    return true;
}

int GroupGraph::n_free_ports(NodeIDType nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Node does not exist");
    }
    const Node& node = nodes.at(nodeID);
    int occupied_ports = 0;
    for (const auto& edge : edges) {
        if (std::get<0>(edge) == nodeID || std::get<2>(edge) == nodeID) {
            occupied_ports++;
        }
    }
    return node.ports.size() - occupied_ports;
}

int GroupGraph::numNodes() const {
    return nodes.size();
}

std::string GroupGraph::printGraph() const {
    std::ostringstream output;
    output << "Nodes:\n";
    for (const auto& entry : nodes) {
        output << "    Node " << entry.first << " (" << entry.second.ntype << ") (" << entry.second.smiles << ") ";
        output<< ": \n        Ports ";
        for (PortType port : entry.second.ports) {
            output << port << " ";
        }
        output << "\n        Hubs  ";
        for (NodeIDType hub : entry.second.hubs) {
            output << hub << " ";
        }
        output << "\n";
    }
    output << "Edges:\n";
    for (const auto& edge : edges) {
        output << "    Edge: " << std::get<0>(edge) << "(" << std::get<1>(edge) << ") -> " 
               << std::get<2>(edge) << "(" << std::get<3>(edge) << ")\n";
    }
    return output.str();
}

std::string GroupGraph::toSmiles() const {
    std::unique_ptr<RDKit::RWMol> molecularGraph(new RDKit::RWMol());
    std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
    int atomCount = 0;

    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Node& node = entry.second;
        std::string smiles = entry.second.smiles;
        std::unique_ptr<RDKit::RWMol> subGraph(RDKit::SmilesToMol(smiles));
        nodePortToAtomIndex[std::to_string(nodeID)] = std::unordered_map<int, int>();
        for (size_t i = 0; i < node.ports.size(); ++i) {
            nodePortToAtomIndex[std::to_string(nodeID)][node.ports[i]] = atomCount + node.hubs[i];
        }
        atomCount += subGraph->getNumAtoms();
    }


    int atomId = -1;
    std::unordered_map<std::string, std::unordered_map<int, int>> nodeSubGraphIndicesToMolecularGraphIndices;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Node& node = entry.second;
        nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)] = std::unordered_map<int, int>();
        std::string smiles = node.smiles;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
        for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*atom)->getIdx()] = atomId;
        }
    }


    atomId = -1;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Node& node = entry.second;
        std::string smiles = node.smiles;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
        for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            RDKit::Atom newAtom = **atom; 
            molecularGraph->addAtom(&newAtom, true);
        }
        for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
            RDKit::Bond newBond = **bond;
            int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getBeginAtomIdx()];
            int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getEndAtomIdx()];
            molecularGraph->addBond(atomIdx1, atomIdx2, newBond.getBondType());
        }
    }


    for (const auto& edge : edges) {
        NodeIDType from = std::get<0>(edge);
        PortType fromPort = std::get<1>(edge);
        NodeIDType to = std::get<2>(edge);
        PortType toPort = std::get<3>(edge);
        int fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
        int toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];

        molecularGraph->addBond(fromAtom, toAtom, RDKit::Bond::SINGLE);
    }
    
    return RDKit::MolToSmiles(*molecularGraph);
}

void GroupGraph::toNautyFormat(int *n, int *m, int *adj) const {
    // Set number of vertices
    *n = nodes.size();

    // Initialize adjacency matrix
    std::fill(adj, adj + (*n) * (*n), 0);
    
    // Fill adjacency matrix from edges
    for (const auto& edge : edges) {
        NodeIDType from = std::get<0>(edge);
        NodeIDType to = std::get<2>(edge);
        adj[from * (*n) + to] = 1;
        adj[to * (*n) + from] = 1; // Assuming undirected graph
    }

    *m = edges.size();
}

std::vector<std::vector<int>> GroupGraph::nodeAut() const {
    return {{}};
}
// std::vector<std::vector<int>> GroupGraph::nodeAut() const {
//     int n = nodes.size();  // Number of vertices
//     int m = SETWORDSNEEDED(n);  // Size of set words for NAUTY
    
//     // Dynamically allocate NAUTY structures
//     DYNALLSTAT(graph, g, g_sz);
//     DYNALLOC2(graph, g, g_sz, n, m, "malloc");
    
//     // Initialize the NAUTY graph to be empty
//     EMPTYGRAPH(g, m, n);

//     // Convert the adjacency matrix to the NAUTY graph representation
//     for (const auto& edge : edges) {
//         int u = std::get<0>(edge);
//         int v = std::get<2>(edge);
//         ADDELEMENT(GRAPHROW(g, u, m), v);
//         ADDELEMENT(GRAPHROW(g, v, m), u);
//     }

//     // NAUTY variables
//     int lab[MAXN], ptn[MAXN], orbits[MAXN];
//     optionblk options_struct;  // Define options structure
//     statsblk stats;

//     // Initialize options manually
//     options_struct.getcanon = FALSE;  // We only need automorphisms, not a canonical form
//     options_struct.defaultptn = TRUE;

//     // Workspace array for nauty
//     static set workspace[160 * MAXN];

//     // Active set (can be nullptr if not needed)
//     set *active = nullptr;

//     // Canonical form (not needed, so we set it to nullptr)
//     graph *canong = nullptr;

//     // Call NAUTY to compute automorphisms
//     nauty(g, lab, ptn, active, orbits, &options_struct, &stats, workspace, 160 * MAXN, m, n, canong);

//     // Collect results
//     std::vector<std::vector<int>> automorphisms;
//     for (int i = 0; i < stats.numorbits; ++i) {
//         std::vector<int> perm(n);
//         for (int j = 0; j < n; ++j) {
//             perm[j] = lab[j];
//         }
//         automorphisms.push_back(perm);
//     }

//     // Free dynamically allocated memory
//     DYNFREE(g, g_sz);

//     return automorphisms;
// }

std::vector<std::vector<std::pair<int, int>>> GroupGraph::edgeAut(const std::vector<std::pair<int, int>>& edge_list) const {
    int n = nodes.size();
    std::vector<std::vector<int>> adj_matrix(n, std::vector<int>(n, 0));
    
    createAdjMatrix(edge_list, adj_matrix);

    std::vector<std::pair<int, int>> edge_list_copy = edge_list;
    std::sort(edge_list_copy.begin(), edge_list_copy.end());

    std::vector<std::vector<std::pair<int, int>>> automorphisms;
    
    do {
        if (areEdgePermutationsEquivalent(edge_list, edge_list_copy)) {
            automorphisms.push_back(edge_list_copy);
        }
    } while (std::next_permutation(edge_list_copy.begin(), edge_list_copy.end()));

    return automorphisms;
}

std::unordered_map<std::string, int> GroupGraph::toVector() const {
    std::unordered_map<std::string, int> hist(nodetypes.size());
    for (const auto& entry : nodes) {
        const Node& node = entry.second;
        std::string ntype = node.ntype;
        hist[ntype] += 1;
    }
    return hist;
}

std::unique_ptr<AtomGraph> GroupGraph::toAtomicGraph() const {
    auto atomGraph = std::make_unique<AtomGraph>();
    const RDKit::PeriodicTable* pt = RDKit::PeriodicTable::getTable();
    std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
    int atomCount = 0;

    if (nodes.size() == 0) {
        throw std::invalid_argument("No nodes in the graph");
    }

    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Node& node = entry.second;
        std::string smiles = entry.second.smiles;
        std::unique_ptr<RDKit::RWMol> subGraph(RDKit::SmilesToMol(smiles));
        nodePortToAtomIndex[std::to_string(nodeID)] = std::unordered_map<int, int>();
        for (size_t i = 0; i < node.ports.size(); ++i) {
            nodePortToAtomIndex[std::to_string(nodeID)][node.ports[i]] = atomCount + node.hubs[i];
        }
        atomCount += subGraph->getNumAtoms();
    }



    int atomId = -1;
    std::unordered_map<std::string, std::unordered_map<int, int>> nodeSubGraphIndicesToMolecularGraphIndices;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Node& node = entry.second;
        nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)] = std::unordered_map<int, int>();
        std::string smiles = node.smiles;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
        for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*atom)->getIdx()] = atomId;
        }
    }

    // Start building the atom graph from mappings defined above
    atomId = -1;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Node& node = entry.second;
        std::string smiles = node.smiles;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
        for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            // RDKit::Atom newAtom = **atom; 
            // molecularGraph->addAtom(&newAtom, true);
            int atomicNumber = (*atom)->getAtomicNum();
            int maxValence = pt->getDefaultValence(atomicNumber);
            atomGraph->addNode((*atom)->getSymbol(), maxValence);

        }
        for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
            // RDKit::Bond newBond = **bond;
            int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getBeginAtomIdx()];
            int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getEndAtomIdx()];
            // molecularGraph->addBond(atomIdx1, atomIdx2, newBond.getBondType());
            atomGraph->addEdge(atomIdx1, atomIdx2);
        }
    }


    for (const auto& edge : edges) {
        NodeIDType from = std::get<0>(edge);
        PortType fromPort = std::get<1>(edge);
        NodeIDType to = std::get<2>(edge);
        PortType toPort = std::get<3>(edge);
        NodeIDType fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
        NodeIDType toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];
        // molecularGraph->addBond(fromAtom, toAtom, RDKit::Bond::SINGLE);
        atomGraph->addEdge(fromAtom, toAtom);
    }

    
    return atomGraph;
}

std::string GroupGraph::serialize() const {
        std::ostringstream oss;
        oss << "{\n  \"nodes\": [\n";
        for (const auto& pair : nodes) {
            const Node& node = pair.second;
            oss << "    {\n      \"id\": " << node.id 
                << ",\n      \"ntype\": \"" << node.ntype
                << "\",\n      \"smiles\": \"" << node.smiles
                << "\",\n      \"ports\": [";
            for (const auto& port : node.ports) {
                oss << port << ",";
            }
            oss.seekp(-1, oss.cur); // remove the last comma
            oss << "],\n      \"hubs\": [";
            for (const auto& hub : node.hubs) {
                oss << hub << ",";
            }
            oss.seekp(-1, oss.cur); // remove the last comma
            oss << "]\n    },\n";
        }
        oss.seekp(-2, oss.cur); // remove the last comma and newline
        oss << "\n  ],\n  \"edges\": [\n";
        for (const auto& edge : edges) {
            oss << "    [" 
                << std::get<0>(edge) << "," 
                << std::get<1>(edge) << "," 
                << std::get<2>(edge) << "," 
                << std::get<3>(edge) << "],\n";
        }
        oss.seekp(-2, oss.cur); // remove the last comma and newline
        oss << "\n  ]\n}";
        return oss.str();
    }

std::string GroupGraph::Canon() const {
        // TODO: Implement canonicalization algorithm, this doesn't work because the it doesn't account for automorphisms
        
        // Step 1: Sort nodes based on their attributes (e.g., id, ntype, smiles)
        std::vector<Node> sortedNodes;
        for (const auto& pair : nodes) {
            sortedNodes.push_back(pair.second);
        }
        std::sort(sortedNodes.begin(), sortedNodes.end(), [](const Node& a, const Node& b) {
            return std::tie(a.ntype, a.smiles, a.id) < std::tie(b.ntype, b.smiles, b.id);
        });

        // Step 2: Sort edges by connected nodes and ports
        std::vector<std::tuple<NodeIDType, PortType, NodeIDType, PortType>> sortedEdges = edges;
        std::sort(sortedEdges.begin(), sortedEdges.end());

        // Step 3: Convert sorted nodes and edges to a canonical SMILES or other format
        std::stringstream ss;
        for (const auto& node : sortedNodes) {
            ss << node.ntype << ":" << node.smiles << ";";
        }
        for (const auto& edge : sortedEdges) {
            ss << std::get<0>(edge) << "-" << std::get<1>(edge) << "-" 
               << std::get<2>(edge) << "-" << std::get<3>(edge) << ";";
        }

        // Step 4: Return the canonical representation
        return ss.str();
    }

//#############################################################################################################

AtomGraph::AtomGraph()
    : nodes(), edges() {}

AtomGraph::AtomGraph(const AtomGraph& other)
    : nodes(other.nodes), edges(other.edges) {}

AtomGraph& AtomGraph::operator=(const AtomGraph& other) {
    if (this != &other) {
        nodes = other.nodes;
        edges = other.edges;
    }
    return *this;
}

bool AtomGraph::operator==(const AtomGraph& other) const {
    // Check if both graphs have the same number of nodes and edges
    if (nodes.size() != other.nodes.size() || edges.size() != other.edges.size()) {
        return false;
    }

    // Create adjacency lists for comparison using unordered_map and unordered_set
    std::unordered_map<NodeIDType, std::unordered_set<NodeIDType>> adjacency_list_1;
    std::unordered_map<NodeIDType, std::unordered_set<NodeIDType>> adjacency_list_2;

    for (const auto& node : edges) {
        adjacency_list_1[node.first] = node.second;
    }

    for (const auto& node : other.edges) {
        adjacency_list_2[node.first] = node.second;
    }

    // Create a vector of node IDs
    std::vector<NodeIDType> nodes_1;
    std::vector<NodeIDType> nodes_2;
    
    for (const auto& node : nodes) {
        nodes_1.push_back(node.first);
    }

    for (const auto& node : other.nodes) {
        nodes_2.push_back(node.first);
    }

    // Try all permutations of node mappings
    std::sort(nodes_1.begin(), nodes_1.end());
    std::sort(nodes_2.begin(), nodes_2.end());

    do {
        std::unordered_map<NodeIDType, NodeIDType> node_mapping;

        for (size_t i = 0; i < nodes_1.size(); ++i) {
            node_mapping[nodes_1[i]] = nodes_2[i];
        }

        bool is_match = true;
        
        for (const auto& entry : adjacency_list_1) {
            NodeIDType node_id_1 = entry.first;
            const auto& neighbors_1 = entry.second;

            NodeIDType node_id_2 = node_mapping[node_id_1];
            const auto& neighbors_2 = adjacency_list_2[node_id_2];

            std::unordered_set<NodeIDType> mapped_neighbors_1;
            for (NodeIDType neighbor : neighbors_1) {
                mapped_neighbors_1.insert(node_mapping[neighbor]);
            }

            if (mapped_neighbors_1 != neighbors_2) {
                is_match = false;
                break;
            }
        }

        if (is_match) {
            return true;
        }

    } while (std::next_permutation(nodes_2.begin(), nodes_2.end()));

    return false;
}

void AtomGraph::addNode(const std::string& type, const unsigned int valency) {
    int id = nodes.size();
    nodes[id] = Node(id, type, valency);
    
}

void AtomGraph::addEdge(NodeIDType src, NodeIDType dst) {
    if (nodes.find(src) == nodes.end() || nodes.find(dst) == nodes.end()) {
        throw std::invalid_argument("Node does not exist");
    }
    if (getFreeValency(src) <= 0 && getFreeValency(dst) <= 0) {
        throw std::invalid_argument("Adding this edge would exceed the valency for the source node and the destination node");
    }
    if (getFreeValency(src) <= 0) {
        throw std::invalid_argument("Adding this edge would exceed the valency for the source node");
    }
    if (getFreeValency(dst) <= 0) {
        throw std::invalid_argument("Adding this edge would exceed the valency for the destination node");
    }
    if (edges[src].find(dst) != edges[src].end()) {
        throw std::invalid_argument("Edge already exists");
    }
    edges[src].insert(dst);
}

int AtomGraph::getFreeValency(NodeIDType nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Node does not exist");
    }
    const Node& node = nodes.at(nodeID);
    int occupied_electrons = edges.find(nodeID) != edges.end() ?  edges.at(nodeID).size() : 0;

    return node.valency - occupied_electrons;
}

std::string AtomGraph::printGraph() const {
    std::ostringstream output;
    output << "Nodes:\n";
    for (const auto& entry : nodes) {
        output << "    Node " << entry.first << " (" << entry.second.ntype << ")" << " Valency: " << entry.second.valency << "\n";
    }
    output << "Edges:\n";
    for (const auto& edge : edges) {
        for (const auto& dst : edge.second) {
            output << "    Edge: " << edge.first << " -> " << dst << "\n";
        }
    }
    return output.str();
}

// std::string AtomGraph::toSmiles(
// ) const {
// }




