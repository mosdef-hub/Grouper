#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <queue>
#include <stack>

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


// Define structures as thread-local using thread_local
thread_local std::vector<setword> g;         // For the nauty graph
thread_local std::vector<int> lab;          // For the label array
thread_local std::vector<int> ptn;          // For the partition array
thread_local std::vector<int> orbits;       // For the orbits array
thread_local DEFAULTOPTIONS_GRAPH(options); // Default nauty options
thread_local statsblk stats;                // Nauty stats structure

// Function to initialize the thread-local nauty structures
void initializeNautyStructures(int n) {
    int m = SETWORDSNEEDED(n);

    // Resize the vectors based on the required size
    g.resize(m * n); // The nauty graph requires m * n words
    lab.resize(n);   // Label array size is n
    ptn.resize(n);   // Partition array size is n
    orbits.resize(n); // Orbit array size is n

    // Clear the graph
    std::fill(g.begin(), g.end(), 0);
}

// Core methods
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

bool GroupGraph::Group::operator==(const Group& other) const {
    return id == other.id &&
           ntype == other.ntype &&
           smarts == other.smarts &&
           hubs == other.hubs;
}

bool GroupGraph::Group::operator!=(const Group& other) const {
    return !(*this == other);
}

std::vector<int> GroupGraph::Group::hubOrbits() const {
    return hubs; // TODO: Implement hub orbits
}

bool GroupGraph::operator==(const GroupGraph& other) const {
    // Check if empty
    if (nodes.empty() && other.nodes.empty()) {
        return true;
    }
    if (nodes.empty() || other.nodes.empty()) {
        return false;
    }

    // Convert GroupGraph to AtomGraph
    std::unique_ptr<AtomGraph> atomGraph1 = this->toAtomicGraph();
    std::unique_ptr<AtomGraph> atomGraph2 = other.toAtomicGraph();

    // Check if the number of nodes and edges are the same
    if (atomGraph1->nodes.size() != atomGraph2->nodes.size()) {
        return false;
    }
    int edgeCount1 = 0;
    int edgeCount2 = 0;
    for (const auto& node : atomGraph1->edges) {
        edgeCount1 += node.second.size();
    }
    for (const auto& node : atomGraph2->edges) {
        edgeCount2 += node.second.size();
    }
    if (edgeCount1 != edgeCount2) {
        return false;
    }

    // Convert AtomGraph to nauty graph
    int n = atomGraph1->nodes.size(); // Assuming the number of nodes is the same for both graphs
    int m = SETWORDSNEEDED(n);

    // Use std::vector instead of DYNALLSTAT and DYNALLOC
    std::vector<setword> g1(m * n, 0); // Initialize graph 1
    std::vector<setword> g2(m * n, 0); // Initialize graph 2
    std::vector<int> lab1(n), ptn1(n), orbits1(n); // Label, partition, and orbits for graph 1
    std::vector<int> lab2(n), ptn2(n), orbits2(n); // Label, partition, and orbits for graph 2
    std::vector<setword> canong1(m * n, 0); // Canonical form for graph 1
    std::vector<setword> canong2(m * n, 0); // Canonical form for graph 2
    setword workspace[160]; // Workspace for nauty

    // Initialize nauty structures
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    // Convert AtomGraph to nauty graph for g1
    EMPTYGRAPH(g1.data(), m, n);
    for (const auto& [id, dst_order] : atomGraph1->edges) {
        if (!dst_order.empty()) {
            for (const auto& [dest,order] : dst_order) {
                int from = id;
                int to = dest;
                ADDONEEDGE(g1.data(), from, to, m);
            }
        }
    }

    // Convert AtomGraph to nauty graph for g2
    EMPTYGRAPH(g2.data(), m, n);
    for (const auto& [id, dst_order] : atomGraph2->edges) {
        if (!dst_order.empty()) {
            for (const auto& [dest,order] : dst_order) {
                int from = id;
                int to = dest;
                ADDONEEDGE(g2.data(), from, to, m);
            }
        }
    }

    // Call nauty to canonicalize the graphs
    options.getcanon = TRUE;
    nauty(g1.data(), lab1.data(), ptn1.data(), nullptr, orbits1.data(), &options, &stats, workspace, 160, m, n, canong1.data());
    nauty(g2.data(), lab2.data(), ptn2.data(), nullptr, orbits2.data(), &options, &stats, workspace, 160, m, n, canong2.data());

    // Compare the canonical forms to determine isomorphism
    if (memcmp(canong1.data(), canong2.data(), sizeof(setword) * m * n) != 0) {
        return false;
    }

    return true;
}

inline bool operator<(const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>& lhs,
                      const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>& rhs) {
    // Use the built-in tuple comparison operator
    return std::tie(lhs) < std::tie(rhs);
}

// Operating methods
void GroupGraph::addNode(
    std::string ntype = "",
    std::string smarts = "",
    std::vector<NodeIDType> hubs = {}
) {
    /*There are 5 ways you can input a node:
        1. ntype, smiles, hubs
        2. ntype, hubs
        3. smiles, hubs
        4. ntype if it already exists
        5. smarts if it already exists
    */

    // Error handling
    if (ntype.empty() && smarts.empty()) {
        throw std::invalid_argument("Either SMARTS or type must be provided");
    }
    for (NodeIDType hub : hubs) {
        if (hub < 0) {
            throw std::invalid_argument("Hub ID must be greater than or equal to 0");
        }
    }

    // Check if the node type already exists
    if (ntype.empty()) {
        ntype = smarts;
    }

    // Initialize id and ports
    int id = nodes.size();
    std::vector<PortType> ports(hubs.size());

    // Case 0: Group type and hubs are provided
    if (!ntype.empty() && nodetypes.find(ntype) != nodetypes.end() && hubs.empty()) {
        hubs = nodetypes[ntype]; // Reuse existing hubs
    }

    // Case 1: Group type (ntype) is provided
    if (!ntype.empty()) {
        // If the ntype (or smarts used as ntype) already exists
        if (nodetypes.find(ntype) != nodetypes.end()) {
            // Check if hub sizes match the existing node type
            if (hubs.empty()) {
                ports = nodetypes[ntype]; // Use the existing ports
            }
            if (nodetypes[ntype].size() != hubs.size()) {
                throw std::invalid_argument("Group type already exists with a different number of hubs");
            }
        } else { // New node type
            std::iota(ports.begin(), ports.end(), 0); // Initialize ports for hubs
            nodetypes[ntype] = ports; // Save the new node type
        }
    }

    // Case 2: smarts is provided (and used as ntype if ntype was empty)
    if (ntype.empty() && !smarts.empty() && nodetypes.find(smarts) != nodetypes.end()) {
        // Check if hub sizes match the existing smarts node
        if (nodetypes[smarts].size() != hubs.size()) {
            throw std::invalid_argument("smarts already exists with a different number of hubs");
        }
    } else if (!smarts.empty() && nodetypes.find(smarts) == nodetypes.end()) {
        std::iota(ports.begin(), ports.end(), 0);
        nodetypes[smarts] = ports; // Save the new smarts as a node type
    }

    // Create the new node and add it to the nodes map
    nodes[id] = Group(ntype, smarts, hubs);
}


bool GroupGraph::addEdge(std::tuple<NodeIDType,PortType> fromNodePort, std::tuple<NodeIDType,PortType>toNodePort, unsigned int bondOrder, bool verbose) {
    NodeIDType from = std::get<0>(fromNodePort);
    PortType fromPort = std::get<1>(fromNodePort);
    NodeIDType to = std::get<0>(toNodePort);
    PortType toPort = std::get<1>(toNodePort);


    // Error handling
    if (numFreePorts(from) <= 0) {
        throw std::invalid_argument("Source node doesn't have enough ports!");
    }
    if (numFreePorts(to) <= 0) {
        throw std::invalid_argument("Destination node doesn't have enough ports!");
    }
    if (from == to) {
        throw std::invalid_argument("Source and destination nodes are the same");
    }
    if (nodes.find(from) == nodes.end()) {
        throw std::invalid_argument("Source node does not exist");
    }
    if (nodes.find(to) == nodes.end()) {
        throw std::invalid_argument("Destination node does not exist");
    }
    if (std::find(nodes[from].ports.begin(), nodes[from].ports.end(), fromPort) == nodes[from].ports.end()) {
        throw std::invalid_argument("Source port does not exist");
    }
    if (std::find(nodes[to].ports.begin(), nodes[to].ports.end(), toPort) == nodes[to].ports.end()) {
        throw std::invalid_argument("Destination port does not exist");
    }
    const std::tuple<NodeIDType, PortType, NodeIDType, PortType, unsigned int> edge = std::make_tuple(from, fromPort, to, toPort, bondOrder);
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

    // Add the edge
    edges.push_back(std::make_tuple(from, fromPort, to, toPort, bondOrder));
    return true;
}

int GroupGraph::numFreePorts(NodeIDType nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Can't calculate numFreePorts because node does not exist");
    }
    const Group& node = nodes.at(nodeID);
    int occupied_ports = 0;
    for (const auto& edge : edges) {
        if (std::get<0>(edge) == nodeID || std::get<2>(edge) == nodeID) {
            occupied_ports++;
        }
    }
    return node.ports.size() - occupied_ports;
}

bool GroupGraph::isPortFree(NodeIDType nodeID, PortType port) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Can't check if port is free because node does not exist");
    }
    if (std::find(nodes.at(nodeID).ports.begin(), nodes.at(nodeID).ports.end(), port) == nodes.at(nodeID).ports.end()) {
        throw std::invalid_argument("Can't check if port is free because port does not exist");
    }
    for (const auto& edge : edges) {
        if ((std::get<0>(edge) == nodeID && std::get<1>(edge) == port) ||
            (std::get<2>(edge) == nodeID && std::get<3>(edge) == port)) {
            return false;
        }
    }
    return true;
}

void GroupGraph::clearEdges() {
    edges.clear();
}

int* GroupGraph::computeEdgeOrbits(
    const std::vector<std::pair<int, int>> edge_list,
    graph* g, int* lab, int* ptn, int* orbits,
    optionblk* options,
    statsblk* stats
) const {
    std::vector<std::vector<int>> edge_list_edge_graph = toEdgeGraph(edge_list);
    int n = edge_list.size();
    int m = SETWORDSNEEDED(n);

    EMPTYGRAPH(g, m, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (edge_list_edge_graph[i][j] == 1) {
                ADDONEEDGE(g, i, j, m);
            }
        }
    }

    densenauty(g, lab, ptn, orbits, options, stats, m, n, NULL);

    int* result = new int[n];
    std::copy(orbits, orbits + n, result);

    return result;
}

// Conversion methods
std::string GroupGraph::printGraph() const {
    std::ostringstream output;
    output << "Nodes:\n";
    for (const auto& entry : nodes) {
        output << "    Group " << entry.first << " (" << entry.second.ntype << ") (" << entry.second.smarts << ") ";
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
               << std::get<2>(edge) << "(" << std::get<3>(edge) << ")  Order: " << std::get<4>(edge) << "\n";
    }
    return output.str();
}

std::string GroupGraph::toSmiles() const {
    std::unique_ptr<RDKit::RWMol> molecularGraph(new RDKit::RWMol());
    std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
    int atomCount = 0;

    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        std::string smarts = entry.second.smarts;
        std::unique_ptr<RDKit::RWMol> subGraph(RDKit::SmartsToMol(smarts));
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
        const Group& node = entry.second;
        nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)] = std::unordered_map<int, int>();
        std::string smarts = node.smarts;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmartsToMol(smarts));
        for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*atom)->getIdx()] = atomId;
        }
    }


    atomId = -1;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        std::string smarts = node.smarts;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmartsToMol(smarts));
        for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            RDKit::Atom newAtom = **atom;
            molecularGraph->addAtom(&newAtom, true);
        }
        for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
            RDKit::Bond newBond = **bond;
            RDKit::Bond::BondType bondType = (*bond)->getBondType();
            int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getBeginAtomIdx()];
            int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getEndAtomIdx()];
            molecularGraph->addBond(atomIdx1, atomIdx2, bondType);
        }
    }


    for (const auto& edge : edges) {
        NodeIDType from = std::get<0>(edge);
        PortType fromPort = std::get<1>(edge);
        NodeIDType to = std::get<2>(edge);
        PortType toPort = std::get<3>(edge);
        unsigned int bondOrder = std::get<4>(edge);
        int fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
        int toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];

        molecularGraph->addBond(fromAtom, toAtom, static_cast<RDKit::Bond::BondType>(bondOrder));
    }

    return RDKit::MolToSmiles(*molecularGraph, true);
}

std::vector<std::vector<int>> GroupGraph::toEdgeGraph(const std::vector<std::pair<int, int>>& edge_list) const {
    int num_edges = edge_list.size();

    // Initialize edge graph as an adjacency matrix where each element is initially 0
    std::vector<std::vector<int>> edge_graph(num_edges, std::vector<int>(num_edges, 0));

    // Iterate through each pair of edges in the edge list
    for (int i = 0; i < num_edges; i++) {
        for (int j = i + 1; j < num_edges; j++) {
            const std::pair<int, int>& edge_i = edge_list[i];
            const std::pair<int, int>& edge_j = edge_list[j];

            // Check if the two edges share a common node
            if (edge_i.first == edge_j.first || edge_i.first == edge_j.second ||
                edge_i.second == edge_j.first || edge_i.second == edge_j.second) {
                // If they share a node, mark them as connected in the edge graph
                edge_graph[i][j] = 1;
                edge_graph[j][i] = 1;  // Symmetric adjacency matrix
            }
        }
    }

    return edge_graph;
}

std::unordered_map<std::string, int> GroupGraph::toVector() const {
    std::unordered_map<std::string, int> hist(nodetypes.size());
    for (const auto& entry : nodes) {
        const Group& node = entry.second;
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
        const Group& node = entry.second;
        std::string smarts = entry.second.smarts;
        std::unique_ptr<RDKit::RWMol> subGraph(RDKit::SmartsToMol(smarts));
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
        const Group& node = entry.second;
        nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)] = std::unordered_map<int, int>();
        std::string smarts = node.smarts;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmartsToMol(smarts));
        for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*atom)->getIdx()] = atomId;
        }
    }

    // Start building the atom graph from mappings defined above
    atomId = -1;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        std::string smarts = node.smarts;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmartsToMol(smarts));
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
            unsigned int bondOrder = (*bond)->getBondTypeAsDouble();
            atomGraph->addEdge(atomIdx1, atomIdx2, bondOrder);
        }
    }


    for (const auto& edge : edges) {
        NodeIDType from = std::get<0>(edge);
        PortType fromPort = std::get<1>(edge);
        NodeIDType to = std::get<2>(edge);
        PortType toPort = std::get<3>(edge);
        NodeIDType fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
        NodeIDType toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];
        unsigned int bondOrder = std::get<4>(edge);
        atomGraph->addEdge(fromAtom, toAtom, bondOrder);
    }


    return atomGraph;
}

std::string GroupGraph::serialize() const {
        std::ostringstream oss;
        oss << "{\n  \"nodes\": [\n";
        for (const auto& pair : nodes) {
            const Group& node = pair.second;
            oss << "    {\n      \"id\": " << node.id
                << ",\n      \"ntype\": \"" << node.ntype
                << "\",\n      \"smarts\": \"" << node.smarts
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

        // Step 1: Sort nodes based on their attributes (e.g., id, ntype, smarts)
        std::vector<Group> sortedNodes;
        for (const auto& pair : nodes) {
            sortedNodes.push_back(pair.second);
        }
        std::sort(sortedNodes.begin(), sortedNodes.end(), [](const Group& a, const Group& b) {
            return std::tie(a.ntype, a.smarts, a.id) < std::tie(b.ntype, b.smarts, b.id);
        });

        // Step 2: Sort edges by connected nodes and ports
        std::vector<std::tuple<NodeIDType, PortType, NodeIDType, PortType, unsigned int>> sortedEdges = edges;
        std::sort(sortedEdges.begin(), sortedEdges.end());

        // Step 3: Convert sorted nodes and edges to a canonical smarts or other format
        std::stringstream ss;
        for (const auto& node : sortedNodes) {
            ss << node.ntype << ":" << node.smarts << ";";
        }
        for (const auto& edge : sortedEdges) {
            ss << std::get<0>(edge) << "-" << std::get<1>(edge) << "-"
               << std::get<2>(edge) << "-" << std::get<3>(edge) << ":" << std::get<4>(edge) << ";";
        }

        // Step 4: Return the canonical representation
        return ss.str();
    }

//#############################################################################################################
//#############################################################################################################
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
    // Check if the number of nodes and edges are the same
    if (this->nodes.size() != other.nodes.size()) {
        return false;
    }
    int edgeCount1 = 0;
    int edgeCount2 = 0;
    for (const auto& node : this->edges) {
        edgeCount1 += node.second.size();
    }
    for (const auto& node : other.edges) {
        edgeCount2 += node.second.size();
    }
    if (edgeCount1 != edgeCount2) {
        return false;
    }

    // Convert AtomGraph to nauty graph
    int n = this->nodes.size(); // Assuming the number of nodes is the same for both graphs
    int m = SETWORDSNEEDED(n);

    // Use std::vector instead of DYNALLSTAT and DYNALLOC
    std::vector<setword> g1(m * n, 0); // Initialize graph 1
    std::vector<setword> g2(m * n, 0); // Initialize graph 2
    std::vector<int> lab1(n), ptn1(n), orbits1(n); // Label, partition, and orbits for graph 1
    std::vector<int> lab2(n), ptn2(n), orbits2(n); // Label, partition, and orbits for graph 2
    std::vector<setword> canong1(m * n, 0); // Canonical form for graph 1
    std::vector<setword> canong2(m * n, 0); // Canonical form for graph 2
    setword workspace[160]; // Workspace for nauty

    // Initialize nauty structures
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    // Convert AtomGraph to nauty graph for g1
    EMPTYGRAPH(g1.data(), m, n);
    for (const auto& [id, dst_order] : this->edges) {
        if (!dst_order.empty()) {
            for (const auto& [dest,order] : dst_order) {
                int from = id;
                int to = dest;
                ADDONEEDGE(g1.data(), from, to, m);
            }
        }
    }

    // Convert AtomGraph to nauty graph for g2
    EMPTYGRAPH(g2.data(), m, n);
    for (const auto& [id, dst_order] : other.edges) {
        if (!dst_order.empty()) {
            for (const auto& [dest,order] : dst_order) {
                int from = id;
                int to = dest;
                ADDONEEDGE(g2.data(), from, to, m);
            }
        }
    }

    // Call nauty to canonicalize the graphs
    options.getcanon = TRUE;
    nauty(g1.data(), lab1.data(), ptn1.data(), nullptr, orbits1.data(), &options, &stats, workspace, 160, m, n, canong1.data());
    nauty(g2.data(), lab2.data(), ptn2.data(), nullptr, orbits2.data(), &options, &stats, workspace, 160, m, n, canong2.data());

    // Compare the canonical forms to determine isomorphism
    if (memcmp(canong1.data(), canong2.data(), sizeof(setword) * m * n) != 0) {
        return false;
    }

    return true;
}

void AtomGraph::addNode(const std::string& type, const unsigned int valency) {
    int id = nodes.size();
    nodes[id] = Atom(id, type, valency);

}

void AtomGraph::addEdge(NodeIDType src, NodeIDType dst, unsigned int order) {
    if (nodes.find(src) == nodes.end() || nodes.find(dst) == nodes.end()) {
        if (nodes.find(src) == nodes.end()) {
            throw std::invalid_argument("Group " + std::to_string(src) + " does not exist");
        }
        else {
            throw std::invalid_argument("Group " + std::to_string(dst) + " does not exist");
        }
    }
    if (getFreeValency(src) <= 0 && getFreeValency(dst) <= 0) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for both nodes");
    }
    if (getFreeValency(src) <= 0) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for the source node");
    }
    if (getFreeValency(dst) <= 0) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for the destination node");
    }
    if (edges[src].find(std::make_pair(dst, order)) != edges[src].end() || edges[dst].find(std::make_pair(src, order)) != edges[dst].end()) {
        throw std::invalid_argument("Edge from " + std::to_string(src) + " to " + std::to_string(dst) + " already exists");
    }
    if (order > 4 || order < 1) {
        throw std::invalid_argument("Bond order of " + std::to_string(order) + " is invalid");
    }
    edges[src].insert(std::make_pair(dst, order));
    edges[dst].insert(std::make_pair(src, order));
}

std::vector<std::vector<std::pair<AtomGraph::NodeIDType, AtomGraph::NodeIDType>>> AtomGraph::substructureSearch(const AtomGraph& query, const std::vector<int>& hubs) const {
    /*
        Returns a list of all subgraph isomorphisms between the query graph and this graph
        format is a list of lists of pairs of node ids where (query_node_id, this_node_id) is a match
    */
    std::vector<std::vector<std::pair<NodeIDType,NodeIDType>>> matches; // To store all matches
    std::unordered_map<NodeIDType, int> queryNeededFreeValency; // To store the number of hubs for each query node

    // Step 0: Pre-process query hubs
    for (const auto& node : query.nodes) {
        queryNeededFreeValency[node.first] = 0;
    }
    for (const auto& h : hubs) {
        queryNeededFreeValency[h]++;
    }
    // Add a count for the bonds that are in the mol graph
    // for (const auto& [nodeid, dstSet] : query.edges) {
    //     int totalBondCount = 0;
    //     for (const auto& [dst, order] : dstSet) {
    //         totalBondCount += order;
    //     }
    //     queryNeededFreeValency[nodeid] += totalBondCount;
    // }
    // std::cout<<"Query Hub Counts: "<<std::endl;
    // for (const auto& [id, count] : queryNeededFreeValency) {
    //     std::cout << id << ": " << count << std::endl;
    // }
    // std::cout << std::endl;

    // std::cout<<"Query Graph: "<<std::endl;
    // std::cout<<query.printGraph()<<std::endl;
    // std::cout<<"This Graph: "<<std::endl;
    // std::cout<<this->printGraph()<<std::endl;

    // std::cout<<"Query Hub Counts: "<<std::endl;
    // for (const auto& [id, count] : queryNeededFreeValency) {
    //     std::cout << id << ": " << count << std::endl;
    // }
    // std::cout << std::endl;

    // Step 1: Pre-filter nodes in the graph based on query node attributes
    std::unordered_map<NodeIDType, std::vector<NodeIDType>> candidateNodes; // Maps query nodes to possible candidates in the main graph
    for (const auto& queryNodePair : query.nodes) {
        const auto& queryNode = queryNodePair.second;
        for (const auto& graphNodePair : nodes) {
            const auto& graphNode = graphNodePair.second;

            // Match based on node type and valency
            if (queryNode.ntype == graphNode.ntype && queryNode.valency <= graphNode.valency) {
                candidateNodes[queryNode.id].push_back(graphNode.id);
            }
        }
    }

    // std::cout << "Candidate Nodes: " << std::endl;
    // for (const auto& [queryNode, candidates] : candidateNodes) {
    //     std::cout << queryNode << ": ";
    //     for (const auto& candidate : candidates) {
    //         std::cout << candidate << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // Step 2: Backtracking function to explore mappings
    std::function<void(std::unordered_map<NodeIDType, NodeIDType>&, std::unordered_set<NodeIDType>&)> backtrack =
        [&](std::unordered_map<NodeIDType, NodeIDType>& currentMapping, std::unordered_set<NodeIDType>& usedNodes) {
            // Debug: Print the current mapping
            // std::cout << "Current Mapping: ";
            // for (const auto& mapping : currentMapping) {
            //     std::cout << "(" << mapping.first << " -> " << mapping.second << ") ";
            // }
            // std::cout << std::endl;
            // If all query nodes are mapped, validate hubs
            if (currentMapping.size() == query.nodes.size()) {
                // Check if the hubs specified match the query node hubs
                for (const auto& [id, count] : queryNeededFreeValency) {
                    NodeIDType graphNodeid = currentMapping[id];
                    // printf("Query Group: %d", id);
                    // printf("Graph Group: %d", graphNodeid);
                    // printf("this->getFreeValency(graphNodeid): %d", this->getFreeValency(graphNodeid));
                    // printf("count: %d", count);
                    // printf("\n");
                    if(this->getFreeValency(graphNodeid) != count) { // Check if number of bonds for query node matches the number of hubs
                        return;
                    }
                }
                // printf("Hubs Matched");
                // Check if the bonds are the same
                for (const auto& [queryNodeid, dstSet] : query.edges) {
                    for (const auto& [dst, order] : dstSet) {
                        NodeIDType graphNodeid = currentMapping[queryNodeid];
                        auto it = edges.find(graphNodeid);
                        if (it == edges.end() || it->second.find(std::make_pair(currentMapping[dst], order)) == it->second.end()) {
                            return;
                        }
                    }
                }
                // printf("Bonds Matched\n");


                // Add the valid match
                std::vector<std::pair<NodeIDType, NodeIDType>> match;
                for (const auto& mapping : currentMapping) {
                    match.push_back(std::make_pair(mapping.first, mapping.second));
                }
                matches.push_back(match);
                return;
            }

            // Select the next unmapped query node
            NodeIDType nextQueryNode = -1;
            for (const auto& queryNodePair : query.nodes) {
                if (currentMapping.find(queryNodePair.first) == currentMapping.end()) {
                    nextQueryNode = queryNodePair.first;
                    break;
                }
            }

            if (nextQueryNode == -1) return; // No unmapped node found (shouldn't happen)

            // Try each candidate node for the selected query node
            for (NodeIDType candidate : candidateNodes[nextQueryNode]) {
                if (usedNodes.count(candidate)) continue;

                // Check edge consistency
                bool valid = true;
                for (const auto& [queryNeighbor, graphNeighbor] : currentMapping) {
                    // Check if the query graph has an edge between nextQueryNode and queryNeighbor
                    bool edgeExists = false;
                    for (const auto& edgePair : query.edges.at(nextQueryNode)) {
                        if (edgePair.first == queryNeighbor) {
                            edgeExists = true;
                            break;
                        }
                    }

                    if (edgeExists) {
                        // Verify the corresponding edge exists in the main graph
                        bool graphEdgeExists = false;
                        for (const auto& edgePair : edges.at(candidate)) {
                            if (edgePair.first == graphNeighbor) {
                                graphEdgeExists = true;
                                break;
                            }
                        }
                        if (!graphEdgeExists) {
                            valid = false;
                            break;
                        }
                    }
                }

                if (!valid) continue;


                // Temporarily map the query node to the candidate
                currentMapping[nextQueryNode] = candidate;
                usedNodes.insert(candidate);

                // Recurse
                backtrack(currentMapping, usedNodes);

                // Backtrack
                currentMapping.erase(nextQueryNode);
                usedNodes.erase(candidate);
            }
        };

    // Step 3: Initialize and start the backtracking process
    std::unordered_map<NodeIDType, NodeIDType> currentMapping; // Maps query node IDs to graph node IDs
    std::unordered_set<NodeIDType> usedNodes; // Tracks already used graph nodes
    backtrack(currentMapping, usedNodes);

    return matches;
}

void AtomGraph::fromSmiles(const std::string& smiles) {
    nodes.clear();
    edges.clear();

    std::unordered_map<std::string, int> standardElementValency = {
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"P", 3}, {"S", 2}, {"Cl", 1}, {"Br", 1}, {"I", 1}
    };
    
    std::stack<NodeIDType> nodeStack; // Stack to handle branching
    std::unordered_map<int, NodeIDType> ringClosures; // Map for ring closure indices
    NodeIDType lastNode = -1;
    int bondOrder = 1; // Default to single bond

    for (size_t i = 0; i < smiles.size(); ++i) {
        char c = smiles[i];

        if (std::isalpha(c)) {
            // Handle atom
            int valency = standardElementValency[std::string(1, c)];
            addNode(std::string(1, c), valency);
            NodeIDType currentNode = nodes.size() - 1;

            // If there's a previous node, add an edge with the current bond order
            if (lastNode != -1) {
                addEdge(lastNode, currentNode, bondOrder);
            }

            lastNode = currentNode;
            bondOrder = 1; // Reset bond order to single after use
        } else if (c == '(') {
            // Start a branch, push the current last node onto the stack
            nodeStack.push(lastNode);
        } else if (c == ')') {
            // End a branch, pop the last node from the stack
            if (!nodeStack.empty()) {
                lastNode = nodeStack.top();
                nodeStack.pop();
            } else {
                throw std::runtime_error("Unmatched closing parenthesis in SMILES string.");
            }
        } else if (std::isdigit(c)) {
            // Handle ring closure
            int ringIndex = c - '0';
            if (ringClosures.count(ringIndex)) {
                // Connect the current node to the ring closure with the current bond order
                addEdge(lastNode, ringClosures[ringIndex], bondOrder);
                ringClosures.erase(ringIndex);
            } else {
                // Store the current node as the ring closure point
                ringClosures[ringIndex] = lastNode;
            }
            bondOrder = 1; // Reset bond order to single after use
        } else if (c == '=') {
            // Set bond order to double
            bondOrder = 2;
        } else {
            // Handle unsupported characters (e.g., invalid SMILES)
            throw std::invalid_argument("Unsupported character in SMILES: " + std::string(1, c));
        }
    }

    // Basic error checking for unclosed rings
    if (!ringClosures.empty()) {
        std::cerr << "Unclosed rings detected: ";
        for (const auto& entry : ringClosures) {
            std::cerr << entry.first << " ";
        }
        std::cerr << std::endl;
        throw std::runtime_error("Unclosed ring detected in SMILES string.");
    }

}


int AtomGraph::getFreeValency(NodeIDType nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Cannot get free valency for non-existent node " + std::to_string(nodeID));
    }
    const Atom& node = nodes.at(nodeID);
    if (edges.find(nodeID) == edges.end()) {
        return node.valency;
    }
    else{
        int occupied_electrons = 0;
        for (const auto& edge : edges.at(nodeID)) {
            occupied_electrons += std::get<1>(edge);
        }
        return node.valency - occupied_electrons;
    }
}

std::string AtomGraph::printGraph() const {
    std::ostringstream output;
    output << "Nodes:\n";
    for (const auto& entry : nodes) {
        output << "    Atom " << entry.first << " (" << entry.second.ntype << ")" << " Valency: " << entry.second.valency << "\n";
    }
    output << "Edges:\n";
    for (const auto& edge : edges) {
        for (const auto& [dst, order] : edge.second) {
            output << "    Edge: " << edge.first << " -> " << dst <<" Order: (" <<order<<")"<<"\n";
        }
    }
    return output.str();
}

std::vector<setword> AtomGraph::toNautyGraph() const {
    // Convert AtomGraph to a nauty graph representation
    int n = nodes.size(); // Number of nodes
    int m = SETWORDSNEEDED(n); // Size of one row of the adjacency matrix in setwords

    // Allocate storage for the graph
    std::vector<setword> g(m * n, 0); // Initialize nauty graph (adjacency matrix)

    // Initialize the nauty graph
    EMPTYGRAPH(g.data(), m, n);

    // Add edges to the graph
    for (const auto& [id, dst_order] : edges) {
        for (const auto& [dest, order] : dst_order) {
            ADDONEEDGE(g.data(), id, dest, m); // Add edge from 'id' to 'dest'
        }
    }

    // Return the nauty graph representation
    return g;
}


std::vector<std::vector<AtomGraph::NodeIDType>> AtomGraph::nodeAut() const {
    int n = nodes.size(); // Number of nodes
    int m = SETWORDSNEEDED(n); // Size of one row of the adjacency matrix in setwords

    // Prepare vectors and workspace
    std::vector<int> lab(n), ptn(n), orbits(n); // Label, partition, and orbits
    std::vector<setword> canong(m * n, 0);      // Canonical form
    setword workspace[160];                    // Workspace for nauty

    // Convert the AtomGraph to a Nauty graph representation
    std::vector<setword> g = this->toNautyGraph();

    // Set up Nauty options
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = true; // Calculate the canonical labeling

    statsblk stats; // Statistics block

    // Run Nauty
    nauty(g.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options, &stats, workspace, 160, m, n, canong.data());

    // Process the automorphism results
    // Orbits represent the equivalence classes of automorphism; map them to a vector of vectors
    std::vector<std::vector<AtomGraph::NodeIDType>> automorphisms(n);
    for (int i = 0; i < n; ++i) {
        automorphisms[orbits[i]].push_back(i);
    }

    // Remove empty entries from automorphisms
    automorphisms.erase(std::remove_if(automorphisms.begin(), automorphisms.end(),
                                       [](const std::vector<int>& v) { return v.empty(); }),
                        automorphisms.end());

    return automorphisms;
}

std::vector<AtomGraph::NodeIDType> AtomGraph::nodeOrbits() const {
    int n = nodes.size(); // Number of nodes
    int m = SETWORDSNEEDED(n); // Size of one row of the adjacency matrix in setwords

    // Prepare vectors and workspace
    std::vector<int> lab(n), ptn(n), orbits(n); // Label, partition, and orbits
    std::vector<setword> canong(m * n, 0);      // Canonical form
    setword workspace[160];                    // Workspace for nauty

    // Convert the AtomGraph to a Nauty graph representation
    std::vector<setword> g = this->toNautyGraph();

    // Set up Nauty options
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = true; // Calculate the canonical labeling

    statsblk stats; // Statistics block

    // Run Nauty
    nauty(g.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options, &stats, workspace, 160, m, n, canong.data());

    return orbits;
}


// std::string AtomGraph::toSmiles(
// ) const {
// }
