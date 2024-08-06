#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>

#include "GroupGraph.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>


GroupGraph::GroupGraph(const std::unordered_map<std::string, std::vector<PortType>>& types)
    : nodeTypes(types) {}

GroupGraph::GroupGraph(const GroupGraph& other)
    : nodes(other.nodes), edges(other.edges), nodeTypes(other.nodeTypes) {}

GroupGraph& GroupGraph::operator=(const GroupGraph& other) {
    if (this != &other) {
        nodes = other.nodes;
        edges = other.edges;
        nodeTypes = other.nodeTypes;
    }
    return *this;
}

void GroupGraph::addNode(NodeID id, const std::string& type) {
    if (nodes.find(id) != nodes.end()) {
        throw std::invalid_argument("Node already exists");
    }
    nodes[id] = {type, nodeTypes.at(type)};
}

bool GroupGraph::addEdge(NodeID from, PortType fromPort, NodeID to, PortType toPort, bool verbose) {
    if (n_free_ports(from) <= 0) {
        if (verbose) {
            std::cout << "Source node doesn't have enough ports!" << std::endl;
        }
        return false;
    }
    if (n_free_ports(to) <= 0) {
        if (verbose) {
            std::cout << "Destination node doesn't have enough ports!" << std::endl;
        }
        return false;
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

int GroupGraph::n_free_ports(NodeID nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Node does not exist");
    }
    const Node& node = nodes.at(nodeID);
    int occupied_ports = 0;
    for (const auto& edge : edges) {
        if (std::get<0>(edge) == nodeID) {
            occupied_ports++;
        }
    }
    return node.ports.size() - occupied_ports;
}

int GroupGraph::numNodes() const {
    return nodes.size();
}

void GroupGraph::printGraph() const {
    std::cout << "Nodes:\n";
    for (const auto& entry : nodes) {
        std::cout << "Node " << entry.first << " (" << entry.second.type << ") : Ports ";
        for (PortType port : entry.second.ports) {
            std::cout << port << " ";
        }
        std::cout << "\n";
    }

    std::cout << "Edges:\n";
    for (const auto& edge : edges) {
        std::cout << "Edge: " << std::get<0>(edge) << "(" << std::get<1>(edge) << ") -> " << std::get<2>(edge) << "(" << std::get<3>(edge) << ")\n";
    }
}

std::string GroupGraph::toSmiles(
    const std::unordered_map<std::string, std::string>& nodeTypeToSmiles,
    const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex
) const {
    std::unique_ptr<RDKit::RWMol> molecularGraph(new RDKit::RWMol());
    std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
    int atomCount = 0;

    for (const auto& entry : nodes) {
        NodeID nodeId = entry.first;
        const Node& node = entry.second;
        std::string smiles = nodeTypeToSmiles.at(node.type);
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
        nodePortToAtomIndex[std::to_string(nodeId)] = std::unordered_map<int, int>();
        for (size_t i = 0; i < node.ports.size(); ++i) {
            nodePortToAtomIndex[std::to_string(nodeId)][node.ports[i]] = atomCount + nodeTypePortToIndex.at(node.type).at(node.ports[i]);
        }
        atomCount += subGraph->getNumAtoms();
    }

    int atomId = -1;
    std::unordered_map<std::string, std::unordered_map<int, int>> nodeSubGraphIndicesToMolecularGraphIndices;
    for (const auto& entry : nodes) {
        NodeID nodeId = entry.first;
        const Node& node = entry.second;
        nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)] = std::unordered_map<int, int>();
        std::string smiles = nodeTypeToSmiles.at(node.type);
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
        for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)][(*atom)->getIdx()] = atomId;
        }
    }

    atomId = -1;
    for (const auto& entry : nodes) {
        NodeID nodeId = entry.first;
        const Node& node = entry.second;
        std::string smiles = nodeTypeToSmiles.at(node.type);
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
        for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            RDKit::Atom newAtom = **atom; 
            molecularGraph->addAtom(&newAtom, true);
        }
        for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
            RDKit::Bond newBond = **bond;
            int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)][(*bond)->getBeginAtomIdx()];
            int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)][(*bond)->getEndAtomIdx()];
            molecularGraph->addBond(atomIdx1, atomIdx2, newBond.getBondType());
        }
    }

    for (const auto& edge : edges) {
        NodeID from = std::get<0>(edge);
        PortType fromPort = std::get<1>(edge);
        NodeID to = std::get<2>(edge);
        PortType toPort = std::get<3>(edge);
        int fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
        int toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];
        molecularGraph->addBond(fromAtom, toAtom, RDKit::Bond::SINGLE);
    }

    std::string smiles = RDKit::MolToSmiles(*molecularGraph);
    
    return smiles;
}



// int main() {

//     std::unordered_map<std::string, std::vector<int>> node_types = {
//         {"N", {0}},
//         {"CO", {0, 1}},
//         {"CC", {0, 1, 2, 3}}
//     };
//     std::unordered_map<int, std::string> int_to_node_type = {
//         {0, "N"},
//         {1, "CO"},
//         {2, "CC"},
//         {3, "CC"},
//         {4, "CC"}
//     };
//     const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex = {
//         {"N", {{0, 0}}},
//         {"CO", {{0, 0}, {1, 0}}},
//         {"CC", {{0, 0}, {1, 0}, {2, 1}, {3, 1}}}
//     };
//     const std::unordered_map<std::string, std::string>& nodeTypeToSmiles = {
//         {"N", "N"},
//         {"CO", "CO"},
//         {"CC", "C=C"}
//     };

//     GroupGraph gG(node_types);
//     gG.addNode(0, "N");
//     gG.addNode(1, "CC");

//     gG.addEdge(0, 0, 1, 0);

//     // gG.printGraph();

//     RDKit::ROMol* moleGraph = gG.toMolecularGraph(nodeTypeToSmiles, nodeTypePortToIndex);

//     std::string smiles = RDKit::MolToSmiles(*moleGraph);
//     std::cout << "SMILES: " << smiles << std::endl;

//     return 0;
// }
