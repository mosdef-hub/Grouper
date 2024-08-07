#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>

#include "dataStructures.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>


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
    if (ports.size() == 0) {
        if (nodetypes.find(ntype) == nodetypes.end() && smiles == "") {
            throw std::invalid_argument("Node type or smiles haven't been defined");
        }
        if (nodetypes.find(smiles) == nodetypes.end() && ntype == "") {
            throw std::invalid_argument("Node type or smiles haven't been defined");
        }
    }
    if (hubs.size() == 0) {
        if (nodetypes.find(ntype) == nodetypes.end() && smiles == "") {
            throw std::invalid_argument("Node type or smiles haven't been defined");
        }
        if (nodetypes.find(smiles) == nodetypes.end() && ntype == "") {
            throw std::invalid_argument("Node type or smiles haven't been defined");
        }
    }
    if (smiles == "") {
        if (nodetypes.find(ntype) != nodetypes.end()) {
            if (nodetypes[ntype] != ports) {
                throw std::invalid_argument("Node type already exists with different ports");
            }
        } 
        if (nodetypes.find(ntype) == nodetypes.end()) {
            id = nodetypes.size();
            nodetypes[ntype] = ports;
        } 
    }
    if (ntype == "") {
        if (nodetypes.find(smiles) != nodetypes.end()) {
            if (nodetypes[smiles] != ports) {
                throw std::invalid_argument("Smiles already exists with different ports");
            }
        }
        else {
            nodetypes[smiles] = ports;
            ntype = smiles;
        }
        
    }
    id = nodes.size();
    nodes[id] = Node{id, ntype, smiles, ports, hubs};
    
}

bool GroupGraph::addEdge(std::tuple<NodeIDType,PortType> fromNodePort, std::tuple<NodeIDType,PortType>toNodePort, bool verbose) {
    NodeIDType from = std::get<0>(fromNodePort);
    PortType fromPort = std::get<1>(fromNodePort);
    NodeIDType to = std::get<0>(toNodePort);
    PortType toPort = std::get<1>(toNodePort);

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

int GroupGraph::n_free_ports(NodeIDType nodeID) const {
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
        std::cout << "NodeType " << entry.first << " (" << entry.second.ntype << ") : Ports ";
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

std::string GroupGraph::toSmiles() const {
    std::unique_ptr<RDKit::RWMol> molecularGraph(new RDKit::RWMol());
    std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
    int atomCount = 0;

    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Node& node = entry.second;
        std::string smiles = entry.second.smiles;
        std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
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

    std::string smiles = RDKit::MolToSmiles(*molecularGraph);
    
    return smiles;
}



// AtomGraph::AtomGraph()
//     : nodes(), edges() {}

// AtomGraph::AtomGraph(const AtomGraph& other)
//     : nodes(other.nodes), edges(other.edges) {}

// AtomGraph& AtomGraph::operator=(const AtomGraph& other) {
//     if (this != &other) {
//         nodes = other.nodes;
//         edges = other.edges;
//     }
//     return *this;
// }

// void AtomGraph::addNode(const std::string& type) {
//     nodes[nodes.size()] = {type};
// }

// bool AtomGraph::addEdge(NodeIDType from, PortType fromPort, NodeIDType to, PortType toPort, bool verbose) {
//     if (n_free_ports(from) <= 0) {
//         if (verbose) {
//             std::cout << "Source node doesn't have enough ports!" << std::endl;
//         }
//         return false;
//     }
//     if (n_free_ports(to) <= 0) {
//         if (verbose) {
//             std::cout << "Destination node doesn't have enough ports!" << std::endl;
//         }
//         return false;
//     }
//     if (nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) {
//         throw std::invalid_argument("Node does not exist");
//     }
//     if (std::find(nodes[from].ports.begin(), nodes[from].ports.end(), fromPort) == nodes[from].ports.end()) {
//         throw std::invalid_argument("Port does not exist");
//     }
//     edges.push_back(std::make_tuple(from, fromPort, to, toPort));
//     return true;
// }

// int AtomGraph::n_free_ports(NodeIDType nodeID) const {
//     if (nodes.find(nodeID) == nodes.end()) {
//         throw std::invalid_argument("Node does not exist");
//     }
//     const Node& node = nodes.at(nodeID);
//     int occupied_ports = 0;
//     for (const auto& edge : edges) {
//         if (std::get<0>(edge) == nodeID) {
//             occupied_ports++;
//         }
//     }
//     return node.ports.size() - occupied_ports;
// }

// int AtomGraph::numNodes() const {
//     return nodes.size();
// }

// void AtomGraph::printGraph() const {
//     std::cout << "Nodes:\n";
//     for (const auto& entry : nodes) {
//         std::cout << "Node " << entry.first << " (" << entry.second.type << ") : Ports ";
//         for (PortType port : entry.second.ports) {
//             std::cout << port << " ";
//         }
//         std::cout << "\n";
//     }

//     std::cout << "Edges:\n";
//     for (const auto& edge : edges) {
//         std::cout << "Edge: " << std::get<0>(edge) << "(" << std::get<1>(edge) << ") -> " << std::get<2>(edge) << "(" << std::get<3>(edge) << ")\n";
//     }
// }

// std::string AtomGraph::toSmiles(
// ) const {
//     std::unique_ptr<RDKit::RWMol> molecularGraph(new RDKit::RWMol());
//     std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
//     int atomCount = 0;

//     for (const auto& entry : nodes) {
//         NodeIDType nodeID = entry.first;
//         const Node& node = entry.second;
//         std::string smiles = nodeTypeToSmiles.at(node.type);
//         std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
//         nodePortToAtomIndex[std::to_string(nodeID)] = std::unordered_map<int, int>();
//         for (size_t i = 0; i < node.ports.size(); ++i) {
//             nodePortToAtomIndex[std::to_string(nodeID)][node.ports[i]] = atomCount + nodeTypePortToIndex.at(node.type).at(node.ports[i]);
//         }
//         atomCount += subGraph->getNumAtoms();
//     }

//     int atomId = -1;
//     std::unordered_map<std::string, std::unordered_map<int, int>> nodeSubGraphIndicesToMolecularGraphIndices;
//     for (const auto& entry : nodes) {
//         NodeIDType nodeID = entry.first;
//         const Node& node = entry.second;
//         nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)] = std::unordered_map<int, int>();
//         std::string smiles = nodeTypeToSmiles.at(node.type);
//         std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
//         for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
//             atomId++;
//             nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*atom)->getIdx()] = atomId;
//         }
//     }

//     atomId = -1;
//     for (const auto& entry : nodes) {
//         NodeIDType nodeID = entry.first;
//         const Node& node = entry.second;
//         std::string smiles = nodeTypeToSmiles.at(node.type);
//         std::unique_ptr<RDKit::ROMol> subGraph(RDKit::SmilesToMol(smiles));
//         for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
//             atomId++;
//             RDKit::Atom newAtom = **atom; 
//             molecularGraph->addAtom(&newAtom, true);
//         }
//         for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
//             RDKit::Bond newBond = **bond;
//             int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getBeginAtomIdx()];
//             int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getEndAtomIdx()];
//             molecularGraph->addBond(atomIdx1, atomIdx2, newBond.getBondType());
//         }
//     }

//     for (const auto& edge : edges) {
//         NodeIDType from = std::get<0>(edge);
//         PortType fromPort = std::get<1>(edge);
//         NodeIDType to = std::get<2>(edge);
//         PortType toPort = std::get<3>(edge);
//         int fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
//         int toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];
//         molecularGraph->addBond(fromAtom, toAtom, RDKit::Bond::SINGLE);
//     }

//     std::string smiles = RDKit::MolToSmiles(*molecularGraph);
    
//     return smiles;
// }
