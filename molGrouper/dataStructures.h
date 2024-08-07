#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <memory>
#include <tuple>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

class GroupGraph {
public:
    using NodeIDType = int;
    using PortType = int;

    // Attributes
    struct Node {
        NodeIDType id;
        std::string ntype;
        std::string smiles;
        std::vector<PortType> ports;
        std::vector<NodeIDType> hubs;
    };
    std::unordered_map<NodeIDType, Node> nodes; ///< Map of node IDs to their respective nodes.
    std::vector<std::tuple<NodeIDType, PortType, NodeIDType, PortType>> edges; ///< List of edges connecting nodes.
    std::unordered_map<std::string, std::vector<PortType>> nodetypes; ///< Map of node types to their respective ports.

    // Methods
    GroupGraph();
    GroupGraph(const GroupGraph& other);
    GroupGraph& operator=(const GroupGraph& other);
    void addNode(
        std::string ntype, 
        std::string smiles, 
        std::vector<PortType> ports,
        std::vector<NodeIDType> hubs
    );
    bool addEdge(
        std::tuple<NodeIDType, PortType> fromNodePort,
        std::tuple<NodeIDType, PortType> toNodePort,
        bool verbose = false
    );
    int n_free_ports(NodeIDType nid) const;
    int numNodes() const;
    void printGraph() const;
    std::string toSmiles() const;
    // std::vector<int> toVector() const;

private:
    // Helper methods or additional private members can be declared here if needed
};

// class AtomGraph {
// public:
//     using NodeIDType = int;

//     struct Node {
//         NodeIDType id;
//         std::string ntype;
//     };

//     std::unordered_map<NodeIDType, Node> nodes; ///< Map of node IDs to their respective nodes.
//     std::vector<std::tuple<NodeIDType, NodeIDType>> edges; ///< List of edges connecting nodes.

//     AtomGraph();
//     AtomGraph(const AtomGraph& other);
//     AtomGraph& operator=(const AtomGraph& other);

//     void addNode(
//         const std::string& ntype = ""
//     );
//     bool addEdge(
//         NodeIDType fromNode,
//         NodeIDType toNode,
//         bool verbose = false
//     );
//     int n_free_ports(NodeIDType nid) const;
//     int numNodes() const;
//     void printGraph() const;

//     std::string toSmiles(
//         const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex
//     ) const;

// private:
//     // Helper methods or additional private members can be declared here if needed
// };

#endif // DATASTRUCTURES_H
