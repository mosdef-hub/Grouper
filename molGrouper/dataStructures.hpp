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
        bool operator==(const Node& other) const;
        Node() : id(0), ntype(""), smiles(""), ports(), hubs() {}

        Node(int id, const std::string& ntype, const std::string& smiles, const std::vector<int>& ports, const std::vector<int>& hubs)
            : id(id), ntype(ntype), smiles(smiles), ports(ports), hubs(hubs) {}
    };
    std::unordered_map<NodeIDType, Node> nodes; ///< Map of node IDs to their respective nodes.
    std::vector<std::tuple<NodeIDType, PortType, NodeIDType, PortType>> edges; ///< List of edges connecting nodes.
    std::unordered_map<std::string, std::vector<PortType>> nodetypes; ///< Map of node types to their respective ports.

    // Methods
    GroupGraph();
    GroupGraph(const GroupGraph& other);
    GroupGraph& operator=(const GroupGraph& other);
    bool operator==(const GroupGraph& other) const;

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
    std::string printGraph() const;
    std::unordered_map<std::string, int> toVector() const;
    std::string toSmiles() const;

private:
    // Helper methods or additional private members can be declared here if needed
};

inline bool operator<(const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType>& lhs,
                      const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType>& rhs);

// Specialize std::hash for GroupGraph::Node
namespace std {
    template <>
    struct hash<GroupGraph::Node> {
        std::size_t operator()(const GroupGraph::Node& node) const {
            std::size_t h1 = std::hash<int>{}(node.id);
            std::size_t h2 = std::hash<std::string>{}(node.ntype);
            std::size_t h3 = std::hash<std::string>{}(node.smiles);
            std::size_t h4 = 0;
            for (const auto& port : node.ports) {
                h4 ^= std::hash<int>{}(port) + 0x9e3779b9 + (h4 << 6) + (h4 >> 2);
            }
            std::size_t h5 = 0;
            for (const auto& hub : node.hubs) {
                h5 ^= std::hash<int>{}(hub) + 0x9e3779b9 + (h5 << 6) + (h5 >> 2);
            }
            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4);
        }
    };
}
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
