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


class AtomGraph {
public:
    using NodeIDType = int;

    struct Node {
        NodeIDType id;
        std::string ntype;
        unsigned int valency;

        // Need to define comparison operators for Node to be used in unordered_map
        bool operator==(const Node& other) const {
            return id == other.id && ntype == other.ntype && valency == other.valency;
        }
        bool operator!=(const Node& other) const {
            return !(*this == other);
        }
        Node() : id(0), ntype(""), valency(0) {}

        Node(int id, const std::string& ntype, const unsigned int valency)
            : id(id), ntype(ntype), valency(valency) {}

    };

    // Custom hash function for the Node struct
    struct NodeHasher {
        std::size_t operator()(const Node& node) const {
            return std::hash<NodeIDType>()(node.id) ^ std::hash<std::string>()(node.ntype) ^ std::hash<unsigned int>()(node.valency);
        }
    };

    std::unordered_map<NodeIDType, Node> nodes; ///< Map of node IDs to their respective nodes.
    std::unordered_map<NodeIDType, std::unordered_set<NodeIDType>> edges; ///< Map of node IDs to sets of connected node IDs.

    AtomGraph();
    AtomGraph(const AtomGraph& other);
    AtomGraph& operator=(const AtomGraph& other);
    bool operator==(const AtomGraph& other) const;

    void addNode(const std::string& ntype = "", unsigned int valency = 0);
    void addEdge(NodeIDType src, NodeIDType dst);
    int getFreeValency(NodeIDType nid) const;
    std::string printGraph() const;

private:
    // Helper methods or additional private members can be declared here if needed
};


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
    std::unique_ptr<AtomGraph> toAtomicGraph() const;

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

#endif // DATASTRUCTURES_H
