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


#include <nauty/nauty.h>
#include <nauty/naututil.h>

#include "autUtils.hpp"


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
    std::vector<std::vector<int>> Aut() const; // Returns the automorphism group of the graph
    void toNautyFormat(int *n, int *m, int *adj) const;

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
        std::vector<NodeIDType> hubs;
        std::vector<PortType> ports;
        bool operator==(const Node& other) const;
        Node() : id(0), ntype(""), smiles(""), hubs(), ports() {}
        Node(int id, const std::string& ntype, const std::string& smiles, const std::vector<int>& hubs)
            : id(id), ntype(ntype), smiles(smiles), hubs(hubs), ports(hubs.size()) {
            std::iota(ports.begin(), ports.end(), 0);
        }

    };
    std::unordered_map<NodeIDType, Node> nodes; ///< Map of node IDs to their respective nodes.
    std::vector<std::tuple<NodeIDType, PortType, NodeIDType, PortType>> edges; ///< List of edges connecting nodes.
    std::unordered_map<std::string, std::vector<PortType>> nodetypes; ///< Map of node types to their respective ports.

    // Core Methods
    GroupGraph();
    GroupGraph(const GroupGraph& other);
    GroupGraph& operator=(const GroupGraph& other);
    bool operator==(const GroupGraph& other) const;
    // Operating methods
    void addNode(
        std::string ntype, 
        std::string smiles, 
        std::vector<NodeIDType> hubs
    );
    bool addEdge(
        std::tuple<NodeIDType, PortType> fromNodePort,
        std::tuple<NodeIDType, PortType> toNodePort,
        bool verbose = false
    );
    int n_free_ports(NodeIDType nid) const;
    int* computeEdgeOrbits(
        const std::vector<std::pair<int, int>> edge_list,
        graph* g, int* lab, int* ptn, int* orbits,
        optionblk* options, statsblk* stats
        ) const;
    void clearEdges();
    // Conversion methods
    std::string printGraph() const;
    std::unordered_map<std::string, int> toVector() const;
    std::string toSmiles() const;
    std::unique_ptr<AtomGraph> toAtomicGraph() const;
    std::string serialize() const;
    std::string Canon() const;
    

private:
    std::vector<std::vector<int>> toEdgeGraph(const std::vector<std::pair<int, int>>& edge_list) const;
    void toNautyGraph(int *n, int *m, int *adj) const;
    int numNodes() const;
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
    template <>
    struct hash<std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType>> {
        std::size_t operator()(const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType>& t) const {
            std::size_t h1 = std::hash<GroupGraph::NodeIDType>{}(std::get<0>(t));
            std::size_t h2 = std::hash<GroupGraph::PortType>{}(std::get<1>(t));
            std::size_t h3 = std::hash<GroupGraph::NodeIDType>{}(std::get<2>(t));
            std::size_t h4 = std::hash<GroupGraph::PortType>{}(std::get<3>(t));
            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
        }
    };
    template <>
    struct hash<GroupGraph> {
        std::size_t operator()(const GroupGraph& graph) const {
            std::size_t h = 0;
            for (const auto& node_pair : graph.nodes) {
                const auto& node = node_pair.second;
                h ^= std::hash<GroupGraph::Node>{}(node) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            for (const auto& edge : graph.edges) {
                h ^= std::hash<std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType>>{}(edge) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
    template <>
    struct hash<std::pair<int, int>> {
        std::size_t operator()(const std::pair<int, int>& p) const {
            std::size_t h1 = std::hash<int>{}(p.first);
            std::size_t h2 = std::hash<int>{}(p.second);
            return h1 ^ (h2 << 1); // or use another combination method
        }
    };
    template <>
    struct hash<std::vector<int>> {
        std::size_t operator()(const std::vector<int>& vec) const {
            std::size_t h = 0;
            for (const auto& elem : vec) {
                h ^= std::hash<int>{}(elem) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
}

#endif // DATASTRUCTURES_H
