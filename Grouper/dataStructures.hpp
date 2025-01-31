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

    struct Atom {
        NodeIDType id;
        std::string ntype;
        unsigned int valency;

        // Need to define comparison operators for Atom to be used in unordered_map
        bool operator==(const Atom& other) const {
            return id == other.id && ntype == other.ntype && valency == other.valency;
        }
        bool operator!=(const Atom& other) const {
            return !(*this == other);
        }
        Atom() : id(0), ntype(""), valency(0) {}

        Atom(int id, const std::string& ntype, const unsigned int valency)
            : id(id), ntype(ntype), valency(valency) {}

        std::string toString() const;

    };

    // Custom hash function for the Atom struct
    struct AtomHasher {
        std::size_t operator()(const Atom& node) const {
            return std::hash<NodeIDType>()(node.id) ^ std::hash<std::string>()(node.ntype) ^ std::hash<unsigned int>()(node.valency);
        }
    };

    // Tuple hasher
    struct PairHasher {
        template <class T1, class T2>
        std::size_t operator()(const std::pair<T1, T2>& t) const {
            return std::hash<T1>{}(t.first) ^ std::hash<T2>{}(t.second);
        }
    };

    std::unordered_map<NodeIDType, Atom> nodes; ///< Map of node IDs to their respective nodes.
    std::unordered_map<NodeIDType, std::unordered_set<std::pair<NodeIDType, unsigned int>, PairHasher>> edges; ///< Map of node IDs to their respective neighbors.

    AtomGraph();
    AtomGraph(const AtomGraph& other);
    AtomGraph& operator=(const AtomGraph& other);
    bool operator==(const AtomGraph& other) const;

    void addNode(const std::string& ntype = "", unsigned int valency = 0);
    void addEdge(NodeIDType src, NodeIDType dst, unsigned int order = 1);
    int getFreeValency(NodeIDType nid) const;
    std::string printGraph() const;
    std::vector<std::vector<NodeIDType>> nodeAut() const;
    std::vector<NodeIDType> nodeOrbits() const;
    std::vector<setword> toNautyGraph() const;
    void fromSmiles(const std::string& smiles);
    std::vector<std::vector<std::pair<AtomGraph::NodeIDType,AtomGraph::NodeIDType>>> substructureSearch(const AtomGraph& query, const std::vector<int>& hubs) const;

private:
};

class GroupGraph {
public:
    using NodeIDType = int;
    using PortType = int;

    // Attributes
    struct Group {
        NodeIDType id;
        std::string ntype;
        std::string smarts;
        std::vector<NodeIDType> hubs;
        std::vector<PortType> ports;
        bool operator==(const Group& other) const;
        bool operator!=(const Group& other) const;
        Group() : ntype(""), smarts(""), hubs(), ports() {}
        Group(const std::string& ntype, const std::string& smarts, const std::vector<int>& hubs)
            : ntype(ntype), smarts(smarts), hubs(hubs), ports(hubs.size()) {
            std::iota(ports.begin(), ports.end(), 0);
        }

        std::vector<int> hubOrbits() const;
        std::string toString() const;

    };
    std::unordered_map<NodeIDType, Group> nodes; ///< Map of node IDs to their respective nodes.
    std::vector<std::tuple<NodeIDType, PortType, NodeIDType, PortType, unsigned int>> edges; ///< List of edges connecting nodes. (srcNodeID, srcPort, dstNodeID, dstPort, bondOrder)
    std::unordered_map<std::string, std::vector<PortType>> nodetypes; ///< Map of node types to their respective ports.

    // Core Methods
    GroupGraph();
    GroupGraph(const GroupGraph& other);
    GroupGraph& operator=(const GroupGraph& other);
    bool operator==(const GroupGraph& other) const;
    // Operating methods
    void addNode(
        std::string ntype,
        std::string smarts,
        std::vector<NodeIDType> hubs
    );
    bool addEdge(
        std::tuple<NodeIDType, PortType> fromNodePort,
        std::tuple<NodeIDType, PortType> toNodePort,
        unsigned int bondOrder = 1,
        bool verbose = false
    );
    int numFreePorts(NodeIDType nid) const;
    int* computeEdgeOrbits(
        const std::vector<std::pair<int, int>> edge_list,
        graph* g, int* lab, int* ptn, int* orbits,
        optionblk* options, statsblk* stats
        ) const;
    void clearEdges();
    bool isPortFree(NodeIDType nodeID, PortType port) const;
    // Conversion methods
    std::string printGraph() const;
    std::unordered_map<std::string, int> toVector() const;
    std::string toSmiles() const;
    std::unique_ptr<AtomGraph> toAtomicGraph() const;
    std::string serialize() const;
    void deserialize(const std::string& state);
    std::string canonize() const;


private:
    std::vector<std::vector<int>> toEdgeGraph(const std::vector<std::pair<int, int>>& edge_list) const;
    void toNautyGraph(int *n, int *m, int *adj) const;
    int numNodes() const;
};

inline bool operator<(const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>& lhs,
                      const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>& rhs);

// Specialize std::hash for GroupGraph::Group, GroupGraph, GroupGraph.edge, AtomGraph::Atom
namespace std {
    template <>
    struct hash<GroupGraph::Group> {
        std::size_t operator()(const GroupGraph::Group& node) const {
            std::size_t h1 = std::hash<int>{}(node.id);
            std::size_t h2 = std::hash<std::string>{}(node.ntype);
            std::size_t h3 = std::hash<std::string>{}(node.smarts);
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
    struct hash<AtomGraph::Atom> {
        std::size_t operator()(const AtomGraph::Atom& atom) const {
            std::size_t h1 = std::hash<int>{}(atom.id);
            std::size_t h2 = std::hash<std::string>{}(atom.ntype);
            std::size_t h3 = std::hash<unsigned int>{}(atom.valency);

            // Combine the individual hashes using XOR and shifting
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
    template <>
    struct hash<std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, 
                        GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>> {
        std::size_t operator()(const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, 
                                                GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>& t) const {
            std::size_t h1 = std::hash<GroupGraph::NodeIDType>{}(std::get<0>(t));
            std::size_t h2 = std::hash<GroupGraph::PortType>{}(std::get<1>(t));
            std::size_t h3 = std::hash<GroupGraph::NodeIDType>{}(std::get<2>(t));
            std::size_t h4 = std::hash<GroupGraph::PortType>{}(std::get<3>(t));
            std::size_t h5 = std::hash<unsigned int>{}(std::get<4>(t));

            return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4);
        }
    };
    template <>
    struct hash<GroupGraph> {
        std::size_t operator()(const GroupGraph& graph) const {
            std::size_t h = 0;
            for (const auto& node_pair : graph.nodes) {
                const auto& node = node_pair.second;
                h ^= std::hash<GroupGraph::Group>{}(node) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            for (const auto& edge : graph.edges) {
                h ^= std::hash<std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, 
                                        GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>>{}(edge) 
                    + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
    template <>
    struct hash<AtomGraph> {
        std::size_t operator()(const AtomGraph& graph) const {
            std::size_t h = 0;

            // Hash the nodes (unordered_map<NodeIDType, Atom>)
            for (const auto& node_pair : graph.nodes) {
                const auto& node = node_pair.second;
                h ^= std::hash<AtomGraph::Atom>{}(node) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }

            // Hash the edges (unordered_map<NodeIDType, unordered_set<std::pair<NodeIDType, unsigned int>>)
            for (const auto& edge_pair : graph.edges) {
                const auto& edge_set = edge_pair.second;
                for (const auto& edge : edge_set) {
                    // Custom hash for std::pair<NodeIDType, unsigned int>
                    std::size_t h1 = std::hash<AtomGraph::NodeIDType>{}(edge.first);
                    std::size_t h2 = std::hash<unsigned int>{}(edge.second);

                    // Combine the two hashes
                    h ^= h1 ^ (h2 << 1) + 0x9e3779b9 + (h << 6) + (h >> 2);
                }
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
