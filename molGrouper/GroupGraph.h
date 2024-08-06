#ifndef GROUPGRAPH_H
#define GROUPGRAPH_H

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
    using NodeID = int;
    using PortType = int;

    struct Node {
        std::string type;
        std::vector<PortType> ports;
    };

    std::unordered_map<NodeID, Node> nodes; ///< Map of node IDs to their respective nodes.
    std::vector<std::tuple<NodeID, PortType, NodeID, PortType>> edges; ///< List of edges connecting nodes.
    std::unordered_map<std::string, std::vector<PortType>> nodeTypes; ///< Map of node types to their respective ports.

    GroupGraph(const std::unordered_map<std::string, std::vector<PortType>>& types);
    GroupGraph(const GroupGraph& other);
    GroupGraph& operator=(const GroupGraph& other);

    void addNode(NodeID id, const std::string& type);
    bool addEdge(NodeID from, PortType fromPort, NodeID to, PortType toPort, bool verbose = false);
    int n_free_ports(NodeID nodeID) const;
    int numNodes() const;
    void printGraph() const;

    std::string toSmiles(
        const std::unordered_map<std::string, std::string>& nodeTypeToSmiles,
        const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex
    ) const;

private:
    // Helper methods or additional private members can be declared here if needed
};

#endif // GROUPGRAPH_H
