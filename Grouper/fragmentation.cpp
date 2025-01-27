#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <nauty/nauty.h>

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <iostream>  // For debugging output


#include "dataStructures.hpp"


struct TupleHash {
    template <typename T1, typename T2, typename T3>
    std::size_t operator ()(const std::tuple<T1, T2, T3>& t) const {
        auto h1 = std::hash<T1>{}(std::get<0>(t)); 
        auto h2 = std::hash<T2>{}(std::get<1>(t)); 
        auto h3 = std::hash<T3>{}(std::get<2>(t)); 

        // Combine the hashes using a better mix (commonly used strategy)
        h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h2
        h1 ^= h3 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h3
        return h1;
    }
};

namespace std {
    template <>
    struct hash<std::tuple<char, char, int>> {
        std::size_t operator()(const std::tuple<char, char, int>& t) const {
            auto h1 = std::hash<char>{}(std::get<0>(t));
            auto h2 = std::hash<char>{}(std::get<1>(t));
            auto h3 = std::hash<int>{}(std::get<2>(t));

            // Combine the hashes using a better mix (commonly used strategy)
            h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h2
            h1 ^= h3 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h3
            return h1;
        }
    };
}

// std::unordered_map<GroupGraph::Node, std::unordered_set<GroupGraph::Node>> determineNodeComposition(const std::unordered_map<std::string, GroupGraph::Node>& nodeDefs){
//     std::unordered_map<GroupGraph::Node, std::unordered_set<GroupGraph::Node>> nodeComposition;
//     for (const auto& nodeDef1 : nodeDefs) {
//         const GroupGraph::Node& nodeData = nodeDef1.second;
//         const std::string& smiles = nodeData.smarts; // Need to use smiles because rdkit can't actually do substructure matching with smarts

//         for (const auto& nodeDef2 : nodeDefs) {
//             const std::string& smarts2 = nodeDef2.first;
//             const GroupGraph::Node& nodeData2 = nodeDef2.second;
//             std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts2));
//             std::vector<RDKit::MatchVectType> matches;
//             RDKit::SubstructMatch(*original, *query, matches);
//             if (!matches.empty()) {
//                 nodeComposition[nodeData].insert(nodeData2);
//             }
//         }
//         if (nodeComposition[nodeData].empty()) {
//             nodeComposition[nodeData].insert(nodeData);
//         }
//     }
//     return nodeComposition;
// }
std::unordered_set<GroupGraph::Node> possibleValencyNode(const GroupGraph::Node& node) {
    /*
        Given a node, return all possible nodes that can be created by adding sampling the hubs
        The upper bound is 2^hubs.size() - 1 but because the atoms in a node can be automorphic to each other this is in practice much smaller
        Ex. C=C, [0,0,1,1] -> C=C, [0,0,1,1], C=C, [0,0,1], C=C, [0,1], C=C, [0,0], C=C, [0], C=C []
        TODO: Implement a more efficient algorithm, currently this is generates all possible subsets, would be more efficient to incorporate the automorphisms
    */ 

    // Resulting set of nodes
    std::unordered_set<GroupGraph::Node> possibleNodes;

    // Get the total number of hubs
    size_t hubCount = node.hubs.size();

    // Generate all subsets of hubs (2^hubCount subsets)
    size_t subsetCount = 1 << hubCount; // 2^hubCount

    for (size_t mask = 0; mask < subsetCount; ++mask) {
        // Create a subset of hubs based on the current mask
        std::vector<int> subset;
        for (size_t i = 0; i < hubCount; ++i) {
            if (mask & (1 << i)) {
                subset.push_back(node.hubs[i]);
            }
        }

        // Create a new node with this subset of hubs
        GroupGraph::Node newNode = node;
        newNode.hubs = subset;
        newNode.ports.resize(subset.size());
        std::iota(newNode.ports.begin(), newNode.ports.end(), 0);

        possibleNodes.insert(newNode);
    }

    return possibleNodes;

}

GroupGraph fragment(
    const std::string& smiles, 
    const std::unordered_set<GroupGraph::Node>& nodeDefs
) {
    GroupGraph groupGraph;

    // Parse the SMILES string to an AtomGraph
    AtomGraph mol;
    mol.fromSmiles(smiles);
    std::cout << mol.printGraph();

    // Calculate possible nodes for each node definition
    std::unordered_map<GroupGraph::Node, std::unordered_set<GroupGraph::Node>> possibleNodes;
    for (const auto& nodeDef : nodeDefs) {
        possibleNodes[nodeDef] = possibleValencyNode(nodeDef);
    }
    // Convert to a vector for sorting
    std::unordered_map<GroupGraph::Node, std::vector<GroupGraph::Node>> possibleNodesVec;
    for (const auto& [node, possibleNodeSet] : possibleNodes) {
        possibleNodesVec[node] = std::vector<GroupGraph::Node>(possibleNodeSet.begin(), possibleNodeSet.end());
    }
    // Sort the nodes by the number of atoms and hubs
    std::vector<GroupGraph::Node> nodesVec;
    for (const auto& [node, possibleNodeSet] : possibleNodes) {
        nodesVec.push_back(node);
    }
    std::sort(nodesVec.begin(), nodesVec.end(), 
        [](const GroupGraph::Node& a, const GroupGraph::Node& b) {
            AtomGraph aGraph;
            aGraph.fromSmiles(a.smarts);
            AtomGraph bGraph;
            bGraph.fromSmiles(b.smarts);
            if (aGraph.nodes.size() != bGraph.nodes.size()) {
                return aGraph.nodes.size() > bGraph.nodes.size();
            }
            return a.hubs.size() > b.hubs.size();
            }
        );
    // Sort the possible nodes by the number of atoms and hubs
    for (auto& [node, nodeCompositionVec] : possibleNodesVec) {
        std::sort(nodeCompositionVec.begin(), nodeCompositionVec.end(), 
        [](const GroupGraph::Node& a, const GroupGraph::Node& b) {
            AtomGraph aGraph;
            aGraph.fromSmiles(a.smarts);
            AtomGraph bGraph;
            bGraph.fromSmiles(b.smarts);
            if (aGraph.nodes.size() != bGraph.nodes.size()) {
                return aGraph.nodes.size() > bGraph.nodes.size();
            }
            return a.hubs.size() > b.hubs.size();
            }
        );
    }

    // flatten the possible nodes
    std::vector<GroupGraph::Node> sortedNodes;
    for (const auto& node : nodesVec) {
        for (const auto& possibleNode : possibleNodesVec[node]) {
            sortedNodes.push_back(possibleNode);
        }
    }
    printf("Sorted possible nodes\n");

    // Recursive function for tree-like fragmentation
    std::function<bool(const std::vector<GroupGraph::Node>&, 
                std::unordered_map<int, GroupGraph::NodeIDType>&, 
                std::unordered_map<int, std::vector<int>>&, 
                std::unordered_map<int, std::string>&, 
                int)>
    attemptFragmentation = [&](const std::vector<GroupGraph::Node>& candidates, 
                                std::unordered_map<int, GroupGraph::NodeIDType>& atomToNodeid, 
                                std::unordered_map<int, std::vector<int>>& atomToPorts, 
                                std::unordered_map<int, std::string>& atomToSmarts, 
                                int startIdx) -> bool {
        if (atomToNodeid.size() == mol.nodes.size()) {
            return true; // Fully fragmented
        }

        for (size_t i = startIdx; i < candidates.size(); ++i) {
            const auto& node = candidates[i];
            AtomGraph query;
            query.fromSmiles(node.smarts);

            // Substructure search
            auto matches = mol.substructureSearch(query, node.hubs);
            if (matches.empty()) {
                continue;
            }

            // Attempt to apply this node
            std::unordered_map<int, GroupGraph::NodeIDType> tempAtomToNodeid = atomToNodeid;
            std::unordered_map<int, std::vector<int>> tempAtomToPorts = atomToPorts;
            std::unordered_map<int, std::string> tempAtomToSmarts = atomToSmarts;
            GroupGraph::NodeIDType currentId = groupGraph.nodes.size();

            bool applied = false;
            for (const auto& match : matches) {
                bool alreadyMatched = false;
                for (const auto& [queryid, molid] : match) {
                    if (tempAtomToNodeid.count(molid)) {
                        alreadyMatched = true;
                        break;
                    }
                }
                if (alreadyMatched) continue;

                // Apply the match
                // // Get the parent node
                auto parent = std::find_if(nodeDefs.begin(), nodeDefs.end(), 
                    [&](const GroupGraph::Node& n) {
                        return n.smarts == node.smarts && n.ntype == node.ntype;
                    }
                );
                if (parent == nodeDefs.end()) {
                    throw std::invalid_argument("Parent node not found.");
                }
                groupGraph.addNode(parent->ntype, parent->smarts, parent->hubs);
                applied = true;
                for (const auto& [queryid, molid] : match) {
                    tempAtomToNodeid[molid] = currentId;
                    tempAtomToSmarts[molid] = node.smarts;

                    // Map ports
                    std::vector<int> hubIndices;
                    for (size_t idx = 0; idx < node.hubs.size(); ++idx) {
                        if (node.hubs[idx] == queryid) {
                            hubIndices.push_back(idx);
                        }
                    }
                    tempAtomToPorts[molid] = hubIndices;
                }
                break;
            }

            if (applied) {
                // Recur with all nodes, including smaller and larger nodes
                if (attemptFragmentation(candidates, tempAtomToNodeid, tempAtomToPorts, tempAtomToSmarts, 0)) {
                    return true;
                }

                // Backtrack
                groupGraph.nodes.erase(groupGraph.nodes.size() - 1);
            }
        }

        return false;
    };

    // Attempt to fragment using the sorted nodes
    std::unordered_map<int, GroupGraph::NodeIDType> atomToNodeid;
    std::unordered_map<int, std::vector<int>> atomToPorts;
    std::unordered_map<int, std::string> atomToSmarts;

    if (!attemptFragmentation(sortedNodes, atomToNodeid, atomToPorts, atomToSmarts, 0)) {
        throw std::invalid_argument("Failed to fully fragment the molecule.");
    }

    printf("Fragmentation successful\n");
    printf("Group graph:\n");
    printf("%s\n", groupGraph.printGraph().c_str());

    // Add edges to the group graph based on molecular connectivity
    for (const auto& [src, dstSet] : mol.edges) {
        for (const auto& [dst, bondOrder] : dstSet) {
            GroupGraph::NodeIDType srcNode = atomToNodeid[src];
            GroupGraph::NodeIDType dstNode = atomToNodeid[dst];
            if (srcNode == dstNode) continue;

            // Find available ports
            int srcPort = -1, dstPort = -1;
            for (int port : atomToPorts[src]) {
                if (groupGraph.isPortFree(srcNode, port)) {
                    srcPort = port;
                    break;
                }
            }
            for (int port : atomToPorts[dst]) {
                if (groupGraph.isPortFree(dstNode, port)) {
                    dstPort = port;
                    break;
                }
            }

            if (srcPort == -1 || dstPort == -1) {
                throw std::invalid_argument("No free ports available for edge.");
            }

            printf("Adding edge: %d-%d\n", srcNode, dstNode);
            printf("Ports: %d-%d\n", srcPort, dstPort);

            groupGraph.addEdge(
                std::make_tuple(srcNode, srcPort),
                std::make_tuple(dstNode, dstPort)
            );
        }
    }

    return groupGraph;
}

