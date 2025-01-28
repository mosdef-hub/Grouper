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
    std::vector<GroupGraph::Node> finalNodeDefinitions;

    // Parse the SMILES string to an AtomGraph
    AtomGraph mol;
    mol.fromSmiles(smiles);
    // std::cout << mol.printGraph();

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
    // printf("Sorted possible nodes\n");

    for(const auto& node : sortedNodes) {
        printf("Node sorted: %s\n", node.smarts.c_str());
    }

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

            printf("Node: %s\n", node.smarts.c_str());
            printf("Hubs: ");
            for (int hub : node.hubs) {
                printf("%d ", hub);
            }
            printf("\n");

            // Substructure search
            auto matches = mol.substructureSearch(query, node.hubs);
            if (matches.empty()) {
                continue;
            }

            std::unordered_map<int, GroupGraph::NodeIDType> tempAtomToNodeid = atomToNodeid;
            std::unordered_map<int, std::vector<int>> tempAtomToPorts = atomToPorts;
            std::unordered_map<int, std::string> tempAtomToSmarts = atomToSmarts;
            GroupGraph::NodeIDType currentId = groupGraph.nodes.size();

            // std::string smarts = node.smarts;

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

                auto parent = std::find_if(
                    nodeDefs.begin(), 
                    nodeDefs.end(),
                    [&](const GroupGraph::Node& n) {
                        return n.smarts == node.smarts && n.ntype == node.ntype;
                    }
                );

                if (parent == nodeDefs.end()) {
                    throw std::invalid_argument("Parent node not found.");
                }
                
                groupGraph.addNode(parent->ntype, parent->smarts, parent->hubs);
                applied = true;


                // Keep track of the added node definition
                finalNodeDefinitions.push_back(node);

                for (const auto& [queryid, molid] : match) {
                    tempAtomToNodeid[molid] = currentId;
                    tempAtomToSmarts[molid] = node.smarts;

                    std::vector<int> hubIndices;
                    for (size_t idx = 0; idx < node.hubs.size(); ++idx) {
                        if (node.hubs[idx] == queryid) {
                            hubIndices.push_back(idx);
                        }
                    }

                    if (!hubIndices.empty()) {
                        tempAtomToPorts[molid] = hubIndices;
                    } else {
                        tempAtomToPorts[molid] = std::vector<int>(parent->hubs.size());
                        std::iota(tempAtomToPorts[molid].begin(), tempAtomToPorts[molid].end(), 0);
                    }
                }
                break;
            }

            if (applied) {
                if (attemptFragmentation(candidates, tempAtomToNodeid, tempAtomToPorts, tempAtomToSmarts, 0)) {
                    atomToNodeid = std::move(tempAtomToNodeid);
                    atomToPorts = std::move(tempAtomToPorts);
                    atomToSmarts = std::move(tempAtomToSmarts);
                    return true;
                }

                groupGraph.nodes.erase(currentId);

                // Remove the last added node definition (backtracking)
                finalNodeDefinitions.pop_back();
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

    // Identify atom to ports
    atomToPorts.clear();
    atomToSmarts.clear();
    atomToNodeid.clear();
    GroupGraph::NodeIDType currentId = -1;
    for (const auto& node : finalNodeDefinitions) {
        // printf("Node: %s\n", node.smarts.c_str());
        AtomGraph query;
        query.fromSmiles(node.smarts);
        auto matches = mol.substructureSearch(query, node.hubs);
        if (matches.empty()) {
            throw std::invalid_argument("Failed to find a match for node definition.");
        }
        for (const auto& match : matches) {
            bool alreadyMatched = false;
            for (const auto& [queryid, molid] : match) {
                if (atomToNodeid.count(molid)) {
                    alreadyMatched = true;
                    break;
                }
            }
            if (alreadyMatched) continue;
            currentId += 1;
            for (const auto& queryid_thisid : match) {
                // printf("queryid: %d, thisid: %d\n", queryid_thisid.first, queryid_thisid.second);
                atomToNodeid[queryid_thisid.second] = currentId;
                atomToSmarts[queryid_thisid.second] = node.smarts;
                std::vector<int> hubIndices;
                auto parent = std::find_if(
                    nodeDefs.begin(), 
                    nodeDefs.end(),
                    [&](const GroupGraph::Node& n) {
                        return n.smarts == node.smarts && n.ntype == node.ntype;
                    }
                );
                for (size_t idx = 0; idx < parent->hubs.size(); ++idx) {
                    const auto& hub = parent->hubs[idx];
                    if (hub == queryid_thisid.first) {
                        hubIndices.push_back(idx);
                    }
                }
                atomToPorts[queryid_thisid.second] = hubIndices;
            }
        }
    }

    

    


    // printf("Fragmentation successful\n");
    // printf("size of atomToNodeid: %d\n", atomToNodeid.size());
    printf("Group graph:\n");
    printf("%s\n", groupGraph.printGraph().c_str());
    for (const auto& [atom, nodeid] : atomToNodeid) {
        printf("Atom %d Node %d", atom, nodeid);
        printf("    Ports: ");
        for (int port : atomToPorts[atom]) {
            printf("%d ", port);
        }
        printf("Smarts: %s", atomToSmarts[atom].c_str());
        printf("Element: %s\n", mol.nodes[atom].ntype.c_str());
    }

    // Add edges to the group graph based on molecular connectivity
    for (const auto& [src, dstSet] : mol.edges) {
        for (const auto& [dst, bondOrder] : dstSet) {
            GroupGraph::NodeIDType srcNode = atomToNodeid[src];
            GroupGraph::NodeIDType dstNode = atomToNodeid[dst];
            if (srcNode == dstNode) continue;

            // Check if the edge already exists
            bool edgeExists = false;
            for (const auto& edge : groupGraph.edges) {
                if ((std::get<0>(edge) == srcNode && std::get<2>(edge) == dstNode) || (std::get<0>(edge) == dstNode && std::get<2>(edge) == srcNode)) {
                    edgeExists = true;
                    break;
                }
            }
            if (edgeExists) continue;
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
                throw std::invalid_argument("No free ports available for " + std::string(atomToSmarts[src]) +  " (" + std::to_string(src) + " )"+ " -> " + std::string(atomToSmarts[dst]) + " (" + std::to_string(dst) + ")");
            }

            // printf("Adding edge: %d-%d\n", srcNode, dstNode);
            // printf("Ports: %d-%d\n", srcPort, dstPort);
            // printf("Smarts: %s-%s\n", atomToSmarts[src].c_str(), atomToSmarts[dst].c_str());

            groupGraph.addEdge(
                std::make_tuple(srcNode, srcPort),
                std::make_tuple(dstNode, dstPort)
            );
        }
    }

    return groupGraph;
}

