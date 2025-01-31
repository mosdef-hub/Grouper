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

// std::unordered_map<GroupGraph::Group, std::unordered_set<GroupGraph::Group>> determineNodeComposition(const std::unordered_map<std::string, GroupGraph::Group>& nodeDefs){
//     std::unordered_map<GroupGraph::Group, std::unordered_set<GroupGraph::Group>> nodeComposition;
//     for (const auto& nodeDef1 : nodeDefs) {
//         const GroupGraph::Group& nodeData = nodeDef1.second;
//         const std::string& smiles = nodeData.smarts; // Need to use smiles because rdkit can't actually do substructure matching with smarts

//         for (const auto& nodeDef2 : nodeDefs) {
//             const std::string& smarts2 = nodeDef2.first;
//             const GroupGraph::Group& nodeData2 = nodeDef2.second;
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
std::unordered_set<GroupGraph::Group> possibleValencyNode(const GroupGraph::Group& node) {
    /*
        Given a node, return all possible nodes that can be created by adding sampling the hubs
        The upper bound is 2^hubs.size() - 1 but because the atoms in a node can be automorphic to each other this is in practice much smaller
        Ex. C=C, [0,0,1,1] -> C=C, [0,0,1,1], C=C, [0,0,1], C=C, [0,1], C=C, [0,0], C=C, [0], C=C []
        TODO: Implement a more efficient algorithm, currently this is generates all possible subsets, would be more efficient to incorporate the automorphisms
    */

    // Resulting set of nodes
    std::unordered_set<GroupGraph::Group> possibleNodes;

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
        GroupGraph::Group newNode = node;
        newNode.hubs = subset;
        newNode.ports.resize(subset.size());
        std::iota(newNode.ports.begin(), newNode.ports.end(), 0);

        possibleNodes.insert(newNode);
    }

    return possibleNodes;

}
// TODO: This should accept nodeDefs which are NodeTraces
// TODO: This should handle Groups which can be SMARTS or SMILES
std::vector<GroupGraph> fragment(
    const std::string& smiles,
    const std::unordered_set<GroupGraph::Group>& nodeDefs
) {
    std::vector<GroupGraph> allFragmentations; // Store all valid fragmentations
    GroupGraph groupGraph;
    std::vector<GroupGraph::Group> finalNodeDefinitions;

    // Parse the SMILES string to an AtomGraph
    AtomGraph mol;
    mol.fromSmiles(smiles);

    // Calculate possible nodes for each node definition
    std::unordered_map<GroupGraph::Group, std::unordered_set<GroupGraph::Group>> possibleNodes;
    for (const auto& nodeDef : nodeDefs) {
        possibleNodes[nodeDef] = possibleValencyNode(nodeDef);
    }

    // Convert to sorted vectors
    std::unordered_map<GroupGraph::Group, std::vector<GroupGraph::Group>> possibleNodesVec;
    for (const auto& [node, possibleNodeSet] : possibleNodes) {
        possibleNodesVec[node] = std::vector<GroupGraph::Group>(possibleNodeSet.begin(), possibleNodeSet.end());
    }

    std::vector<GroupGraph::Group> nodesVec;
    for (const auto& [node, possibleNodeSet] : possibleNodes) {
        nodesVec.push_back(node);
    }

    std::sort(nodesVec.begin(), nodesVec.end(),
        [](const GroupGraph::Group& a, const GroupGraph::Group& b) {
            AtomGraph aGraph, bGraph;
            aGraph.fromSmiles(a.smarts); // fromSmarts?
            bGraph.fromSmiles(b.smarts); // fromSmarts?
            if (aGraph.nodes.size() != bGraph.nodes.size()) {
                return aGraph.nodes.size() > bGraph.nodes.size();
            }
            return a.hubs.size() > b.hubs.size();
        }
    );

    for (auto& [node, nodeCompositionVec] : possibleNodesVec) {
        std::sort(nodeCompositionVec.begin(), nodeCompositionVec.end(),
        [](const GroupGraph::Group& a, const GroupGraph::Group& b) {
            AtomGraph aGraph, bGraph;

            aGraph.fromSmiles(a.smarts); // fromSmarts
            bGraph.fromSmiles(b.smarts); // fromSmarts
            if (aGraph.nodes.size() != bGraph.nodes.size()) {
                return aGraph.nodes.size() > bGraph.nodes.size();
            }
            return a.hubs.size() > b.hubs.size();
        });
    }

    // Flatten the possible nodes
    std::vector<GroupGraph::Group> sortedNodes;
    for (const auto& node : nodesVec) {
        for (const auto& possibleNode : possibleNodesVec[node]) {
            sortedNodes.push_back(possibleNode);
        }
    }

    // Recursive function to explore all valid fragmentations
    std::function<void(const std::vector<GroupGraph::Group>&,
                       std::unordered_map<int, GroupGraph::NodeIDType>,
                       std::unordered_map<int, std::vector<int>>,
                       std::unordered_map<int, std::string>,
                       int,
                       GroupGraph)>
    attemptFragmentation = [&](const std::vector<GroupGraph::Group>& candidates,
                               std::unordered_map<int, GroupGraph::NodeIDType> atomToNodeid,
                               std::unordered_map<int, std::vector<int>> atomToPorts,
                               std::unordered_map<int, std::string> atomToSmarts,
                               int startIdx,
                               GroupGraph currentGraph) {
        if (atomToNodeid.size() == mol.nodes.size()) {
            // Add edges to the group graph based on molecular connectivity
            for (const auto& [src, dst, bondOrder] : mol.edges) {
                GroupGraph::NodeIDType srcNode = atomToNodeid[src];
                GroupGraph::NodeIDType dstNode = atomToNodeid[dst];
                if (srcNode == dstNode) continue;

                // Check if the edge already exists
                bool edgeExists = false;
                for (const auto& edge : currentGraph.edges) {
                    if ((std::get<0>(edge) == srcNode && std::get<2>(edge) == dstNode) || (std::get<0>(edge) == dstNode && std::get<2>(edge) == srcNode)) {
                        edgeExists = true;
                        break;
                    }
                }
                if (edgeExists) continue;
                // Find available ports
                int srcPort = -1, dstPort = -1;
                for (int port : atomToPorts[src]) {
                    if (currentGraph.isPortFree(srcNode, port)) {
                        srcPort = port;
                        break;
                    }
                }
                for (int port : atomToPorts[dst]) {
                    if (currentGraph.isPortFree(dstNode, port)) {
                        dstPort = port;
                        break;
                    }
                }

                if (srcPort == -1 || dstPort == -1) {
                    // printf("Failed to find free ports for edge: %d-%d\n%s\n", src, dst, currentGraph.printGraph().c_str());
                    return;
                    // throw std::invalid_argument("No free ports available for " + std::string(atomToSmarts[src]) +  " ( " + std::to_string(src) + " )"+ " -> " + std::string(atomToSmarts[dst]) + " ( " + std::to_string(dst) + " )");
                }

                    // printf("GroupGraph : %s\n", currentGraph.printGraph().c_str());
                    // printf("AtomtoPorts:\n");
                    // for (const auto& [atom, ports] : atomToPorts) {
                    //     printf("%d: ", atom);
                    //     printf ("Ports: ");
                    //     for (const auto& port : ports) {
                    //         printf("%d ", port);
                    //     }
                    //     printf("Smarts: %s ", atomToSmarts[atom].c_str());
                    //     printf(" element %s", mol.nodes[atom].ntype.c_str());
                    //     printf("NodeID: %d\n", atomToNodeid[atom]);
                    //     printf("\n");
                    // }
                    // printf("Edge: %d-%d, Ports: %d-%d\n", src, dst, srcPort, dstPort);
                    // printf("Smarts: %s-%s\n", atomToSmarts[src].c_str(), atomToSmarts[dst].c_str());

                    // printf("Bond Order: %d\n", bondOrder);

                currentGraph.addEdge(
                    std::make_tuple(srcNode, srcPort),
                    std::make_tuple(dstNode, dstPort),
                    bondOrder
                );
            }
            allFragmentations.push_back(currentGraph); // Store this valid fragmentation
            // printf("Added fragmentation\n");
            return;
        }

        for (size_t i = startIdx; i < candidates.size(); ++i) {
            const auto& node = candidates[i];
            AtomGraph query;
            query.fromSmiles(node.smarts); // fromSMARTS

            auto matches = mol.substructureSearch(query, node.hubs);
            if (matches.empty()) {
                continue;
            }

            std::unordered_map<int, GroupGraph::NodeIDType> tempAtomToNodeid = atomToNodeid;
            std::unordered_map<int, std::vector<int>> tempAtomToPorts = atomToPorts;
            std::unordered_map<int, std::string> tempAtomToSmarts = atomToSmarts;
            GroupGraph tempGraph = currentGraph;  // Copy to maintain independent states
            GroupGraph::NodeIDType currentId = tempGraph.nodes.size();

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
                    [&](const GroupGraph::Group& n) {
                        return n.smarts == node.smarts && n.ntype == node.ntype;
                    }
                );

                if (parent == nodeDefs.end()) {
                    throw std::invalid_argument("Parent node not found.");
                }

                tempGraph.addNode(parent->ntype, parent->smarts, parent->hubs);

                for (const auto& [queryid, molid] : match) {
                    tempAtomToNodeid[molid] = currentId;
                    tempAtomToSmarts[molid] = node.smarts;

                    std::vector<int> hubIndices;
                    for (size_t idx = 0; idx < parent->hubs.size(); ++idx) {
                        if (parent->hubs[idx] == queryid) {
                            hubIndices.push_back(idx);
                        }
                    }

                    tempAtomToPorts[molid] = hubIndices;
                }
                currentId++;

                // Recursive call to explore further
                attemptFragmentation(candidates, tempAtomToNodeid, tempAtomToPorts, tempAtomToSmarts, i + 1, tempGraph);
            }
        }
    };

    // Initial attempt to fragment using the sorted nodes
    std::unordered_map<int, GroupGraph::NodeIDType> atomToNodeid;
    std::unordered_map<int, std::vector<int>> atomToPorts;
    std::unordered_map<int, std::string> atomToSmarts;

    attemptFragmentation(sortedNodes, atomToNodeid, atomToPorts, atomToSmarts, 0, groupGraph);

    if (allFragmentations.empty()) {
        throw std::invalid_argument("Failed to fully fragment the molecule.");
    }

    return allFragmentations;
}
