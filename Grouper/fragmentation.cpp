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
    */


    std::unordered_set<GroupGraph::Node> possibleNodes;
    AtomGraph nodeGraph;
    nodeGraph.fromSmiles(node.smarts);

    // Get node orbits
    std::vector<int> orbits = nodeGraph.nodeOrbits();

    // Group hubs by their orbits
    std::unordered_map<int, std::vector<int>> orbitGroups;
    for (int hub : node.hubs) {
        orbitGroups[orbits[hub]].push_back(hub);
    }

    // Generate all combinations of hubs for each orbit
    std::vector<std::vector<std::vector<int>>> orbitCombinations;
    for (const auto& [orbit, hubs] : orbitGroups) {
        std::vector<std::vector<int>> combinations;

        // Generate subsets of hubs in the current orbit
        size_t subsetCount = hubs.size();
        for (size_t n = 1; n <= subsetCount; ++n) { // Avoid empty subset
            std::function<void(size_t, std::vector<int>)> generateSubset = 
                [&](size_t index, std::vector<int> currentSubset) {
                    if (currentSubset.size() == n) {
                        combinations.push_back(currentSubset);
                        return;
                    }
                    for (size_t i = index; i < hubs.size(); ++i) {
                        currentSubset.push_back(hubs[i]);
                        generateSubset(i + 1, currentSubset);
                        currentSubset.pop_back();
                    }
                };
            generateSubset(0, {});
        }

        orbitCombinations.push_back(combinations);
    }

    // Combine combinations across orbits
    std::vector<std::vector<int>> combinedHubs = {{}};
    for (const auto& combinations : orbitCombinations) {
        std::vector<std::vector<int>> newCombinations;

        for (const auto& existing : combinedHubs) {
            for (const auto& subset : combinations) {
                std::vector<int> newCombination = existing;
                newCombination.insert(newCombination.end(), subset.begin(), subset.end());
                newCombinations.push_back(newCombination);
            }
        }

        combinedHubs = std::move(newCombinations);
    }

    // Create nodes with the combined hub configurations
    for (const auto& hubs : combinedHubs) {
        GroupGraph::Node newNode = node;
        newNode.hubs = hubs;
        possibleNodes.insert(newNode);
    }

    // add the null node
    GroupGraph::Node nullNode = GroupGraph::Node(node.ntype, node.smarts, {});
    possibleNodes.insert(nullNode);

    return possibleNodes;
}

GroupGraph fragment(
    const std::string& smiles, 
    const std::unordered_set<GroupGraph::Node>& nodeDefs
) {

    GroupGraph groupGraph;

    printf("Got here\n");

    // Parse the SMILES string to an AtomGraph
    AtomGraph mol;
    mol.fromSmiles(smiles);

    std::cout<<mol.printGraph();

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
    // // print all possible nodes
    // for (const auto& [node, possibleNodes] : possibleNodes) {
    //     printf("Node: %s\n", node.smarts.c_str());
    //     for (const auto& possibleNode : possibleNodes) {
    //         printf("    %s\n", possibleNode.smarts.c_str());
    //         for (const auto& hub : possibleNode.hubs) {
    //             printf("        %d\n", hub);
    //         }
    //     }
    // }

    // Sort the possible nodes by the number of atoms and hubs
    for (auto& [node, nodeCompositionVec] : possibleNodesVec) {
        std::sort(nodeCompositionVec.begin(), nodeCompositionVec.end(), 
        [](const GroupGraph::Node& a, const GroupGraph::Node& b) {
            AtomGraph aGraph;
            aGraph.fromSmiles(a.smarts);
            AtomGraph bGraph;
            bGraph.fromSmiles(b.smarts);
            if (aGraph.nodes.size() != bGraph.nodes.size()) {
                return aGraph.nodes.size() < bGraph.nodes.size();
            }
            return a.hubs.size() < b.hubs.size();
            }
        );
    }
    printf("Sorted possible nodes\n");

    // Maps to keep track of atom indices and corresponding information
    std::unordered_map<int, GroupGraph::NodeIDType> atomToNodeid; // Maps atom indices to node IDs
    std::unordered_map<int, std::vector<int>> atomToPorts;  // atom to matching port indices
    std::unordered_map<int, std::string> atomToSmarts;  // Maps atom indices to SMARTS patterns
    GroupGraph::NodeIDType nodeId = 0;

    // Determine if nodes can be made by composing other nodes
    // std::unordered_map<GroupGraph::Node, std::unordered_set<GroupGraph::Node>> nodeComposition = determineNodeComposition(nodeDefs);
    for (const auto& [node, nodeCompositionsVec] : possibleNodesVec) {
        for (const auto& compNodeDef : nodeCompositionsVec) {
            const std::string& smarts = compNodeDef.smarts;
            AtomGraph query;
            query.fromSmiles(smarts);

            printf("Query: %s\n", smarts.c_str());
            
            std::vector<std::vector<std::pair<AtomGraph::NodeIDType,AtomGraph::NodeIDType>>> matches = mol.substructureSearch(query, compNodeDef.hubs);

            if (matches.empty()) {
                continue;
            }

            // print al; the matches idenitified
            for (const auto& match : matches) {
                printf("Match: ");
                for (const auto& [queryid, thisid] : match) {
                    printf("(%d, %d) ", queryid, thisid);
                }
                printf("\n");
            }

            for (const auto& match : matches) {
                // match is a vector of pairs (query atom index, mol atom index)
                bool alreadyMatched = false;
                for (const auto& queryid_thisid : match) {
                    if (atomToNodeid.count(queryid_thisid.second)) {
                        alreadyMatched = true;
                        break;
                    }
                }
                if (alreadyMatched) {
                    continue;
                }


                // Add the node to the graph
                groupGraph.addNode(
                    node.ntype, 
                    node.smarts, 
                    node.hubs
                );
                GroupGraph::NodeIDType currentId = nodeId++;

                for (const auto& queryid_thisid : match) {
                    atomToNodeid[queryid_thisid.second] = currentId;
                    atomToSmarts[queryid_thisid.second] = smarts;  // Track which SMARTS pattern the atom belongs to
                    // check if atom is associated with multiple ports
                    std::vector<int> hubIndices;
                    for (size_t idx = 0; idx < node.hubs.size(); ++idx) {
                        const auto& hub = node.hubs[idx];
                        if (hub == queryid_thisid.first) {
                            hubIndices.push_back(idx);
                        }
                    }
                    atomToPorts[queryid_thisid.second] = hubIndices;
                }
            }
        }
    }

    // print all maps assoicated with atom index
    for (const auto& [atomIdx, nodeid] : atomToNodeid) {
        printf("atomIdx: %d ports: ", atomIdx);
        for (const auto& port : atomToPorts[atomIdx]) {
            printf("%d ", port);
        }
        printf(" nodeid: %d", atomToNodeid[atomIdx]);
        printf(" element: %s", mol.nodes[atomIdx].ntype.c_str());
        printf(" smarts: %s", atomToSmarts[atomIdx].c_str());
        printf("\n");
    }

    // If any atoms are not matched to a node, throw an error
    std::vector<int> unmatchedAtoms;
    std::string unmatchedAtomsStr = "";
    for (const auto& [nodeid, atom] : mol.nodes) {
        if (!atomToNodeid.count(nodeid)) {
            unmatchedAtoms.push_back(nodeid);
        }
    }
    if (!unmatchedAtoms.empty()) { // TODO: Maybe should be a warning instead of an error
        for (const auto& atomIdx : unmatchedAtoms) {
            unmatchedAtomsStr += std::to_string(atomIdx) + " ";
        }
        std::cout<<groupGraph.printGraph();
        throw std::invalid_argument("Atoms not matched to any node: " + unmatchedAtomsStr);
    }


    printf("Got here\n");
    printf("Group Graph: \n");
    printf("%s\n", groupGraph.printGraph().c_str());

    // Add edges based on connectivity in the original molecule
    for (const auto& [nodeid, dstSet] : mol.edges) {
        int beginAtomIdx = nodeid;
        for (const auto& [dst, order] : dstSet) {
            int endAtomIdx = dst;


            GroupGraph::NodeIDType fromNode = atomToNodeid[beginAtomIdx];
            GroupGraph::NodeIDType toNode = atomToNodeid[endAtomIdx];
            std::string fromSmarts = atomToSmarts[beginAtomIdx];
            std::string toSmarts = atomToSmarts[endAtomIdx];


            // Check if both atoms are part of the same SMARTS subgraph
            if (fromNode != toNode) {
                // Check if the edge already exists
                bool edgeExists = false;
                for (const auto& edge : groupGraph.edges) {
                    if ((std::get<0>(edge) == fromNode && std::get<2>(edge) == toNode) || (std::get<0>(edge) == toNode && std::get<2>(edge) == fromNode)) {
                        edgeExists = true;
                        break;
                    }
                }
                if (!edgeExists) {
                    int src_port = -1;
                    int dst_port = -1;
                    // select first available port for both nodes
                    if (groupGraph.n_free_ports(fromNode) > 0) {
                        src_port = atomToPorts[beginAtomIdx][0];
                    }
                    if (groupGraph.n_free_ports(toNode) > 0) {
                        dst_port = atomToPorts[endAtomIdx][0];
                    }
                    printf("Adding edge from %s to %s\n", fromSmarts.c_str(), toSmarts.c_str());
                    printf("Adding edge from %d to %d\n", fromNode, toNode);
                    printf("    from port: %d\n", src_port);
                    printf("    to port: %d\n", dst_port);

                    if (src_port == -1 || dst_port == -1) {
                        throw std::invalid_argument("No free ports available for edge, occured while adding edge from " + fromSmarts + " to " + toSmarts);
                    }

                    groupGraph.addEdge(
                        std::make_tuple(fromNode, src_port), 
                        std::make_tuple(toNode, dst_port)
                    );
                }
            }
        }
    }
    return groupGraph;
}