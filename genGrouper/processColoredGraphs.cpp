#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>
#include <numeric>
#include <exception>
#include <iterator>
#include <set>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

#include "dataStructures.hpp"
#include "colorPermutations.cuh"
#include "debugTools.hpp"
#include <iterator>
#include <mutex>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/GraphMol.h>


void generate_all_colorings(
    const std::unordered_map<std::pair<int, int>, int>& possible_edge_colors, 
    std::vector<std::vector<int>>& all_colorings, 
    std::vector<int>& current_coloring, 
    size_t edge_index) 
{
    if (edge_index == possible_edge_colors.size()) {
        all_colorings.push_back(current_coloring);
        return;
    }

    auto it = possible_edge_colors.begin();
    std::advance(it, edge_index);
    int num_colors = it->second;

    for (int color = 0; color < num_colors; ++color) {
        current_coloring[edge_index] = color;
        generate_all_colorings(possible_edge_colors, all_colorings, current_coloring, edge_index + 1);
    }
}

// Helper function to apply automorphisms to colorings
std::vector<std::vector<int>> apply_edge_automorphisms(
    const std::vector<std::vector<int>>& colorings, 
    const std::vector<std::vector<std::pair<int, int>>>& automorphisms) 
{
    std::unordered_set<std::vector<int>> unique_colorings;

    for (const auto& coloring : colorings) {
        for (const auto& automorphism : automorphisms) {
            std::vector<int> permuted_coloring(coloring.size());
            for (size_t i = 0; i < automorphism.size(); ++i) {
                // Apply the edge automorphism to the coloring
                int src_index = automorphism[i].first;
                // int dst_index = automorphism[i].second;
                permuted_coloring[i] = coloring[src_index];
            }
            unique_colorings.insert(permuted_coloring);
        }
    }

    return std::vector<std::vector<int>>(unique_colorings.begin(), unique_colorings.end());
}


std::pair<int, int> color_to_ports(int color, const std::vector<int>& src_ports, const std::vector<int>& dst_ports) {
    int src_port = src_ports[color / dst_ports.size()];
    int dst_port = dst_ports[color % dst_ports.size()];
    return {src_port, dst_port};
}

bool check_max_bond_not_exceeded(
    const std::vector<std::pair<int,int>>& edge_list,
    const std::vector<int>& colors, 
    const std::unordered_map<std::string, std::vector<int>>& node_types, 
    const std::unordered_map<int, std::string>& int_to_node_type) {

    std::vector<int> node_bond_count;
    for (int i = 0; i < colors.size(); ++i) {
        std::string node_type = int_to_node_type.at(colors[i]);
        int n_ports = node_types.at(node_type).size();
        node_bond_count.push_back(n_ports);
    }
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        int src_color = colors[src];
        int dst_color = colors[dst];
        node_bond_count[src] -= 1;
        node_bond_count[dst] -= 1;
    }
    for (const auto& count : node_bond_count) {
        if (count < 0) {
            return false;
        }
    }
    return true;
}


// Helper function to find the representative of an orbit (union-find technique)
int findRepresentative(int edgeIndex, std::vector<int>& edgeOrbit) {
    if (edgeOrbit[edgeIndex] != edgeIndex) {
        edgeOrbit[edgeIndex] = findRepresentative(edgeOrbit[edgeIndex], edgeOrbit);
    }
    return edgeOrbit[edgeIndex];
};

// Helper function to unite two edge orbits (union-find technique)
void unionOrbits(int edgeA, int edgeB, std::vector<int>& edgeOrbit) {
    int repA = findRepresentative(edgeA, edgeOrbit);
    int repB = findRepresentative(edgeB, edgeOrbit);
    if (repA != repB) {
        edgeOrbit[repA] = repB; // Unite the two orbits
    }
};

std::vector<std::vector<int>> computeEdgeOrbits(
    const std::vector<std::pair<int, int>>& edges,  // List of edges (pairs of vertices)
    const std::vector<std::vector<int>>& automorphisms // List of automorphisms (vertex permutations)
) {
    int numEdges = edges.size();
    int numAutomorphisms = automorphisms.size();
    
    // Initialize edge orbits: Initially, each edge is in its own orbit
    std::vector<int> edgeOrbit(numEdges);
    for (int i = 0; i < numEdges; ++i) {
        edgeOrbit[i] = i; // Each edge is its own orbit representative
    }

    // Apply each automorphism to all edges
    for (const auto& automorphism : automorphisms) {
        for (int i = 0; i < numEdges; ++i) {
            int v1 = edges[i].first;
            int v2 = edges[i].second;

            // Apply automorphism to both vertices of the edge
            int permutedV1 = automorphism[v1];
            int permutedV2 = automorphism[v2];

            // Find the edge (permutedV1, permutedV2) in the edge list
            for (int j = 0; j < numEdges; ++j) {
                if ((edges[j].first == permutedV1 && edges[j].second == permutedV2) ||
                    (edges[j].first == permutedV2 && edges[j].second == permutedV1)) {
                    // Unite the two edges' orbits
                    unionOrbits(i, j, edgeOrbit);
                    break;
                }
            }
        }
    }

    // Group edges by their orbit representative
    std::unordered_map<int, std::vector<int>> orbitGroups;
    for (int i = 0; i < numEdges; ++i) {
        int rep = findRepresentative(i, edgeOrbit);
        orbitGroups[rep].push_back(i);
    }

    // Collect the results into a vector of edge orbits
    std::vector<std::vector<int>> edgeOrbits;
    for (const auto& orbitGroup : orbitGroups) {
        edgeOrbits.push_back(orbitGroup.second);
    }

    return edgeOrbits;
}

// // Function to generate all non-isomorphic colored graphs
// std::vector<GroupGraph> generate_non_isomorphic_colored_graphs(
//     const std::vector<std::pair<int, int>>& edge_list,
//     const std::unordered_map<int, std::string>& int_to_node_type,
//     const std::unordered_map<int, std::string>& int_to_smiles,
//     const std::unordered_map<std::string, std::vector<int>>& node_types,
//     const std::unordered_map<std::string, std::vector<int>>& node_type_to_hub,
//     const std::vector<int> node_colors) {
//     std::unordered_set<GroupGraph> unique_graphs;
//     std::unordered_map<std::pair<int, int>, int> possible_edge_colors;

//     // Calculate the possible edge colors for each edge
//     for (const auto& edge : edge_list) {
//         int src = edge.first;
//         int dst = edge.second;
//         int n_src_ports = node_types.at(int_to_node_type.at(node_colors[src])).size();
//         int n_dst_ports = node_types.at(int_to_node_type.at(node_colors[dst])).size();
//         int n_colors = n_src_ports * n_dst_ports;
//         possible_edge_colors[{src, dst}] = n_colors;
//     }

//     // Prepare to generate all possible colorings
//     std::vector<std::vector<int>> all_colorings;
//     std::vector<int> current_coloring(possible_edge_colors.size(), 0);
//     generate_all_colorings(possible_edge_colors, all_colorings, current_coloring, 0);

//     // Create a graph object to compute automorphisms
//     GroupGraph gG;
//     // Add nodes to the graph
//     for (std::vector<int>::size_type i = 0; i < node_colors.size(); ++i) {
//         gG.addNode(
//             int_to_node_type.at(node_colors[i]), 
//             int_to_smiles.at(node_colors[i]), 
//             node_types.at(int_to_node_type.at(node_colors[i])), 
//             node_type_to_hub.at(int_to_node_type.at(node_colors[i]))
//         );
//     }

//     // Compute the edge automorphisms of the graph
//     std::vector<std::vector<std::pair<int, int>>> automorphisms = gG.edgeAut(edge_list);

//     // Apply automorphisms to filter out redundant colorings
//     std::vector<std::vector<int>> unique_colorings = apply_edge_automorphisms(all_colorings, automorphisms);
//     // Iterate over unique colorings
//     std::vector<GroupGraph> result;
//     for (const auto& coloring : unique_colorings) {
//         GroupGraph gG;

//         // Add nodes to the graph
//         for (std::vector<int>::size_type i = 0; i < node_colors.size(); ++i) {
//             gG.addNode(
//                 int_to_node_type.at(node_colors[i]), 
//                 int_to_smiles.at(node_colors[i]), 
//                 node_types.at(int_to_node_type.at(node_colors[i])), 
//                 node_type_to_hub.at(int_to_node_type.at(node_colors[i]))
//             );
//         }
//         // Add edges with the current coloring
//         size_t edge_index = 0;
//         try {
//             for (const auto& edge : edge_list) {
//                 int src = edge.first;
//                 int dst = edge.second;
//                 int color = coloring[edge_index++];
//                 std::pair<int,int> colorPort = color_to_ports(color, node_types.at(int_to_node_type.at(node_colors[src])), node_types.at(int_to_node_type.at(node_colors[dst])));
//                 int srcPort = colorPort.first;
//                 int dstPort = colorPort.second;

//                 gG.addEdge({src, srcPort}, {dst, dstPort});
//             }
//         } catch (const std::exception& e) {
//             continue;
//         }

//         // Check if the graph is unique considering permutations
//         if (unique_graphs.find(gG) == unique_graphs.end()) {
//             unique_graphs.insert(gG);
//             result.push_back(gG);
//         }
//     }

//     return result;
// }


// Assume GroupGraph and related types are defined elsewhere

// Initialize logger
// Logger& logger = Logger::getInstance();

// void initializeLogger() {
//     logger.setLogLevel(Logger::LogLevel::DEBUG);
//     logger.enableFileLogging("log.txt");
// }

void process_nauty_output(
    const std::string& line, 
    const std::unordered_set<GroupGraph::Node>& node_defs,
    std::unordered_set<GroupGraph>* graph_basis,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    bool verbose
) {
    
    // initializeLogger();

    std::vector<GroupGraph> group_graphs_list;
    std::unordered_set<std::string> canon_set;

    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        return;
    }

    std::string node_description_str = line.substr(0, split_pos);
    std::string edge_description_str = line.substr(split_pos + 2);


    std::istringstream node_desc_iss(node_description_str);
    std::vector<std::string> node_description(std::istream_iterator<std::string>{node_desc_iss},
                                              std::istream_iterator<std::string>());

    std::istringstream edge_desc_iss(edge_description_str);
    std::vector<std::string> edge_description(std::istream_iterator<std::string>{edge_desc_iss},
                                              std::istream_iterator<std::string>());

    // Convert edge descriptions to integer pairs
    std::vector<std::pair<int, int>> edge_list;
    try {
        for (size_t i = 0; i < edge_description.size(); i += 2) {
            edge_list.emplace_back(std::stoi(edge_description[i]), std::stoi(edge_description[i + 1]));
        }
    } catch (const std::exception& e) {
        return;
    }

    int n_vertices;
    try {
        n_vertices = std::stoi(node_description[0]);
    } catch (const std::exception& e) {
        return;
    }

    std::vector<int> colors;
    try {
        for (size_t i = 2; i < node_description.size(); ++i) {
            colors.push_back(std::stoi(node_description[i]));
        }
    } catch (const std::exception& e) {
        return;
    }

    try {
        if (node_defs.size() < static_cast<size_t>(*std::max_element(colors.begin(), colors.end()) + 1)) {
            throw std::runtime_error("Number of nodes in node_defs does not match the number of nodes in the nauty_output_file...");
        }
    } catch (const std::exception& e) {
        throw;  // Re-throwing the exception after logging
    }


    // Initialize data structures
    std::vector<std::pair<int, int>> non_colored_edge_list;
    std::unordered_map<int, std::string> int_to_node_type;
    std::unordered_map<std::string, std::vector<int>> node_types;
    std::unordered_map<std::string, std::vector<int>> node_type_to_hub;
    std::unordered_map<int, std::string> int_to_smiles;
    std::unordered_map<std::string, std::string> type_to_smiles;

    // Populate node types and smiles
    for (const auto& node : node_defs) {
        node_types[node.ntype] = node.ports;
        type_to_smiles[node.ntype] = node.smiles;
    }

    for (const auto& node : node_defs) {
        for (const auto& h : node.hubs) {
            node_type_to_hub[node.ntype].push_back(h);
        }
    }

    // Map integers to node types and smiles
    int i = 0;
    for (const auto& [node_type, ports] : node_types) {
        int_to_node_type[i] = node_type;
        int_to_smiles[i] = type_to_smiles[node_type];
        i++;
    }

    // Create a histogram of node types for positive constraints
    std::unordered_map<std::string, int> node_hist;
    for (const auto& node : node_defs) {
        node_hist[node.ntype] = 0;
    }
    for (const auto& c : colors) {
        node_hist.at(int_to_node_type.at(c)) += 1;
    }

    // Check if any node has too many bonds
    if (!check_max_bond_not_exceeded(edge_list, colors, node_types, int_to_node_type)) {
        return;
    }


    // Check positive constraints
    for (const auto& [node_type, count] : node_hist) {
        if (positiveConstraints.find(node_type) == positiveConstraints.end()) {
            continue;
        }
        if (count < positiveConstraints.at(node_type)) {
            return;
        }
    }

    // Add nodes to the graph
    GroupGraph gG;
    try {
        for (int i = 0; i < n_vertices; ++i) {
            gG.addNode(
                int_to_node_type.at(colors[i]), 
                int_to_smiles.at(colors[i]), 
                node_types.at(int_to_node_type.at(colors[i])), 
                node_type_to_hub.at(int_to_node_type.at(colors[i]))
            );
        }
    } catch (const std::exception& e) {
        return;
    }

    // Prepare to generate all possible colorings
    std::vector<std::vector<int>> all_colorings;
    std::unordered_map<std::pair<int, int>, int> possible_edge_colors;

    // Calculate the possible edge colors for each edge
    try {
        for (const auto& edge : edge_list) {
            int src = edge.first;
            int dst = edge.second;
            int n_src_ports = node_types.at(int_to_node_type.at(colors[src])).size();
            int n_dst_ports = node_types.at(int_to_node_type.at(colors[dst])).size();
            int n_colors = n_src_ports * n_dst_ports;
            possible_edge_colors[{src, dst}] = n_colors;
        }
    } catch (const std::exception& e) {
        return;
    }


    // Placeholder for generate_all_colorings function
    // You need to ensure this function is defined and handles the current_coloring size appropriately
    std::vector<int> current_coloring(possible_edge_colors.size(), 0);
    generate_all_colorings(possible_edge_colors, all_colorings, current_coloring, 0);

    // Compute the edge automorphisms of the graph
    std::vector<std::vector<std::pair<int, int>>> edge_automorphisms;
    try {
        edge_automorphisms = gG.edgeAut(edge_list);
    } catch (const std::exception& e) {
        return;
    }

    // Placeholder for apply_edge_automorphisms_gpu function
    // Ensure this function is defined and correctly integrates with your GPU processing
    std::vector<std::vector<int>> unique_colorings = apply_edge_automorphisms_gpu(all_colorings, edge_automorphisms);

    // Iterate over unique colorings
    for (const auto& coloring : unique_colorings) {
        // Remove edges from the graph
        gG.clearEdges();
        // Add edges with the current coloring
        size_t edge_index = 0;
        try {
            for (const auto& edge : edge_list) {
                int src = edge.first;
                int dst = edge.second;
                int color = coloring[edge_index++];
                std::pair<int, int> colorPort = color_to_ports(color, node_types.at(int_to_node_type.at(colors[src])), node_types.at(int_to_node_type.at(colors[dst])));
                int srcPort = colorPort.first;
                int dstPort = colorPort.second;

                gG.addEdge({src, srcPort}, {dst, dstPort});
            }
        } catch (const std::exception& e) {
        }

        // Check if the graph is unique considering permutations
        try {
            std::string smiles = gG.toSmiles();
            if (canon_set.find(smiles) == canon_set.end()) {
                canon_set.insert(smiles);
                graph_basis->insert(gG);
            } else {
            }
        } catch (const std::exception& e) {
            continue;
        }
    }

    // logger.log("Exiting process_nauty_output function.", Logger::INFO);
}

    
    // std::vector<GroupGraph> graphs = generate_non_isomorphic_colored_graphs(
    //     edge_list, 
    //     int_to_node_type, 
    //     int_to_smiles, 
    //     node_types, 
    //     node_type_to_hub,
    //     colors
    // );

    // for (const auto& graph : graphs) {
    //     if (canon_set.find(graph.toSmiles()) == canon_set.end()) {
    //         canon_set.insert(graph.toSmiles());
    //         graph_basis->insert(graph);
    //     }
    // }
    


