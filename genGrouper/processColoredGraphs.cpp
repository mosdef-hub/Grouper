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

#include "dataStructures.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
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
                int dst_index = automorphism[i].second;
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

// Function to generate all non-isomorphic colored graphs
std::vector<GroupGraph> generate_non_isomorphic_colored_graphs(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::unordered_map<int, std::string>& int_to_node_type,
    const std::unordered_map<int, std::string>& int_to_smiles,
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    const std::unordered_map<std::string, std::vector<int>>& node_type_to_hub,
    const std::vector<int> node_colors
) {
    std::unordered_set<GroupGraph> unique_graphs;
    std::unordered_map<std::pair<int, int>, int> possible_edge_colors;

    // Calculate the possible edge colors for each edge
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        int n_src_ports = node_types.at(int_to_node_type.at(node_colors[src])).size();
        int n_dst_ports = node_types.at(int_to_node_type.at(node_colors[dst])).size();
        int n_colors = n_src_ports * n_dst_ports;
        possible_edge_colors[{src, dst}] = n_colors;
    }

    // Prepare to generate all possible colorings
    std::vector<std::vector<int>> all_colorings;
    std::vector<int> current_coloring(possible_edge_colors.size(), 0);
    generate_all_colorings(possible_edge_colors, all_colorings, current_coloring, 0);

    // Create a graph object to compute automorphisms
    GroupGraph gG;
    // Add nodes to the graph
    for (int i = 0; i < node_colors.size(); ++i) {
        gG.addNode(
            int_to_node_type.at(node_colors[i]), 
            int_to_smiles.at(node_colors[i]), 
            node_types.at(int_to_node_type.at(node_colors[i])), 
            node_type_to_hub.at(int_to_node_type.at(node_colors[i]))
        );
    }

    // Compute the edge automorphisms of the graph
    std::vector<std::vector<std::pair<int, int>>> automorphisms = gG.edgeAut(edge_list);
    for (const auto& automorphism : automorphisms) {
        for (const auto& edge : automorphism) {
            std::cout << edge.first << " " << edge.second << " | ";
        }
        std::cout << std::endl;
    }

    // Apply automorphisms to filter out redundant colorings
    std::vector<std::vector<int>> unique_colorings = apply_edge_automorphisms(all_colorings, automorphisms);
    // Iterate over unique colorings
    std::vector<GroupGraph> result;
    for (const auto& coloring : unique_colorings) {
        GroupGraph gG;

        // Add nodes to the graph
        for (int i = 0; i < node_colors.size(); ++i) {
            gG.addNode(
                int_to_node_type.at(node_colors[i]), 
                int_to_smiles.at(node_colors[i]), 
                node_types.at(int_to_node_type.at(node_colors[i])), 
                node_type_to_hub.at(int_to_node_type.at(node_colors[i]))
            );
        }
        // Add edges with the current coloring
        size_t edge_index = 0;
        try {
            for (const auto& edge : edge_list) {
                int src = edge.first;
                int dst = edge.second;
                int color = coloring[edge_index++];
                std::pair<int,int> colorPort = color_to_ports(color, node_types.at(int_to_node_type.at(node_colors[src])), node_types.at(int_to_node_type.at(node_colors[dst])));
                int srcPort = colorPort.first;
                int dstPort = colorPort.second;

                gG.addEdge({src, srcPort}, {dst, dstPort});
            }
        } catch (const std::exception& e) {
            continue;
        }

        // Check if the graph is unique considering permutations
        if (unique_graphs.find(gG) == unique_graphs.end()) {
            unique_graphs.insert(gG);
            result.push_back(gG);
        }
    }

    return result;
}


std::unordered_set<std::string> process_nauty_output(
    const std::string& line, 
    const std::unordered_set<GroupGraph::Node>& node_defs,
    bool verbose
) {

    std::vector<GroupGraph> group_graphs_list;
    std::unordered_set<std::string> graph_basis;

    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        return {};
    }

    std::string node_description_str = line.substr(0, split_pos);
    std::string edge_description_str = line.substr(split_pos + 2);

    std::istringstream node_desc_iss(node_description_str);
    std::vector<std::string> node_description(std::istream_iterator<std::string>{node_desc_iss},
                                              std::istream_iterator<std::string>());

    std::istringstream edge_desc_iss(edge_description_str);
    std::vector<std::string> edge_description(std::istream_iterator<std::string>{edge_desc_iss},
                                              std::istream_iterator<std::string>());

    std::vector<std::pair<int, int>> edge_list;
    for (size_t i = 0; i < edge_description.size(); i += 2) {
        edge_list.emplace_back(std::stoi(edge_description[i]), std::stoi(edge_description[i + 1]));
    }

    

    int n_vertices = std::stoi(node_description[0]);
    // int n_edges = std::stoi(node_description[1]);
    std::vector<int> colors;
    for (size_t i = 2; i < node_description.size(); ++i) {
        colors.push_back(std::stoi(node_description[i]));
    }

    // Error checking
    // for (const auto& color : colors) {
    //     std::cout << "colors: " << color << " ";
    // }
    if (node_defs.size() < static_cast<size_t>(*std::max_element(colors.begin(), colors.end()) + 1)) {
        throw std::runtime_error("Number of nodes in node_defs does not match the number of nodes in the nauty_output_file...");
    }


    // Actual function starts here
    std::vector<std::pair<int, int>> non_colored_edge_list;
    std::unordered_map<int, std::string> int_to_node_type;
    std::unordered_map<std::string, std::vector<int>> node_types;
    std::unordered_map<std::string, std::vector<int>> node_type_to_hub;
    std::unordered_map<int, std::string> int_to_smiles;
    std::unordered_map<std::string, std::string> type_to_smiles;

    for (const auto& node : node_defs) {
        node_types[node.ntype] = node.ports;
        type_to_smiles[node.ntype] = node.smiles;
    }

    for (const auto& node: node_defs) {
        for (const auto& h : node.hubs) {
            node_type_to_hub[node.ntype].push_back(h);
        }
    }

    int i = 0;
    for (const auto& [node_type, ports] : node_types) {
        int_to_node_type[i] = node_type;
        int_to_smiles[i] = type_to_smiles[node_type];
        i++;
    }


    // Add nodes to the graph
    GroupGraph gG;
    for (int i = 0; i < n_vertices; ++i) {
        gG.addNode(
            int_to_node_type.at(colors[i]), 
            int_to_smiles.at(colors[i]), 
            node_types.at(int_to_node_type.at(colors[i])), 
            node_type_to_hub.at(int_to_node_type.at(colors[i]))
        );
    }



    std::vector<GroupGraph> graphs = generate_non_isomorphic_colored_graphs(
        edge_list, 
        int_to_node_type, 
        int_to_smiles, 
        node_types, 
        node_type_to_hub,
        colors
    );

    // Convert the graphs to smiles
    for (const auto& graph : graphs) {
        graph_basis.insert(graph.toSmiles());
    }


    return graph_basis;
}

