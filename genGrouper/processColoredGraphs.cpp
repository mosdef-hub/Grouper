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

#include "dataStructures.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>


void get_sensible_port_combos(
    const std::vector<std::vector<std::pair<int, int>>>& v, 
    const GroupGraph& gG, 
    const std::vector<std::pair<int, int>>& edge_list, 
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    std::unordered_set<std::string>& smiles_set, 
    const bool verbose) {

        
    if (v.empty()) {
        return;
    }

    GroupGraph temp_gG;
    std::string smiles;

    // Calculate the total size of the Cartesian product
    auto size = std::accumulate(v.begin(), v.end(), 1,
                                [](size_t s, const std::vector<std::pair<int, int>>& sub) { return s * sub.size(); });
    if (verbose) {
        std::cout << "\nSize of Cartesian product: " << size << std::endl;
    }

    // Initialize counters for each vector of pairs
    std::vector<size_t> counters(v.size(), 0);

    for (auto i = 0; i < size; ++i) {
        std::vector<std::pair<int, int>> combination;

        for (size_t j = 0; j < v.size(); ++j) {
            combination.push_back(v[j][counters[j]]);
        }

        try {
            temp_gG = gG;
            for (size_t j = 0; j < edge_list.size(); ++j) {
                temp_gG.addEdge(
                    std::make_tuple(edge_list[j].first, combination[j].first),
                    std::make_tuple(edge_list[j].second, combination[j].second)
                    );
            }

            // convert to smiles
            smiles = temp_gG.toSmiles();
            // Convert to hash
            smiles_set.insert(smiles);
        } 
        catch (const std::exception& e) {
            if (verbose) {
                std::cerr << "Couldn't produce graph from vcolg output: " << e.what() << std::endl;
            }
        }

        // Update counters for the next combination
        for (size_t j = v.size(); j-- > 0;) {
            if (++counters[j] < v[j].size()) {
                break;
            }
            counters[j] = 0;
        }
    }
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
    //     std::cout << color << " ";
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



    GroupGraph gG;
    for (int i = 0; i < n_vertices; ++i) {
        gG.addNode(int_to_node_type.at(colors[i]), int_to_smiles.at(colors[i]), node_types.at(int_to_node_type.at(colors[i])), node_type_to_hub.at(int_to_node_type.at(colors[i])));
    }

    std::map<int, std::vector<std::pair<int, int>>> edge_index_to_edge_color;
    for (size_t i = 0; i < edge_list.size(); ++i) {
        int src = edge_list[i].first;
        int dst = edge_list[i].second;
        std::vector<int> src_ports = node_types.at(int_to_node_type.at(colors[src]));
        std::vector<int> dst_ports = node_types.at(int_to_node_type.at(colors[dst]));
        non_colored_edge_list.push_back({src, dst});
        std::vector<std::pair<int, int>> port_combinations;
        for (int sp : src_ports) {
            for (int dp : dst_ports) {
                port_combinations.push_back({sp, dp});
            }
        }
        edge_index_to_edge_color[i] = port_combinations;
    }

    // Use a graph library to create a graph from non_colored_edge_list and check degrees
    std::map<int, int> node_degrees; // Placeholder for node degrees
    for (const auto& edge : non_colored_edge_list) {
        node_degrees[edge.first]++;
        node_degrees[edge.second]++;
    }
    for (const auto& [node, degree] : node_degrees) {
        if (degree > node_types.at(int_to_node_type.at(colors[node])).size()) {
            return {};
        }
    }

    std::vector<std::vector<std::pair<int, int>>> edge_colors;
    for (const auto& [index, port_combinations] : edge_index_to_edge_color) {
        edge_colors.push_back(port_combinations);
    }


    // Generate and filter port combinations
    get_sensible_port_combos(
        edge_colors,
        gG,
        non_colored_edge_list, 
        node_types, 
        graph_basis,
        verbose
    );


    return graph_basis;
}

