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
#include "GroupGraph.cpp"

std::vector<std::vector<std::pair<int, int>>> cartesian_product(const std::vector<std::vector<std::pair<int, int>>>& v) {
    std::vector<std::vector<std::pair<int, int>>> result;
    if (v.empty()) {
        return result;
    }

    // Calculate the total size of the Cartesian product
    auto size = std::accumulate(v.begin(), v.end(), 1,
                                [](size_t s, const std::vector<std::pair<int, int>>& sub) { return s * sub.size(); });

    // Resize the result vector to hold all combinations
    result.resize(size, std::vector<std::pair<int, int>>(v.size()));

    // Initialize counters for each vector of pairs
    std::vector<size_t> counters(v.size(), 0);

    // Populate the result vector with combinations
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < v.size(); ++j) {
            result[i][j] = v[j][counters[j]];
        }

        // Update counters for the next combination
        for (size_t j = v.size(); j-- > 0;) {
            if (++counters[j] < v[j].size()) {
                break;
            }
            counters[j] = 0;
        }
    }

    return result;
}

std::vector<GroupGraph> process_nauty_graph_vcolg_output(const std::string& line, 
                                                         const std::unordered_map<std::string, std::vector<int>>& node_types,
                                                         const std::unordered_map<int, std::string>& int_to_node_type,
                                                         bool verbose = false) {
    std::vector<GroupGraph> group_graphs_list;

    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        return group_graphs_list;
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
    int n_edges = std::stoi(node_description[1]);
    std::vector<int> colors;
    for (size_t i = 2; i < node_description.size(); ++i) {
        colors.push_back(std::stoi(node_description[i]));
    }

    std::vector<std::pair<int, int>> non_colored_edge_list;

    GroupGraph gG(node_types);
    for (int i = 0; i < n_vertices; ++i) {
        gG.addNode(i, int_to_node_type.at(colors[i]));
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
    for (const auto& [index, combinations] : edge_index_to_edge_color) {
        edge_colors.push_back(combinations);
    }

    auto port_combos = cartesian_product(edge_colors);
    for (const auto& port_combo : port_combos) {
        try {
            GroupGraph temp_gG = gG;
            for (size_t j = 0; j < edge_list.size(); ++j) {
                temp_gG.addEdge(
                    edge_list[j].first, port_combo[j].first,
                    edge_list[j].second, port_combo[j].second
                    );
            }
            group_graphs_list.push_back(temp_gG);
        } catch (const std::exception& e) {
            if (verbose) {
                std::cerr << "Couldn't produce graph from vcolg output: " << e.what() << std::endl;
            }
        }
    }
    return group_graphs_list;
}


// int main() {

//     std::unordered_map<std::string, std::vector<int>> node_types = {
//         {"N", {0}},
//         {"CO", {0, 1}},
//         {"CC", {0, 1, 2, 3}}
//     };
//     std::unordered_map<int, std::string> int_to_node_type = {
//         {0, "N"},
//         {1, "CO"},
//         {2, "CC"},
//         {3, "CC"},
//         {4, "CC"}
//     };
//     std::ifstream input_file("vcolg_out.txt");
//     if (!input_file.is_open()) {
//         std::cerr << "Error opening input file." << std::endl;
//         return 1;
//     }
//     std::string line;
//     while (std::getline(input_file, line)) {
//         if (line.empty()) {
//             continue; // Skip empty lines
//         }
//         auto result = process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, true);
//         for (const auto& gG : result) {
//             std::cout << gG.calculateMemoryUsage() << std::endl;
//         }
//         std::cout << "Processed line. Result size: " << result.size() << std::endl;
//     }
//     input_file.close();

//     return 0;
// }
