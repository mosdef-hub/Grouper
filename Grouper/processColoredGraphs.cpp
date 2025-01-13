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
#include <boost/functional/hash.hpp>

#include "dataStructures.hpp"
#include "debugTools.cpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/GraphMol.h>


struct hash_pair {
    std::size_t operator()(const std::pair<int, int>& p) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        return seed;
    }
};


std::pair<int, int> color_to_ports(int color, const std::vector<int>& src_ports, const std::vector<int>& dst_ports) {
    int src_port = src_ports[color / dst_ports.size()];
    int dst_port = dst_ports[color % dst_ports.size()];
    return {src_port, dst_port};
}

std::tuple< int, std::vector<int>, std::vector<std::pair<int, int>> > parse_nauty_graph_line(const std::string& line, const std::unordered_set<GroupGraph::Node>& node_defs) {
    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        throw std::runtime_error("Invalid nauty output line...");
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

    if (node_defs.size() < static_cast<size_t>(*std::max_element(colors.begin(), colors.end()) + 1)) {
        throw std::runtime_error("Number of nodes in node_defs does not match the number of nodes in the nauty_output_file...");
    }

    return {n_vertices, colors, edge_list};
}

bool check_max_bond_not_exceeded(
    const std::vector<std::pair<int,int>>& edge_list,
    const std::vector<int>& colors, 
    const std::unordered_map<std::string, std::vector<int>>& node_types, 
    const std::unordered_map<int, std::string>& int_to_node_type) {


    std::vector<int> node_bond_count;
    for (std::size_t i = 0; i < colors.size(); ++i) {
        std::string node_type = int_to_node_type.at(colors[i]);
        int n_ports = node_types.at(node_type).size();
        node_bond_count.push_back(n_ports);
    }
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        // int src_color = colors[src];
        // int dst_color = colors[dst];
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

// Function to filter orbits based on node-node color pairs
std::vector<std::unordered_set<std::pair<int, int>, hash_pair>> filterOrbits(
    const std::vector<std::pair<int, int>>& edge_list,
    const int* edge_orbits,  // Pointer to array of orbit assignments
    const std::vector<int>& colors) 
{
    std::vector<std::unordered_set<std::pair<int, int>, hash_pair>> edge_orbits_filtered;
    std::unordered_map<std::pair<int, int>, std::unordered_set<std::pair<int, int>, hash_pair>, hash_pair> color_orbit_map;
    int num_edges = edge_list.size();

    // Map each node-node color pair to its associated edges
    for (int i = 0; i < num_edges; i++) {
        int src_color = colors[edge_list[i].first];
        int dst_color = colors[edge_list[i].second];
        std::pair<int, int> color_pair = std::make_pair(src_color, dst_color);
        color_orbit_map[color_pair].insert(edge_list[i]);
    }

    // Create new orbits based on node-node color pairs
    for (auto const& [color_pair, edges] : color_orbit_map) {
        edge_orbits_filtered.push_back(edges);
    }

    return edge_orbits_filtered;
}

std::vector<std::vector<int>> generateOrbitCombinations(
    const std::unordered_set<std::pair<int, int>, hash_pair>& edge_orbit, 
    const std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair>& available_colors) 
{
    std::vector<std::vector<int>> all_combinations;
    std::vector<int> all_colors;

    std::vector<int> edge_colors = available_colors.at(*edge_orbit.begin());
    int n_edges = edge_orbit.size();

    // recursive function to generate all possible combinations
    std::function<void(std::vector<int>, int, std::vector<int>, int)> generate_combinations = [&](std::vector<int> elements, std::size_t n_selections, std::vector<int> current_combination = {} , std::size_t i = 0) {
        if (current_combination.size() == n_selections) {
            all_combinations.push_back(current_combination);
            return;
        }
        for (std::size_t j = i; j < elements.size(); ++j) {
            current_combination.push_back(elements[j]);
            generate_combinations(elements, n_selections, current_combination, j + 1);
            current_combination.pop_back();
        }
    };

    generate_combinations(edge_colors, n_edges, {}, 0);

    return all_combinations;
}

// A helper function to compute all non-isomorphic colorings
std::vector<std::vector<int>> generateNonAutomorphicEdgeColorings(
    std::vector<std::pair<int, int>>& edge_list, // Edges in the graph
    std::vector<std::unordered_set<std::pair<int, int>, hash_pair>>& edge_orbits, // Edge orbits
    std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair>& available_colors // Available colors for each edge
) {
    std::vector<int> current_coloring(edge_list.size(), -1);  // Initialize coloring (all edges uncolored)
    std::vector<std::vector<int>> all_colorings;  // Store all valid colorings
    std::vector<std::vector<std::vector<int>>> colored_orbits;  // Store all valid colorings for each orbit

    // Generate possible colorings for each orbit
    for (const auto& orbit : edge_orbits) {
        std::vector<std::vector<int>> orbit_colorings = generateOrbitCombinations(orbit, available_colors);
        colored_orbits.push_back(orbit_colorings);
    }

    // Helper function to recursively combine colorings from orbits
    std::function<void(int)> combineOrbits = [&](std::size_t i) {
        if (i == colored_orbits.size()) {
            // If we have processed all orbits, store the current coloring
            all_colorings.push_back(current_coloring);
            return;
        }

        // Loop through all colorings of the current orbit
        for (const auto& coloring : colored_orbits[i]) {
            // Apply the coloring to the corresponding edges in the current orbit
            int color_idx = 0;
            for (const auto& edge : edge_orbits[i]) {
                // Find the index of the edge in the edge_list
                auto it = std::find(edge_list.begin(), edge_list.end(), edge);
                if (it != edge_list.end()) {
                    int edge_index = std::distance(edge_list.begin(), it);
                    current_coloring[edge_index] = coloring[color_idx];
                }
                ++color_idx;
            }

            // Recursively process the next orbit
            combineOrbits(i + 1);

            // Undo the coloring for backtracking
            color_idx = 0;
            for (const auto& edge : edge_orbits[i]) {
                auto it = std::find(edge_list.begin(), edge_list.end(), edge);
                if (it != edge_list.end()) {
                    int edge_index = std::distance(edge_list.begin(), it);
                    current_coloring[edge_index] = -1;
                }
                ++color_idx;
            }
        }
    };

    // Start combining from the first orbit
    combineOrbits(0);

    return all_colorings;
}



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
    bool verbose,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats // Pass nauty structures
) {
    // initializeLogger();
    // logger.log("This is a debug message", Logger::LogLevel::DEBUG);

    // Process the nauty output line
    auto [n_vertices, colors, edge_list] = parse_nauty_graph_line(line, node_defs);

    // Actual function starts here
    std::vector<std::pair<int, int>> non_colored_edge_list;
    std::unordered_map<int, std::string> int_to_node_type;
    std::unordered_map<std::string, std::vector<int>> node_types;
    std::unordered_map<std::string, std::vector<int>> node_type_to_hub;
    std::unordered_map<int, std::string> int_to_smiles;
    std::unordered_map<std::string, std::string> type_to_smiles;
    std::vector<GroupGraph> group_graphs_list;
    std::unordered_set<std::string> canon_set;


    // Create necessary maps
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
    
    // Create a histogram of node types for positive constraints
    if (positiveConstraints.size() > 0) {
        std::unordered_map<std::string, int> node_hist;
        for (const auto& node : node_defs) {
            node_hist[node.ntype] = 0;
        }
        for (const auto& c : colors) {
            node_hist[int_to_node_type.at(c)] += 1;
        }
        
        // Check if the number of nodes of each type is less than the number listed in the positive constraints
        for (const auto& [node_type, count] : node_hist) {
            if (positiveConstraints.find(node_type) == positiveConstraints.end()) {
                continue;
            }
            if (count < positiveConstraints.at(node_type)) {
                return;
            }
        }
    }

    // Check if any node has too many bonds
    if (!check_max_bond_not_exceeded(edge_list, colors, node_types, int_to_node_type)) {
        return;
    }

    // Add nodes to the graph to create node colored graph
    GroupGraph gG;
    for (int i = 0; i < n_vertices; ++i) {
        gG.addNode(
            int_to_node_type.at(colors[i]), 
            int_to_smiles.at(colors[i]), 
            node_type_to_hub.at(int_to_node_type.at(colors[i]))
        );
    }

    // Generate all possible colorings for the edges 
    // TODO: this isn't correct since edge colors need to be different for each node-node pair (red(2),blue(2))=(0,1,2,3) is not the same as (blue(2), blue(2)) = (4,5,6,7)
    std::vector<std::vector<int>> all_colorings;
    std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> possible_edge_colors;
    std::vector<int> current_coloring(possible_edge_colors.size(), 0);
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        int n_src_ports = node_types.at(int_to_node_type.at(colors[src])).size();
        int n_dst_ports = node_types.at(int_to_node_type.at(colors[dst])).size();
        int n_colors = n_src_ports * n_dst_ports;
        std::vector<int> e_colors(n_colors);
        std::iota(e_colors.begin(), e_colors.end(), 0);
        possible_edge_colors[{src, dst}] = e_colors;
    }


    // Compute edge orbits using nauty
    int* edge_orbits = gG.computeEdgeOrbits(edge_list, g, lab, ptn, orbits, options, stats);

    // Filter orbits based on node-node color pairs
    std::vector<std::unordered_set<std::pair<int, int>, hash_pair>> edge_orbits_filtered = filterOrbits(edge_list, edge_orbits, colors);

    // Compute all non-automorphic colorings
    std::vector<std::vector<int>> unique_colorings = generateNonAutomorphicEdgeColorings(edge_list, edge_orbits_filtered, possible_edge_colors);

    // Iterate over unique colorings
    for (const auto& coloring: unique_colorings){
        gG.clearEdges();
        // Add edges with the current coloring
        size_t edge_index = 0;
        try {
            for (const auto& edge : edge_list) {
                int src = edge.first;
                int dst = edge.second;
                int color = coloring[edge_index++];
                std::pair<int,int> colorPort = color_to_ports(color, node_types.at(int_to_node_type.at(colors[src])), node_types.at(int_to_node_type.at(colors[dst])));
                int srcPort = colorPort.first;
                int dstPort = colorPort.second;

                gG.addEdge({src, srcPort}, {dst, dstPort});
            }
        } catch (const std::exception& e) {
            // std::cout << "Error adding edge: " << e.what() << std::endl;
            continue;
        }

    //  Check if the graph is unique considering permutations
        if (canon_set.find(gG.toSmiles()) == canon_set.end()) {
            canon_set.insert(gG.toSmiles());
            graph_basis->insert(gG);
        }
    }
}

