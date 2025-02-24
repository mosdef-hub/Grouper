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

std::tuple< int, std::vector<int>, std::vector<std::pair<int, int>> > parse_nauty_graph_line(const std::string& line, const std::unordered_set<GroupGraph::Group>& node_defs) {
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
    std::vector<std::unordered_set<std::pair<int, int>, hash_pair>> edge_orbits_filtered; // Filtered edge orbits with 
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

std::vector<std::vector<std::pair<int, int>>> generateOrbitCombinations(
    const std::unordered_set<std::pair<int, int>, hash_pair>& edge_orbit,
    const std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, hash_pair>& available_colors)
{
    std::vector<std::vector<std::pair<int, int>>> all_combinations;

    // Get possible colors for one edge (assumes same options for entire orbit)
    std::vector<std::pair<int, int>> edge_colors = available_colors.at(*edge_orbit.begin());
    int n_edges = edge_orbit.size();

    // printf("Generating combinations for orbit of size %lu\n", n_edges);
    // printf("Available colors: ");
    // for(const auto& color : edge_colors) {
    //     printf("(%d, %d) ", color.first, color.second);
    // }
    // printf("\n");
    // printf("edge_orbit size: %lu\n", edge_orbit.size());

    // Recursive function to generate all possible combinations
    std::function<void(std::vector<std::pair<int, int>>, int, std::vector<std::pair<int, int>>, int)> generate_combinations =
        [&](std::vector<std::pair<int, int>> elements, std::size_t n_selections, std::vector<std::pair<int, int>> current_combination = {}, std::size_t i = 0)
    {
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

    // printf("Generated %lu combinations for orbit of size %lu\n", all_combinations.size(), n_edges);
    // for(const auto& combination : all_combinations) {
    //     for(const auto& color : combination) {
    //         printf("(%d, %d) ", color.first, color.second);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    return all_combinations;
}

std::vector<std::vector<std::pair<int, int>>> generateNonAutomorphicEdgeColorings(
    const std::vector<std::pair<int, int>>& edge_list, // Edges in the graph
    const std::vector<std::unordered_set<std::pair<int, int>, hash_pair>>& edge_orbits, // Edge orbits
    const std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, hash_pair>& available_colors // Available port pairings for each edge
) {
    std::pair<int, int> default_pair = {-1, -1};
    std::vector<std::vector<std::pair<int, int>>> all_colorings; // Store all valid colorings
    std::vector<std::vector<std::vector<std::pair<int, int>>>> colored_orbits; // Store all valid colorings for each orbit
    std::vector<std::pair<int, int>> current_coloring(edge_list.size(), default_pair); // Initialize coloring (all edges uncolored)

    // Precompute a mapping from edge to its index in edge_list
    std::unordered_map<std::pair<int, int>, int, hash_pair> edge_index_map;
    for (size_t i = 0; i < edge_list.size(); ++i) {
        edge_index_map[edge_list[i]] = i;
    }

    // Generate possible colorings for each orbit
    for (const auto& orbit : edge_orbits) {
        std::vector<std::vector<std::pair<int, int>>> orbit_colorings = generateOrbitCombinations(orbit, available_colors);
        colored_orbits.push_back(orbit_colorings);
    }

    // Convert edge orbits to a vector for easier indexing
    std::vector<std::vector<std::pair<int, int>>> edge_orbits_vec;
    for (const auto& orbit : edge_orbits) {
        edge_orbits_vec.push_back(std::vector<std::pair<int, int>>(orbit.begin(), orbit.end()));
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
            for (size_t color_idx = 0; color_idx < edge_orbits_vec[i].size(); ++color_idx) {
                const auto& edge = edge_orbits_vec[i][color_idx];
                if (edge_index_map.find(edge) != edge_index_map.end()) {
                    int edge_index = edge_index_map[edge];
                    current_coloring[edge_index] = coloring[color_idx];
                }
            }

            // Recursively process the next orbit
            combineOrbits(i + 1);

            // Undo the coloring for backtracking
            for (size_t color_idx = 0; color_idx < edge_orbits_vec[i].size(); ++color_idx) {
                const auto& edge = edge_orbits_vec[i][color_idx];
                if (edge_index_map.find(edge) != edge_index_map.end()) {
                    int edge_index = edge_index_map[edge];
                    current_coloring[edge_index] = default_pair;
                }
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
    const std::unordered_set<GroupGraph::Group>& node_defs,
    std::unordered_set<GroupGraph>* graph_basis,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    const std::unordered_map<std::string, std::vector<int>>& type_to_hub_orbits,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats // Pass nauty structures
) {
    // initializeLogger();
    // logger.log("This is a debug message", Logger::LogLevel::DEBUG);

    // Process the nauty output line
    auto [n_vertices, colors, edge_list] = parse_nauty_graph_line(line, node_defs);

    // Actual function starts here
    std::unordered_map<int, std::string> int_to_node_type;
    std::unordered_map<std::string, std::vector<int>> node_types;
    std::unordered_map<std::string, std::vector<int>> node_type_to_hub;
    std::unordered_map<int, std::string> int_to_pattern;
    std::unordered_map<std::string, std::string> type_to_pattern;
    std::vector<GroupGraph> group_graphs_list;
    std::unordered_set<std::string> canon_set;


    // Create necessary maps
    for (const auto& node : node_defs) {
        node_types[node.ntype] = node.ports;
        type_to_pattern[node.ntype] = node.pattern;
    }
    for (const auto& node: node_defs) {
        for (const auto& h : node.hubs) {
            node_type_to_hub[node.ntype].push_back(h);
        }
    }
    int i = 0;
    for (const auto& [node_type, ports] : node_types) {
        int_to_node_type[i] = node_type;
        int_to_pattern[i] = type_to_pattern[node_type];
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
            int_to_pattern.at(colors[i]),
            node_type_to_hub.at(int_to_node_type.at(colors[i]))
        );
    }

    // Calculate node orbits using nauty then compute edge orbits
    auto [nodeOrbits, edgeOrbits] = gG.computeOrbits(
        edge_list, colors,
        g, lab, ptn, orbits, options, stats // Pass nauty structures
    );

    // Generate all possible colorings for the edges
    std::vector<std::vector<int>> all_colorings;
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int,int>>, hash_pair> possible_edge_colors; // Possible edge colors for each edge (srcNode, dstNode) -> (srcPort1, dstPort1)...(srcPortN, dstPortN)
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        std::vector<int> srcPorts = node_types.at(int_to_node_type.at(colors[src]));
        std::vector<int> dstPorts = node_types.at(int_to_node_type.at(colors[dst]));
        std::vector<std::pair<int,int>> e_colors;
        for (int i = 0; i < srcPorts.size(); i++) { // For now just use all possible colors TODO: change this to use hub orbits to determine reduced set of colors
            for (int j = 0; j < dstPorts.size(); j++) {
                e_colors.push_back({srcPorts[i], dstPorts[j]});
            }
        }
        possible_edge_colors[{src, dst}] = e_colors;
    }

    // Change the edge orbits to a vector of sets for generating colorings
    std::vector<std::unordered_set<std::pair<int, int>, hash_pair>> edge_orbits_vector;
    std::unordered_set<int> unique_orbits;
    for (const auto& orbit : edgeOrbits) {
        unique_orbits.insert(orbit);
    }
    // Put edge orbits in same order into a set
    for (int i = 0; i < unique_orbits.size(); i++) {
        edge_orbits_vector.push_back({});
    }
    std::unordered_map<int,int> orbit_to_index;
    int index = 0;
    for (const auto& orbit : unique_orbits) {
        if (orbit_to_index.find(orbit) == orbit_to_index.end()) {
            orbit_to_index[orbit] = index;
            index++;
        }
    }
    std::unordered_map<std::pair<int, int>, int, hash_pair> edge_to_index;
    for (int i = 0; i < edge_list.size(); i++) {
        edge_to_index[edge_list[i]] = i;
    }
    for (const auto& edge : edge_list) {
        int edge_index = edge_to_index[edge];  // Get unique edge index
        edge_orbits_vector[orbit_to_index[edgeOrbits[edge_index]]].insert(edge);
    }

    // Compute all non-automorphic colorings
    std::vector<std::vector<std::pair<int, int>>> unique_colorings = generateNonAutomorphicEdgeColorings(edge_list, edge_orbits_vector, possible_edge_colors);

    // Iterate over unique colorings
    for (const auto& coloring: unique_colorings){
        gG.clearEdges();
        int edge_index = 0;
        try {
            for (const auto& edge : edge_list) {
                int src = edge.first;
                int dst = edge.second;
                auto [srcPort, dstPort] = coloring[edge_index];
                gG.addEdge({src, srcPort}, {dst, dstPort});
                edge_index++;
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

