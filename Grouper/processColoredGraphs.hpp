
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <map>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <sstream>
#include <iterator>
#include <tuple>

#include <libpq-fe.h>

#ifdef __cplusplus
    #define _Thread_local thread_local
#else
    #define _Thread_local _Thread_local  // Use C's definition
#endif



#include "nauty/nauty.h"

#include "dataStructures.hpp"

class GroupGraph; // Forward declaration of GroupGraph

void generate_all_colorings(
    const std::unordered_map<std::pair<int, int>, int>& possible_edge_colors, 
    std::vector<std::vector<int>>& all_colorings, 
    std::vector<int>& current_coloring, 
    size_t edge_index
);

std::vector<std::vector<int>> apply_automorphisms(
    const std::vector<std::vector<int>>& colorings, 
    const std::vector<std::vector<int>>& automorphisms
);

std::vector<GroupGraph> generate_non_isomorphic_colored_graphs(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::unordered_map<int, std::string>& int_to_node_type,
    const std::unordered_map<std::string, std::string>& int_to_smarts,
    const std::unordered_map<std::string, std::vector<std::string>>& node_types,
    const std::unordered_map<std::string, std::string>& node_type_to_hub
);

std::vector<std::vector<int>> compute_non_isomorphic_colorings(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<std::unordered_set<std::pair<int, int>>>& edge_orbits,
    const std::unordered_map<std::pair<int, int>, std::vector<int>>& available_colors
);

std::vector<std::vector<int>> apply_edge_automorphisms(
    const std::vector<std::vector<int>>& colorings, 
    const std::vector<std::vector<std::pair<int, int>>>& automorphisms
);

std::pair<int, int> color_to_ports(int color, const std::vector<int>& src_ports, const std::vector<int>& dst_ports);

std::vector<std::unordered_set<std::pair<int, int>>> compute_edge_orbits(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<std::vector<int>>& automorphisms
);

void enumerate_colorings(
    size_t orbit_idx,
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<std::unordered_set<std::pair<int, int>>>& edge_orbits,
    const std::unordered_map<std::pair<int, int>, std::vector<int>>& available_colors,
    std::vector<int>& current_coloring,
    std::vector<std::vector<int>>& all_colorings
);

std::tuple<int, std::vector<int>, std::vector<std::pair<int, int>>> parse_nauty_graph_line(
    const std::string& line, 
    const std::unordered_set<GroupGraph::Node>& node_defs
);

// Function to process nauty output
void process_nauty_output(
    const std::string& line, 
    const std::unordered_set<GroupGraph::Node>& node_defs,
    std::unordered_set<GroupGraph>* graph_basis,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    bool verbose,
    graph* g,
    int* lab,
    int* ptn,
    int* orbits,
    optionblk* options,
    statsblk* stats
);

