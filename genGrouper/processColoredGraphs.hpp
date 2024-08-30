

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
    const std::unordered_map<std::string, std::string>& int_to_smiles,
    const std::unordered_map<std::string, std::vector<std::string>>& node_types,
    const std::unordered_map<std::string, std::string>& node_type_to_hub
);

// Function to process nauty output
void process_nauty_output(
    const std::string& line, 
    const std::unordered_set<GroupGraph::Node>& node_defs,
    std::unordered_set<std::string>* graph_basis,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    bool verbose
);

