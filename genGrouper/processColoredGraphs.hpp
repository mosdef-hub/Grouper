#ifndef GENERATE_H
#define GENERATE_H

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

// Function to get sensible port combinations
void get_sensible_port_combos(
    const std::vector<std::vector<std::pair<int, int>>>& v, 
    const GroupGraph& gG, 
    const std::vector<std::pair<int, int>>& edge_list, 
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    std::unordered_set<std::string>& smiles_set, 
    const bool verbose
);

// Function to process nauty output
std::unordered_set<std::string> process_nauty_output(
    const std::string& line, 
    const std::unordered_set<GroupGraph::Node>& node_defs,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    bool verbose
);

#endif // GENERATE_H
