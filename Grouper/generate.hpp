#pragma once

#include <unordered_set>
#include <string>
#include "dataStructures.hpp"

void update_progress(int current, int total);

std::unordered_set<GroupGraph> exhaustiveGenerate(
    int n_nodes, 
    std::unordered_set<GroupGraph::Node> node_defs, 
    std::string nauty_path,
    std::string input_file_path,
    int num_procs = 32,
    std::unordered_map<std::string, int> positiveConstraints = {},
    std::unordered_set<std::string> negativeConstraints = {},
    std::string config_path = "",
    bool verbose = false
);
