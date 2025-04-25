#pragma once

#include <unordered_set>
#include <string>
#include "dataStructures.hpp"

void update_progress(int current, int total);

std::unordered_set<GroupGraph> exhaustiveGenerate(
    int n_nodes, 
    std::unordered_set<GroupGraph::Group> node_defs, 
    int num_procs = -1,
    std::string vcolg_output_file = "",
    std::unordered_map<std::string, int> positiveConstraints = {},
    std::unordered_set<std::string> negativeConstraints = {},
    std::string config_path = ""
);

std::unordered_set<GroupGraph> randomGenerate(
    int n_nodes, 
    const std::unordered_set<GroupGraph::Group>& node_defs,
    int num_graphs = 100,
    int num_procs = -1,
    const std::unordered_map<std::string, int>& positiveConstraints = {},
    const std::unordered_set<std::string>& negativeConstraints = {}
);