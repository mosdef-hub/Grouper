#ifndef GENERATE_HPP
#define GENERATE_HPP

#include <unordered_set>
#include <string>
#include <vector>
#include <utility>
#include "dataStructures.hpp"

std::unordered_set<GroupGraph> exhaustiveGenerate(
    int n_nodes,
    std::unordered_set<GroupGraph::Group> node_defs,
    int num_procs,
    std::string vcolg_output_file,
    std::unordered_map<std::string, int> positiveConstraints,
    std::unordered_set<std::string> negativeConstraints,
    std::string config_path
);

std::unordered_set<GroupGraph> randomGenerate(
    int n_nodes,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    int num_graphs,
    int num_procs,
    const std::unordered_map<std::string, int>& positiveConstraints,
    const std::unordered_set<std::string>& negativeConstraints
);

std::unordered_set<GroupGraph> randomSample(
    int n_nodes,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    int num_graphs,
    int num_procs,
    const std::unordered_map<std::string, int>& positiveConstraints,
    const std::unordered_set<std::string>& negativeConstraints,
    double extra_edge_prob,
    const std::string& color_strategy,
    int max_attempts,
    long long seed,
    bool show_progress
);

#endif // GENERATE_HPP
