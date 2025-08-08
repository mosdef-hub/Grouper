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

#endif // GENERATE_HPP
