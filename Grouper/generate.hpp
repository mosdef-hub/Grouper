#ifndef GENERATE_HPP
#define GENERATE_HPP

#include <unordered_set>
#include <string>
#include <vector>
#include <utility>
#include "dataStructures.hpp"

extern thread_local size_t isomorphism_checks;

std::pair<std::unordered_set<GroupGraph>, size_t> exhaustiveGenerate(
    int n_nodes,
    std::unordered_set<GroupGraph::Group> node_defs,
    int num_procs,
    std::string vcolg_output_file,
    std::unordered_map<std::string, int> positiveConstraints,
    std::unordered_set<std::string> negativeConstraints,
    std::string config_path
);

std::pair<std::unordered_set<GroupGraph>, size_t> randomGenerate(
    int n_nodes,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    int num_graphs,
    int num_procs,
    const std::unordered_map<std::string, int>& positiveConstraints,
    const std::unordered_set<std::string>& negativeConstraints
);

#endif // GENERATE_HPP
