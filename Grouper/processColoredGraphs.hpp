
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

#include <nauty/nauty.h>

#include "dataStructures.hpp"

class GroupGraph; // Forward declaration of GroupGraph


std::tuple<int, std::vector<int>, std::vector<std::pair<int, int>>> parse_nauty_graph_line(
    const std::string& line,
    const std::unordered_set<GroupGraph::Group>& node_defs
);

// Function to process nauty output
void process_nauty_output(
    const std::string& line,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    std::unordered_set<GroupGraph>* graph_basis,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats // Pass nauty structures
);
