
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

// Hash for std::vector<setword> — used as the key type for canon-based dedup.
struct LocalBasisHash {
    std::size_t operator()(const std::vector<setword>& v) const noexcept {
        // Same mixing pattern as boost::hash_combine, written out so this
        // header doesn't depend on boost.
        std::size_t seed = 0;
        for (auto x : v) {
            seed ^= std::hash<setword>{}(x) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Per-thread accumulator for exhaustive_generate's parallel section.
// Keyed by canonical-form vector so the same molecule produced by
// different vcolg lines collapses into a single entry — replaces the
// previous {per-line canon_set + per-thread unordered_set<GroupGraph>}
// pair, which kept duplicate GroupGraph copies for every cross-line
// occurrence of a molecule. Memory now scales with unique outputs
// rather than with total work performed.
using LocalBasis = std::unordered_map<std::vector<setword>, GroupGraph, LocalBasisHash>;


std::tuple<int, std::vector<int>, std::vector<std::pair<int, int>>> parse_nauty_graph_line(
    const std::string& line,
    const std::unordered_set<GroupGraph::Group>& node_defs
);

// Function to process nauty output
void process_nauty_output(
    const std::string& line,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    LocalBasis* graph_basis,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats // Pass nauty structures
);
