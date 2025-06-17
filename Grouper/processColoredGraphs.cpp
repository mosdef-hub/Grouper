#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>
#include <numeric>
#include <exception>
#include <iterator>
#include <set>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <boost/functional/hash.hpp>

#include "dataStructures.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/GraphMol.h>


struct hash_pair {
    std::size_t operator()(const std::pair<int, int>& p) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, p.first);
        boost::hash_combine(seed, p.second);
        return seed;
    }
};

struct hash_vector {
    std::size_t operator()(const std::vector<setword>& v) const {
        std::size_t seed = 0;
        for (const auto& i : v) {
            boost::hash_combine(seed, i);
        }
        return seed;
    }
};
std::tuple< int, std::vector<int>, std::vector<std::pair<int, int>> > parse_nauty_graph_line(const std::string& line, const std::unordered_set<GroupGraph::Group>& node_defs) {
    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        throw std::runtime_error("Invalid nauty output line...");
    }

    std::string node_description_str = line.substr(0, split_pos);
    std::string edge_description_str = line.substr(split_pos + 2);

    std::istringstream node_desc_iss(node_description_str);
    std::vector<std::string> node_description(std::istream_iterator<std::string>{node_desc_iss},
                                              std::istream_iterator<std::string>());

    std::istringstream edge_desc_iss(edge_description_str);
    std::vector<std::string> edge_description(std::istream_iterator<std::string>{edge_desc_iss},
                                              std::istream_iterator<std::string>());

    std::vector<std::pair<int, int>> edge_list;
    for (size_t i = 0; i < edge_description.size(); i += 2) {
        edge_list.emplace_back(std::stoi(edge_description[i]), std::stoi(edge_description[i + 1]));
    }



    int n_vertices = std::stoi(node_description[0]);
    // int n_edges = std::stoi(node_description[1]);
    std::vector<int> colors;
    for (size_t i = 2; i < node_description.size(); ++i) {
        colors.push_back(std::stoi(node_description[i]));
    }

    if (node_defs.size() < static_cast<size_t>(*std::max_element(colors.begin(), colors.end()) + 1)) {
        throw std::runtime_error("Number of nodes in node_defs does not match the number of nodes in the nauty_output_file...");
    }

    return {n_vertices, colors, edge_list};
}

//*****************************************************************************
// Data structure for representing the automorphism group.
// Each element is a permutation on the edge indices.
// In a full application you would obtain this group from nauty.
struct EdgeGroup {
    // perms holds one permutation per group element.
    // Each permutation maps an edge index i (from 0 to num_edges-1)
    // to a new index.
    std::vector<std::vector<int>> perms;
};

//*****************************************************************************
// Lexicographical comparison between two colorings (vectors of ints).
// Returns true if vector 'a' is lexicographically less than vector 'b'.
bool lex_compare(const std::vector<int>& a, const std::vector<int>& b) {
    size_t n = a.size();
    for (size_t i = 0; i < n; ++i) {
        if (a[i] < b[i])
            return true;
        if (a[i] > b[i])
            return false;
    }
    return false; // They are equal.
}

//*****************************************************************************
// Given a coloring and a permutation (on edge indices), produce the permuted coloring.
std::vector<int> apply_permutation(const std::vector<int>& coloring, 
                                   const std::vector<int>& perm) {
    assert(coloring.size() == perm.size());
    std::vector<int> permuted(coloring.size());
    for (size_t i = 0; i < coloring.size(); ++i) {
        permuted[i] = coloring[perm[i]];
    }
    return permuted;
}


bool check_max_bond_not_exceeded(
    const std::vector<std::pair<int,int>>& edge_list,
    const std::vector<int>& colors,
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    const std::unordered_map<int, std::string>& int_to_node_type) {


    std::vector<int> node_bond_count;
    for (std::size_t i = 0; i < colors.size(); ++i) {
        std::string node_type = int_to_node_type.at(colors[i]);
        int n_ports = node_types.at(node_type).size();
        node_bond_count.push_back(n_ports);
    }
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        // int src_color = colors[src];
        // int dst_color = colors[dst];
        node_bond_count[src] -= 1;
        node_bond_count[dst] -= 1;
    }
    for (const auto& count : node_bond_count) {
        if (count < 0) {
            return false;
        }
    }
    return true;
}

// Function to filter orbits based on node-node color pairs
std::vector<std::unordered_set<std::pair<int, int>, hash_pair>> filterOrbits(
    const std::vector<std::pair<int, int>>& edge_list,
    const int* edge_orbits,  // Pointer to array of orbit assignments
    const std::vector<int>& colors)
{
    std::vector<std::unordered_set<std::pair<int, int>, hash_pair>> edge_orbits_filtered; // Filtered edge orbits with
    std::unordered_map<std::pair<int, int>, std::unordered_set<std::pair<int, int>, hash_pair>, hash_pair> color_orbit_map;
    int num_edges = edge_list.size();

    // Map each node-node color pair to its associated edges
    for (int i = 0; i < num_edges; i++) {
        int src_color = colors[edge_list[i].first];
        int dst_color = colors[edge_list[i].second];
        std::pair<int, int> color_pair = std::make_pair(src_color, dst_color);
        color_orbit_map[color_pair].insert(edge_list[i]);
    }

    // Create new orbits based on node-node color pairs
    for (auto const& [color_pair, edges] : color_orbit_map) {
        edge_orbits_filtered.push_back(edges);
    }

    return edge_orbits_filtered;
}

std::vector<std::vector<int>> generateOrbitCombinations(
    const std::unordered_set<std::pair<int, int>, hash_pair>& edge_orbit,
    const std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair>& available_colors)
{
    std::vector<std::vector<int>> all_combinations;
    std::vector<int> all_colors;

    std::vector<int> edge_colors = available_colors.at(*edge_orbit.begin());
    int n_edges = edge_orbit.size();

    // recursive function to generate all possible combinations
    std::function<void(std::vector<int>, int, std::vector<int>, int)> generate_combinations = [&](std::vector<int> elements, std::size_t n_selections, std::vector<int> current_combination = {} , std::size_t i = 0) {
        if (current_combination.size() == n_selections) {
            all_combinations.push_back(current_combination);
            return;
        }
        for (std::size_t j = i; j < elements.size(); ++j) {
            current_combination.push_back(elements[j]);
            generate_combinations(elements, n_selections, current_combination, j + 1);
            current_combination.pop_back();
        }
    };
    generate_combinations(edge_colors, n_edges, {}, 0);
    return all_combinations;
}

/**
 * A helper routine that tests whether a complete edge-coloring is canonical
 * with respect to automorphisms.  Here we use a simple analogue of vcolg's
 * ismax: for each orbit (a set of symmetric edges), we first determine the set
 * of indices (according to a fixed ordering of edge_list) and then require that
 * the color at the minimum index is at least as large as the colors on every
 * other edge in the same orbit.
 *
 * @param coloring   Complete coloring (one color per edge)
 * @param edge_list  The ordered list of edges.
 * @param edge_orbits  The automorphism orbits (each a set of edges).
 * @param edge_index_map A precomputed map from an edge to its index in edge_list.
 *
 * @return true if this coloring passes the maximality test.
 */
bool is_maximal_edge_coloring(
    const std::vector<int>& coloring,
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<std::unordered_set<std::pair<int, int>, hash_pair>>& edge_orbits,
    const std::unordered_map<std::pair<int, int>, int, hash_pair>& edge_index_map)
{
    for (const auto& orbit : edge_orbits) {
        // Collect the indices (in the fixed ordering) of the edges in this orbit.
        std::vector<int> indices;
        for (const auto& edge : orbit) {
            auto it = edge_index_map.find(edge);
            if(it != edge_index_map.end()){
                indices.push_back(it->second);
            }
        }
        // If no index was found for this orbit, skip it.
        if (indices.empty()) continue;
        // Determine the canonical representative: here we simply choose the smallest index.
        int canon_index = *std::min_element(indices.begin(), indices.end());
        // For each edge in the orbit, we require that the color of the canonical edge
        // is at least as high as that of every other edge.
        for (int idx : indices) {
            if (coloring[canon_index] < coloring[idx])
                return false; // Reject: a symmetric edge received a "better" color.
        }
    }
    return true;
}

/**
 * Recursively assigns colors to the edges (like the scan function in vcolg).
 * When all edges are colored, it tests maximality (mimicking trythisone).
 *
 * @param edge_index   The current index in edge_list to assign a color.
 * @param edge_list    The complete ordered list of edges.
 * @param available_colors A mapping that gives for each edge the vector of allowed colors.
 * @param edge_orbits    Automorphism orbits of edges.
 * @param edge_index_map Precomputed map from each edge to its index in edge_list.
 * @param current_coloring The partial coloring (length equals edge_list.size()).
 * @param all_colorings  In/out: accumulates acceptable complete colorings.
 */
void scan_edge_coloring(
    int edge_index,
    const std::vector<std::pair<int, int>>& edge_list,
    const std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair>& available_colors,
    const std::vector<std::unordered_set<std::pair<int, int>, hash_pair>>& edge_orbits,
    const std::unordered_map<std::pair<int, int>, int, hash_pair>& edge_index_map,
    std::vector<int>& current_coloring,
    std::vector<std::vector<int>>& all_colorings)
{
    if (edge_index == static_cast<int>(edge_list.size())) {
        // All edges are colored. Now test maximality.
        if (is_maximal_edge_coloring(current_coloring, edge_list, edge_orbits, edge_index_map)) {
            all_colorings.push_back(current_coloring);
        }
        return;
    }

    // Get the current edge.
    const auto& edge = edge_list[edge_index];

    // Retrieve the allowed colors for this edge.
    auto it = available_colors.find(edge);
    if (it == available_colors.end()) {
        // If no colors are available (this should not happen if input is well-formed),
        // then simply skip this edge.
        scan_edge_coloring(edge_index + 1, edge_list, available_colors,
                           edge_orbits, edge_index_map, current_coloring, all_colorings);
        return;
    }

    const std::vector<int>& colors_for_edge = it->second;
    // Try each allowed color.
    for (int color : colors_for_edge) {
        current_coloring[edge_index] = color;
        scan_edge_coloring(edge_index + 1, edge_list, available_colors,
                           edge_orbits, edge_index_map, current_coloring, all_colorings);
    }
    // Undo the assignment (backtrack)
    current_coloring[edge_index] = -1;
}

/**
 * Generates non-automorphic edge colorings using a recursive backtracking
 * search that is modeled on the logic in vcolg.
 *
 * @param edge_list    The list of edges in the graph.
 * @param edge_orbits  The automorphism orbits of edges.
 * @param available_colors A mapping from each edge to the list of allowed colors.
 * @return A vector of all acceptable complete colorings.
 */
std::vector<std::vector<int>> generateNonAutomorphicEdgeColorings(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<std::unordered_set<std::pair<int, int>, hash_pair>>& edge_orbits,
    const std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair>& available_colors)
{
    std::vector<std::vector<int>> all_colorings;
    // Initialize the current coloring: one entry per edge; -1 means "uncolored."
    std::vector<int> current_coloring(edge_list.size(), -1);

    // Precompute a mapping from an edge to its index in edge_list.
    std::unordered_map<std::pair<int, int>, int, hash_pair> edge_index_map;
    for (size_t i = 0; i < edge_list.size(); ++i) {
        edge_index_map[edge_list[i]] = static_cast<int>(i);
    }

    // Start the recursive backtracking search.
    scan_edge_coloring(0, edge_list, available_colors, edge_orbits, edge_index_map,
                       current_coloring, all_colorings);

    return all_colorings;
}


//*****************************************************************************
// Full group-based maximality test for an edge-coloring.
// For each automorphism (permutation) in the group, apply it to the current
// coloring and compare the result with the original using lexicographical order.
// The coloring is accepted if and only if none of its images are lexicographically
// greater than the coloring itself.
bool full_maximality_test(const std::vector<int>& coloring, const EdgeGroup& group) {
    for (const auto& perm : group.perms) {
        std::vector<int> permuted = apply_permutation(coloring, perm);
        // If the current coloring is lexicographically less than one of its images,
        // then it is not canonical.
        if (lex_compare(coloring, permuted)) {
            return false;
        }
    }
    return true;
}

//*****************************************************************************
// Recursive backtracking search for edge colorings.
// (This routine is analogous to vcolg's 'scan' procedure.)
// 'allowed_colors' is a vector where allowed_colors[i] holds the list of allowed 
// colors for edge i (assuming we fix the ordering of edges as given in edge_list).
// When a full assignment is reached, the routine uses the full group-based maximality test.
void scan_edge_coloring_full(int index,
                             int num_edges,
                             const std::vector<std::vector<int>>& allowed_colors,
                             const EdgeGroup& group,
                             std::vector<int>& current_coloring,
                             std::function<void(const std::vector<int>&)> coloring_callback) {
    if (index == num_edges) {
        // A complete edge coloring has been built.
        if (full_maximality_test(current_coloring, group)) {
            coloring_callback(current_coloring);
        }
        return;
    }
    // Try each allowed color for the edge at 'index'
    for (int color : allowed_colors[index]) {
        current_coloring[index] = color;
        scan_edge_coloring_full(index + 1, num_edges, allowed_colors, group, 
                                current_coloring, coloring_callback);
    }
    // Backtrack: mark current position as uncolored.
    current_coloring[index] = -1;
}

//*****************************************************************************
// Full group-based function to generate non-automorphic edge colorings.
// 'edge_list' is the fixed ordering of edges (each edge is a pair of node indices).
// 'available_colors' maps each edge (as a pair) to its available colors.
// 'group' is the full automorphism group (each element permutes the indices of edge_list).
//
// In this approach, we assign colors one edge at a time and, when a full assignment
// is reached, we check that no automorphic image of that assignment is lexicographically greater.
void generateNonAutomorphicEdgeColorings_Full(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair>& available_colors,
    const EdgeGroup& group,
    std::function<void(const std::vector<int>&)> coloring_callback) {
    int num_edges = edge_list.size();
    // Precompute the list of allowed colors for each edge (using the ordering in edge_list).
    std::vector<std::vector<int>> allowed_colors(num_edges);
    for (int i = 0; i < num_edges; ++i) {
        auto it = available_colors.find(edge_list[i]);
        if (it != available_colors.end()) {
            allowed_colors[i] = it->second;
        } else {
            // If an edge is missing available colors, we could choose to throw an error or simply
            // use an empty list. Here, we throw.
            throw std::runtime_error("Available colors not specified for one or more edges.");
        }
    }

    // Prepare a vector for the current coloring: -1 means uncolored.
    std::vector<int> current_coloring(num_edges, -1);

    // Begin the recursive search.
    scan_edge_coloring_full(0, num_edges, allowed_colors, group, current_coloring, coloring_callback);
}

//***************************************************************************
// Global variable for collecting vertex automorphism generators from nauty.
thread_local std::vector<std::vector<int>> vertex_aut_generators;

//***************************************************************************
// Nauty's userautomproc callback for collecting vertex automorphism generators.
// 'perm' is an array of length n representing an automorphism of the vertex set.
void vertex_autom_callback(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n) {
    // Save the generator. Note: copy only the first n elements.
    std::vector<int> gen(perm, perm + n);
    vertex_aut_generators.push_back(gen);
}
//***************************************************************************
// Given an edge list and a vertex permutation (generator), compute the induced
// permutation on the edge set. Edges are assumed to be undirected and normalized
// as (min(u,v), max(u,v)).
std::vector<int> induceEdgePermutation(
    const std::vector<std::pair<int,int>>& edge_list,
    const std::vector<int>& vertex_perm)
{
    int num_edges = static_cast<int>(edge_list.size());
    std::vector<int> induced(num_edges, -1);

    // Build a map from normalized edge to its index.
    std::map<std::pair<int,int>, int> edgeIndex;
    for (int i = 0; i < num_edges; ++i) {
        int u = edge_list[i].first;
        int v = edge_list[i].second;
        int a = std::min(u, v), b = std::max(u, v);
        edgeIndex[{a, b}] = i;
    }
    // For each edge, compute its image under the vertex permutation.
    for (int i = 0; i < num_edges; ++i) {
        int u = edge_list[i].first;
        int v = edge_list[i].second;
        int nu = vertex_perm[u];
        int nv = vertex_perm[v];
        int a = std::min(nu, nv), b = std::max(nu, nv);
        auto it = edgeIndex.find({a, b});
        if (it == edgeIndex.end())
            throw std::runtime_error("Edge image not found in edge index map.");
        induced[i] = it->second;
    }
    return induced;
}

//***************************************************************************
// Function to obtain the generators of the edge automorphism group.
// 'edge_list' is the fixed ordering of edges (each an undirected edge given as (u,v)).
// 'g' is a nauty graph on the vertices of your graph (constructed prior to calling).
// 'n' is the number of vertices; 'm' is the number of words per vertex (typically from SETWORDSNEEDED(n)).
// The graph 'g' is used as input to nauty to compute the vertex automorphism group.
// After running nauty, the vertex generators are used to induce edge permutations.
EdgeGroup obtainEdgeAutomorphismGenerators(
    const std::vector<std::pair<int,int>>& edge_list,
    const std::vector<int>& node_colors,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats
)
{
    int n = node_colors.size();
    int m = SETWORDSNEEDED(n);
    setword workspace[160]; // Nauty workspace

    // Initialize graph structure
    EMPTYGRAPH(g, m, n);
    for (const auto& edge : edge_list) ADDONEEDGE(g, edge.first, edge.second, m);

    // Sort nodes by color and initialize `lab` and `ptn`
    std::vector<std::pair<int, int>> color_sorted_nodes;
    for (int i = 0; i < n; ++i) color_sorted_nodes.emplace_back(node_colors[i], i);
    std::sort(color_sorted_nodes.begin(), color_sorted_nodes.end());

    for (int i = 0; i < n; ++i) lab[i] = color_sorted_nodes[i].second;
    for (int i = 0; i < n - 1; ++i) ptn[i] = (color_sorted_nodes[i].first == color_sorted_nodes[i + 1].first) ? 1 : 0;
    ptn[n - 1] = 0;

    // We want to collect the vertex automorphisms.
    options->userautomproc = vertex_autom_callback;
    options->getcanon = FALSE;
    options->defaultptn = FALSE;

    // Run nauty on the vertex graph.
    densenauty(g, lab, ptn, orbits, options, stats, m, n, workspace);
    
    // At this point, vertex_aut_generators has been filled via the callback.
    // For debugging or checking, you can examine vertex_aut_generators.size().

    // Build the edge automorphism group from the vertex generators.
    EdgeGroup edgeGroup;
    for (const auto& vertexGen : vertex_aut_generators) {
        // Compute the induced permutation on edges.
        std::vector<int> edgePerm = induceEdgePermutation(edge_list, vertexGen);
        edgeGroup.perms.push_back(edgePerm);
    }
    
    // Clear the temporary storage for future calls.
    vertex_aut_generators.clear();
    
    return edgeGroup;
}

//*****************************************************************************
void processColoring(
    const std::vector<int>& coloring,
    const std::vector<std::pair<int, int>>& edge_list,
    const std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, hash_pair>& color_to_port_pair,
    const std::unordered_map<int, std::string>& int_to_node_type,
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    const std::vector<int>& colors,
    std::unordered_set<std::string>& canon_set,
    std::unordered_set<GroupGraph>* graph_basis,
    GroupGraph& gG
) {
    gG.clearEdges();
    size_t edge_index = 0;
    bool all_edges_added = true;
    
    for (const auto& edge : edge_list) {
        // Ensure canonical order for undirected edge:
        int s = edge.first;
        int t = edge.second;

        int color = coloring[edge_index++];
        std::pair<int,int> colorPort = color_to_port_pair.at(edge)[color];
        int sPort = colorPort.first;
        int tPort = colorPort.second;
        bool added;

        added = gG.addEdge({s, sPort}, {t, tPort}, 1, false);
        if (!added)
            added = gG.addEdge({s, tPort}, {t, sPort}, 1, false);

        if (!added) {
            all_edges_added = false;
            break;
        }
    }
    
    if (!all_edges_added) return;
    
    std::string smiles = gG.toSmiles();
    if (canon_set.insert(smiles).second) {
        graph_basis->insert(gG);
    }
}
//*****************************************************************************

void process_nauty_output(
    const std::string& line,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    std::unordered_set<GroupGraph>* graph_basis,
    const std::unordered_map<std::string, int> positiveConstraints,
    const std::unordered_set<std::string> negativeConstraints,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats // Pass nauty structures
) {
    // logMemoryUsage("Start process_nauty_output");
    // Process the nauty output line
    auto [n_vertices, colors, edge_list] = parse_nauty_graph_line(line, node_defs);

    // Actual function starts here
    std::unordered_map<int, std::string> int_to_node_type;
    std::unordered_map<std::string, std::vector<int>> node_types;
    std::unordered_map<std::string, std::vector<int>> node_type_to_hub;
    std::unordered_map<int, std::string> int_to_pattern;
    std::unordered_map<std::string, std::string> node_type_to_pattern_type;
    std::unordered_map<std::string, std::string> type_to_pattern;
    std::vector<GroupGraph> group_graphs_list;
    // std::unordered_set<std::vector<setword>, hash_vector> canon_set;
    std::unordered_set<std::string> canon_set;

    // Create necessary maps
    for (const auto& node : node_defs) {
        node_types[node.ntype] = node.ports;
        type_to_pattern[node.ntype] = node.pattern;
        node_type_to_pattern_type[node.ntype] = node.patternType;
    }
    for (const auto& node: node_defs) {
        for (const auto& h : node.hubs) {
            node_type_to_hub[node.ntype].push_back(h);
        }
    }
    int i = 0;
    for (const auto& [node_type, ports] : node_types) {
        int_to_node_type[i] = node_type;
        int_to_pattern[i] = type_to_pattern[node_type];
        i++;
    }
    // Create a histogram of node types for positive constraints
    if (positiveConstraints.size() > 0) {
        std::unordered_map<std::string, int> node_hist;
        for (const auto& node : node_defs) {
            node_hist[node.ntype] = 0;
        }
        for (const auto& c : colors) {
            node_hist[int_to_node_type.at(c)] += 1;
        }

        // Check if the number of nodes of each type is less than the number listed in the positive constraints
        for (const auto& [node_type, count] : node_hist) {
            if (positiveConstraints.find(node_type) == positiveConstraints.end()) {
                continue;
            }
            if (count < positiveConstraints.at(node_type)) {
                return;
            }
        }
    }
    // Check if any node has too many bonds
    if (!check_max_bond_not_exceeded(edge_list, colors, node_types, int_to_node_type)) {
        return;
    }
    // Add nodes to the graph to create node colored graph
    GroupGraph gG;
    for (int i = 0; i < n_vertices; ++i) {
        gG.addNode(
            int_to_node_type.at(colors[i]),
            int_to_pattern.at(colors[i]),
            node_type_to_hub.at(int_to_node_type.at(colors[i])),
            node_type_to_pattern_type.at(int_to_node_type.at(colors[i]))
        );
    }
    // Get degree of each node
    std::unordered_map<int, int> node_degree;
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        node_degree[src]++;
        node_degree[dst]++;
    }
    // Compute port representatives for each node
    std::unordered_map<int, std::vector<std::vector<int>>> node_port_representatives;
    std::unordered_map<std::pair<int, int>, std::vector<int>, hash_pair> possible_edge_colors;
    std::unordered_map<std::pair<int, int>,std::vector<std::pair<int, int>>,hash_pair> color_to_port_pair;
    GroupGraph::Group tmp_g;
    for (const auto& [node, degree] : node_degree) {
        tmp_g.ntype = int_to_node_type.at(colors[node]);
        tmp_g.hubs = node_type_to_hub.at(int_to_node_type.at(colors[node]));
        tmp_g.pattern = int_to_pattern.at(colors[node]);
        node_port_representatives[node] = tmp_g.getPossibleAttachments(degree);
    }

    // lambda function to flatten the port representatives
    auto flatten_ports = [](const std::vector<std::vector<int>>& reps) {
        std::unordered_set<int> ports;
        for (const auto& group : reps) {
            ports.insert(group.begin(), group.end());
        }
        return std::vector<int>(ports.begin(), ports.end());
    };
    // Compute possible edge colors
    for (const auto& [src, dst] : edge_list) {
        const auto& src_reps = node_port_representatives.at(src);
        const auto& dst_reps = node_port_representatives.at(dst);
    
        std::vector<int> src_ports = flatten_ports(src_reps);
        std::vector<int> dst_ports = flatten_ports(dst_reps);
    
        std::vector<int> e_colors;
        std::vector<std::pair<int, int>> port_pairs;
        int color_index = 0;
    
        for (int i : src_ports) {
            for (int j : dst_ports) {
                e_colors.push_back(color_index++);
                port_pairs.push_back({i, j});
            }
        }
        possible_edge_colors[{src, dst}] = e_colors;
        color_to_port_pair[{src, dst}] = port_pairs;
    }

    // Generate edge automorphism group
    EdgeGroup edgeGroup = obtainEdgeAutomorphismGenerators(
        edge_list, colors, 
        g, lab, ptn, orbits, options, stats // Pass nauty structures
    );
    // logMemoryUsage("After obtainEdgeAutomorphismGenerators");
    // std::vector<std::vector<int>> unique_colorings = generateNonAutomorphicEdgeColorings_Full(edge_list, possible_edge_colors, edgeGroup);
    generateNonAutomorphicEdgeColorings_Full(
        edge_list, possible_edge_colors, edgeGroup,
        [&](const std::vector<int>& coloring) {
            processColoring(
                coloring,
                edge_list,
                color_to_port_pair,
                int_to_node_type,
                node_types,
                colors,
                canon_set,
                graph_basis,
                gG
            );
        }
    );
}
