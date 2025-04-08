#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <queue>
#include <stack>
#include <nlohmann/json.hpp>

#include "dataStructures.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>
#include <RDGeneral/types.h>

#include <nauty/nauty.h>
#include <nauty/naututil.h>

#define MAX_EDGES 100

// Define structures as thread-local using thread_local
thread_local std::vector<setword> g;         // For the nauty graph
thread_local std::vector<int> lab;          // For the label array
thread_local std::vector<int> ptn;          // For the partition array
thread_local std::vector<int> orbits;       // For the orbits array
thread_local DEFAULTOPTIONS_GRAPH(options); // Default nauty options
thread_local statsblk stats;                // Nauty stats structure

struct hash_vector {
    std::size_t operator()(const std::vector<setword>& v) const {
        std::size_t seed = 0;
        for (int i : v) {
            boost::hash_combine(seed, i);
        }
        return seed;
    }
};

struct NautyEdgeData {
    int num_edges;
    int edges[MAX_EDGES][2];
    int edge_orbits[MAX_EDGES];
};

// Each thread gets its own pointer to avoid conflicts
thread_local NautyEdgeData* edge_data_ptr = nullptr;

// Function to initialize the thread-local nauty structures
void initializeNautyStructures(int n) {
    int m = SETWORDSNEEDED(n);

    // Resize the vectors based on the required size
    g.resize(m * n); // The nauty graph requires m * n words
    lab.resize(n);   // Label array size is n
    ptn.resize(n);   // Partition array size is n
    orbits.resize(n); // Orbit array size is n

    // Clear the graph
    std::fill(g.begin(), g.end(), 0);
}

// Function to create a rdkit mol whether smarts or smiles
std::unique_ptr<RDKit::ROMol> createMol(const std::string& pattern, bool isSmarts) {
    return std::unique_ptr<RDKit::ROMol>(isSmarts ? RDKit::SmartsToMol(pattern) : RDKit::SmilesToMol(pattern));
}

// Core methods
GroupGraph::GroupGraph()
    : nodes(), edges(), nodetypes() {}

GroupGraph::GroupGraph(const GroupGraph& other)
    : nodes(other.nodes), edges(other.edges), nodetypes(other.nodetypes) {}

GroupGraph& GroupGraph::operator=(const GroupGraph& other) {
    if (this != &other) {
        nodes = other.nodes;
        edges = other.edges;
        nodetypes = other.nodetypes;
    }
    return *this;
}

bool GroupGraph::Group::operator==(const Group& other) const {
    return ntype == other.ntype &&
           pattern == other.pattern &&
           hubs == other.hubs;
}

bool GroupGraph::Group::operator!=(const Group& other) const {
    return !(*this == other);
}

GroupGraph::Group::Group(const std::string& ntype, const std::string& pattern, const std::vector<int>& hubs, const bool isSmarts){

    // Error handling
    if (ntype.empty() && pattern.empty()) {
        throw std::invalid_argument("Either SMARTS or type must be provided");
    }
    if (hubs.empty()) {
        throw std::invalid_argument("Hubs must be provided");
    }
    for (int hub : hubs) {
        if (hub < 0) {
            throw std::invalid_argument("Hub ID must be greater than or equal to 0");
        }
    }
    // Validate the Group inputs are possible
    AtomGraph atomGraph;
    if (isSmarts) {
        atomGraph.fromSmarts(pattern);
    }
    else {
        atomGraph.fromSmiles(pattern);
    }
    std::unordered_map<int, int> atomFreeValency;
    for (int hub : hubs) {
        if (hub > static_cast<int>(atomGraph.nodes.size()) - 1) {
            throw std::invalid_argument("Hub ID "+ std::to_string(hub) +" is greater than the number of atoms in the group");
        }
    }

    for (const auto& [id, node] : atomGraph.nodes) {
        atomFreeValency[id] = node.valency;
    }
    for (int hub : hubs) {
        atomFreeValency[hub] -= 1;
    }
    for (const auto& [id, valency] : atomFreeValency) {
        if (valency < 0) {
            throw std::invalid_argument("Atom "+std::to_string(id)+" has a negative valency after adding the hub");
        }
    }
    // Validate pattern is a valid SMARTS or SMILES string
    std::unique_ptr<RDKit::ROMol> mol = createMol(pattern, isSmarts);
    if (!mol) {
        throw std::invalid_argument("Invalid SMARTS or SMILES: " + pattern + " provided");
    }
    this->ntype = ntype;
    this->pattern = pattern;
    this->hubs = hubs;
    this->isSmarts = isSmarts;
    std::vector<PortType> ports(hubs.size());
    std::iota(ports.begin(), ports.end(), 0);
    this->ports = ports;
}

std::vector<int> GroupGraph::Group::hubOrbits() const {

    AtomGraph atomGraph;
    if (isSmarts) {
        atomGraph.fromSmarts(pattern);
    }
    else {
        atomGraph.fromSmiles(pattern);
    }
    int n = atomGraph.nodes.size();
    if (n == 0) return {}; // Edge case: no atoms

    // Step 1: Convert group to a nauty-compatible atomic graph
    int m = SETWORDSNEEDED(n);
    std::vector<setword> adj(n * m, 0);

    // Map atoms to nauty node indices
    std::unordered_map<int, int> atom_to_nauty;
    int index = 0;
    for (const auto& [id, node] : atomGraph.nodes) {
        atom_to_nauty[id] = index++;
    }

    // Add edges for atomic connectivity
    for (const auto& [src, dst, order] : atomGraph.edges) {
        int from = atom_to_nauty[src];
        int to = atom_to_nauty[dst];
        ADDONEEDGE(adj.data(), from, to, m);
    }

    // Step 2: Compute atom orbits using nauty
    std::vector<int> lab(n), ptn(n), orbits(n);
    std::vector<setword> canong(n);
    DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    options.getcanon = TRUE;

    densenauty(adj.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats, m, n, canong.data());

    // Step 3: Map orbits back to hubs
    std::vector<int> hub_orbits;
    for (int hub : hubs) {
        if (atom_to_nauty.find(hub) != atom_to_nauty.end()) {
            hub_orbits.push_back(orbits[atom_to_nauty[hub]]);
        }
    }

    return hub_orbits;
}

std::vector<std::vector<int>> GroupGraph::Group::getPossibleAttachments(int degree) const {
    std::vector<std::vector<int>> possible_attachments;

    if (degree > hubs.size()) {
        return possible_attachments; // Not enough ports for the given degree
    }

    // Generate all possible attachment sets using combinations
    std::vector<int> indices(hubs.size());
    std::iota(indices.begin(), indices.end(), 0); // Fill with 0, 1, ..., hubs.size()-1

    std::vector<int> combination;
    std::function<void(int, int)> generate_combinations = [&](int start, int count) {
        if (count == 0) {
            possible_attachments.push_back(combination);
            return;
        }
        for (size_t i = start; i <= indices.size() - count; ++i) {
            combination.push_back(indices[i]);
            generate_combinations(i + 1, count - 1);
            combination.pop_back();
        }
    };

    generate_combinations(0, degree); // Generate all subsets of size `degree`

    // Filter out isomorphic attachment sets
    AtomGraph atomGraph;
    if (isSmarts) {
        atomGraph.fromSmarts(pattern);
    } else {
        atomGraph.fromSmiles(pattern);
    }

    // Store unique canonical forms
    std::unordered_set<std::vector<setword>, hash_vector> unique_canon_attachments;
    std::vector<std::vector<int>> non_isomorphic_attachments;

    for (const auto& attachment : possible_attachments) {
        AtomGraph modifiedGraph = atomGraph;

        // Add dummy attachment atoms to these positions
        for (int port_idx : attachment) {
            modifiedGraph.addNode("I"); // Add dummy attachment atom
            modifiedGraph.addEdge(hubs[port_idx], modifiedGraph.nodes.size() - 1); // Connect to the hub
        }

        // Compute the canonical form using Nauty
        std::vector<setword> nauty_graph = modifiedGraph.toNautyGraph();
        std::vector<int> lab(modifiedGraph.nodes.size()), ptn(modifiedGraph.nodes.size()), orbits(modifiedGraph.nodes.size());
        std::vector<setword> canong(modifiedGraph.nodes.size());

        // Sort nodes by color and initialize `lab` and `ptn`
        int n = modifiedGraph.nodes.size();
        std::vector<std::string> node_colors(modifiedGraph.nodes.size());
        for (int i = 0; i < modifiedGraph.nodes.size(); ++i) node_colors[i] = modifiedGraph.nodes[i].ntype;

        std::unordered_map<std::string, int> color_to_index;
        int color_index = 0;
        for (const auto [id, atom] : modifiedGraph.nodes) {
            if (color_to_index.find(atom.ntype) == color_to_index.end()){
                color_to_index[atom.ntype] = color_index;
                color_index++;
            }
        }
        std::vector<std::pair<int, int>> color_sorted_nodes;
        for (int i = 0; i < n; ++i) color_sorted_nodes.emplace_back(color_to_index[node_colors[i]], i);
        std::sort(color_sorted_nodes.begin(), color_sorted_nodes.end());
        for (int i = 0; i < n; ++i) lab[i] = color_sorted_nodes[i].second;
        for (int i = 0; i < n - 1; ++i) ptn[i] = (color_sorted_nodes[i].first == color_sorted_nodes[i + 1].first) ? 1 : 0;
        ptn[n - 1] = 0;

        DEFAULTOPTIONS_GRAPH(options);
        statsblk stats;
        options.getcanon = TRUE;
        options.defaultptn = FALSE;
        densenauty(nauty_graph.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats, SETWORDSNEEDED(modifiedGraph.nodes.size()), modifiedGraph.nodes.size(), canong.data());

        // Keep only unique (non-isomorphic) attachment sets
        if (unique_canon_attachments.insert(canong).second) {
            non_isomorphic_attachments.push_back(attachment);
        }
    }
    return non_isomorphic_attachments;
}

std::string GroupGraph::Group::toString() const {
    std::ostringstream output;
    output << "Group " <<" (" << ntype << ") (" << pattern << ") ";
    output << ": \n    Ports ";
    for (PortType port : ports) {
        output << port << " ";
    }
    output << "\n    Hubs  ";
    for (NodeIDType hub : hubs) {
        output << hub << " ";
    }
    return output.str();
}

bool GroupGraph::operator==(const GroupGraph& other) const {
    // Check if empty
    if (nodes.empty() && other.nodes.empty()) {
        return true;
    }
    if (nodes.empty() || other.nodes.empty()) {
        return false;
    }

    // Convert GroupGraph to AtomGraph
    std::unique_ptr<AtomGraph> atomGraph1 = this->toAtomicGraph();
    std::unique_ptr<AtomGraph> atomGraph2 = other.toAtomicGraph();

    // Check if the number of nodes and edges are the same
    if (atomGraph1->nodes.size() != atomGraph2->nodes.size()) {
        return false;
    }
    // int edgeCount1 = 0;
    // int edgeCount2 = 0;
    // for (const auto& node : atomGraph1->edges) {
    //     edgeCount1 += node.second.size();
    // }
    // for (const auto& node : atomGraph2->edges) {
    //     edgeCount2 += node.second.size();
    // }
    // if (edgeCount1 != edgeCount2) {
    //     return false;
    // }
    // Check if the number of edges are the same
    if (atomGraph1->edges.size() != atomGraph2->edges.size()) {
        return false;
    }

    // Convert AtomGraph to nauty graph
    int n = atomGraph1->nodes.size(); // Assuming the number of nodes is the same for both graphs
    int m = SETWORDSNEEDED(n);

    // Use std::vector instead of DYNALLSTAT and DYNALLOC
    std::vector<setword> g1(m * n, 0); // Initialize graph 1
    std::vector<setword> g2(m * n, 0); // Initialize graph 2
    std::vector<int> lab1(n), ptn1(n), orbits1(n); // Label, partition, and orbits for graph 1
    std::vector<int> lab2(n), ptn2(n), orbits2(n); // Label, partition, and orbits for graph 2
    std::vector<setword> canong1(m * n, 0); // Canonical form for graph 1
    std::vector<setword> canong2(m * n, 0); // Canonical form for graph 2
    setword workspace[160]; // Workspace for nauty

    // Initialize nauty structures
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    // Convert AtomGraph to nauty graph for g1
    EMPTYGRAPH(g1.data(), m, n);
    // for (const auto& [id, dst_order] : atomGraph1->edges) {
    //     if (!dst_order.empty()) {
    //         for (const auto& [dest,order] : dst_order) {
    //             int from = id;
    //             int to = dest;
    //             ADDONEEDGE(g1.data(), from, to, m);
    //         }
    //     }
    // }
    for (const auto& [src, dst, order]: atomGraph1->edges) {
        ADDONEEDGE(g1.data(), src, dst, m);
    }

    // Convert AtomGraph to nauty graph for g2
    EMPTYGRAPH(g2.data(), m, n);
    // for (const auto& [id, dst_order] : atomGraph2->edges) {
    //     if (!dst_order.empty()) {
    //         for (const auto& [dest,order] : dst_order) {
    //             int from = id;
    //             int to = dest;
    //             ADDONEEDGE(g2.data(), from, to, m);
    //         }
    //     }
    // }
    for (const auto& [src, dst, order]: atomGraph2->edges) {
        ADDONEEDGE(g2.data(), src, dst, m);
    }

    // Call nauty to canonicalize the graphs
    options.getcanon = TRUE;
    nauty(g1.data(), lab1.data(), ptn1.data(), nullptr, orbits1.data(), &options, &stats, workspace, 160, m, n, canong1.data());
    nauty(g2.data(), lab2.data(), ptn2.data(), nullptr, orbits2.data(), &options, &stats, workspace, 160, m, n, canong2.data());

    // Compare the canonical forms to determine isomorphism
    if (memcmp(canong1.data(), canong2.data(), sizeof(setword) * m * n) != 0) {
        return false;
    }

    return true;
}

inline bool operator<(const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>& lhs,
                      const std::tuple<GroupGraph::NodeIDType, GroupGraph::PortType, GroupGraph::NodeIDType, GroupGraph::PortType, unsigned int>& rhs) {
    // Use the built-in tuple comparison operator
    return std::tie(lhs) < std::tie(rhs);
}

// Operating methods
void GroupGraph::addNode(
    std::string ntype = "",
    std::string pattern = "",
    std::vector<NodeIDType> hubs = {},
    bool isSmarts
) {
    /*There are 2 ways you can input a node:
        0. ntype, pattern, hubs
        1. ntype if it already exists
    */

    // Error handling
    if (ntype.empty()) {
        throw std::invalid_argument("Group type must be provided");
    }
    // Case 0: Group type, pattern, and hubs are provided
    if (!ntype.empty() && !pattern.empty() && !hubs.empty()) {
        // Error handling
        for (int hub : hubs) {
            if (hub < 0) {
                throw std::invalid_argument("Hub ID must be greater than or equal to 0");
            }
        }
        for (const auto& entry : nodes) {
            if (entry.second.ntype == ntype && entry.second.pattern != pattern && entry.second.hubs != hubs) {
                throw std::invalid_argument("Group type already exists with different SMARTS/SMILES or hubs");
            }
        }
        int id = nodes.size();
        nodes[id] = Group(ntype, pattern, hubs, isSmarts);
        nodetypes[ntype] = hubs;
    }
    // Case 1: Group type (ntype) is provided
    else if (!ntype.empty() && pattern.empty() && hubs.empty()) {
        if (nodetypes.find(ntype) == nodetypes.end()) {
            throw std::invalid_argument("Group type does not exist yet, please provide SMARTS and hubs");
        }
        int id = nodes.size();
        std::vector<int> hubs = nodetypes[ntype];
        for (const auto& [i, node] : nodes) {
            if (node.ntype == ntype) {
                hubs = node.hubs;
                pattern = node.pattern;
                break;
            }
        }
        nodes[id] = Group(ntype, pattern, hubs, isSmarts);
    }
    else {
        std::string hubs_str = "";
        for (int hub : hubs) {
            hubs_str += std::to_string(hub) + " ";
        }
        throw std::invalid_argument("Invalid input for add_node ntype: " + ntype + " SMARTS/SMILES: " + pattern + " hubs: " + hubs_str);
    }
}

bool GroupGraph::addEdge(std::tuple<NodeIDType,PortType> fromNodePort, std::tuple<NodeIDType,PortType>toNodePort, unsigned int bondOrder, bool verbose) {
    NodeIDType from = std::get<0>(fromNodePort);
    PortType fromPort = std::get<1>(fromNodePort);
    NodeIDType to = std::get<0>(toNodePort);
    PortType toPort = std::get<1>(toNodePort);


    // Error handling
    if (numFreePorts(from) <= 0) {
        throw std::invalid_argument("Source node doesn't have enough ports!");
    }
    if (numFreePorts(to) <= 0) {
        throw std::invalid_argument("Destination node doesn't have enough ports!");
    }
    if (from == to) {
        throw std::invalid_argument("Source and destination nodes are the same");
    }
    if (nodes.find(from) == nodes.end()) {
        throw std::invalid_argument("Source node does not exist");
    }
    if (nodes.find(to) == nodes.end()) {
        throw std::invalid_argument("Destination node does not exist");
    }
    if (std::find(nodes[from].ports.begin(), nodes[from].ports.end(), fromPort) == nodes[from].ports.end()) {
        throw std::invalid_argument("Source port does not exist");
    }
    if (std::find(nodes[to].ports.begin(), nodes[to].ports.end(), toPort) == nodes[to].ports.end()) {
        throw std::invalid_argument("Destination port does not exist");
    }
    const std::tuple<NodeIDType, PortType, NodeIDType, PortType, unsigned int> edge = std::make_tuple(from, fromPort, to, toPort, bondOrder);
    if (std::find(edges.begin(), edges.end(), edge) != edges.end()) {
            throw std::invalid_argument("Edge already exists");
    }
    for (const auto& existingEdge : edges) {
        if (std::get<0>(existingEdge) == from && std::get<1>(existingEdge) == fromPort) {
            throw std::invalid_argument("Source port already in use");
        }
        if (std::get<2>(existingEdge) == to && std::get<3>(existingEdge) == toPort) {
            throw std::invalid_argument("Destination port already in use");
        }
    }

    // Add the edge
    edges.insert(std::make_tuple(from, fromPort, to, toPort, bondOrder));
    // edges.insert(std::make_tuple(to, toPort, from, fromPort, bondOrder)); // Uncomment this line to make the graph undirected
    return true;
}

int GroupGraph::numFreePorts(NodeIDType nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Can't calculate numFreePorts because node does not exist");
    }
    const Group& node = nodes.at(nodeID);
    int occupied_ports = 0;
    for (const auto& edge : edges) {
        if (std::get<0>(edge) == nodeID || std::get<2>(edge) == nodeID) {
            occupied_ports++;
        }
    }
    return node.ports.size() - occupied_ports;
}

bool GroupGraph::isPortFree(NodeIDType nodeID, PortType port) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Can't check if port is free because node does not exist");
    }
    if (std::find(nodes.at(nodeID).ports.begin(), nodes.at(nodeID).ports.end(), port) == nodes.at(nodeID).ports.end()) {
        throw std::invalid_argument("Can't check if port is free because port does not exist");
    }
    for (const auto& edge : edges) {
        if ((std::get<0>(edge) == nodeID && std::get<1>(edge) == port) ||
            (std::get<2>(edge) == nodeID && std::get<3>(edge) == port)) {
            return false;
        }
    }
    return true;
}

void GroupGraph::clearEdges() {
    edges.clear();
}

// Union-Find Data Structure
struct DisjointSet {
    std::vector<int> parent;

    DisjointSet(int n) : parent(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }

    void unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX != rootY) parent[rootY] = rootX;
    }
};

void update_edge_orbits(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n) {

    if (!edge_data_ptr) return; // Ensure pointer is set

    int num_edges = edge_data_ptr->num_edges;
    int (*nauty_edges)[2] = edge_data_ptr->edges;
    int* edge_orbits = edge_data_ptr->edge_orbits;

    DisjointSet uf(num_edges); // Disjoint-set to track edge orbit merging

    for (int i = 0; i < num_edges; i++) {
        int v1 = nauty_edges[i][0], v2 = nauty_edges[i][1];
        int new_v1 = perm[v1], new_v2 = perm[v2];

        for (int j = 0; j < num_edges; j++) {
            if ((nauty_edges[j][0] == new_v1 && nauty_edges[j][1] == new_v2) ||
                (nauty_edges[j][0] == new_v2 && nauty_edges[j][1] == new_v1)) {

                // Merge edge orbits using Union-Find
                uf.unite(i, j);
            }
        }
    }

    // Assign final edge orbit values
    for (int i = 0; i < num_edges; i++) {
        edge_orbits[i] = uf.find(i); // Assign orbit representatives
    }
}

std::pair< std::vector<int>, std::vector<int> > GroupGraph::computeOrbits(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<int>& node_colors,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats
    // int* num_edges, int nauty_edges[][2], int edge_orbits[]
) const {
    int n = nodes.size();
    int m = SETWORDSNEEDED(n);
    setword workspace[160]; // Nauty workspace

    // Allocate a local edge data struct for this thread
    NautyEdgeData edge_data;
    edge_data.num_edges = edge_list.size();
    if (edge_data.num_edges > MAX_EDGES) {
        throw std::runtime_error("Too many edges; increase MAX_EDGES.");
    }

    for (size_t i = 0; i < edge_list.size(); i++) {
        edge_data.edges[i][0] = edge_list[i].first;
        edge_data.edges[i][1] = edge_list[i].second;
        edge_data.edge_orbits[i] = i;  // Initially, each edge is its own orbit
    }

    // Set the thread-local pointer for update_edge_orbits
    edge_data_ptr = &edge_data;

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

    // Configure Nauty options
    options->getcanon = FALSE;
    options->defaultptn = FALSE;
    options->userautomproc = update_edge_orbits;

    // Run Nauty
    densenauty(g, lab, ptn, orbits, options, stats, m, n, workspace);

    // Convert node orbit array to vector
    std::vector<int> node_orbits(n), edge_orbits_vec(edge_data.num_edges);
    for (int i = 0; i < n; ++i) node_orbits[i] = orbits[i];
    for (int i = 0; i < edge_data.num_edges; ++i) edge_orbits_vec[i] = edge_data.edge_orbits[i];


    // Clear thread-local pointer
    edge_data_ptr = nullptr;

    return {node_orbits, edge_orbits_vec};
}

std::pair< std::vector<int>, std::vector<int> > GroupGraph::computeOrbits(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<int>& node_colors
) const {
    int n = nodes.size();
    int m = SETWORDSNEEDED(n);
    setword workspace[160]; // Nauty workspace
    graph g[m * n];         // For the nauty graph
    int lab[n];             // For the label array
    int ptn[n];             // For the partition array
    int orbits[n];          // For the orbits array
    DEFAULTOPTIONS_GRAPH(options); // Default nauty options
    statsblk stats;                // Nauty stats structure


    // Allocate a local edge data struct for this thread
    NautyEdgeData edge_data;
    edge_data.num_edges = edge_list.size();
    if (edge_data.num_edges > MAX_EDGES) {
        throw std::runtime_error("Too many edges; increase MAX_EDGES.");
    }

    for (size_t i = 0; i < edge_list.size(); i++) {
        edge_data.edges[i][0] = edge_list[i].first;
        edge_data.edges[i][1] = edge_list[i].second;
        edge_data.edge_orbits[i] = i;  // Initially, each edge is its own orbit
    }

    // Set the thread-local pointer for update_edge_orbits
    edge_data_ptr = &edge_data;

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

    // Configure Nauty options
    options.getcanon = FALSE;
    options.defaultptn = FALSE;
    options.userautomproc = update_edge_orbits;

    // Run Nauty
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, workspace);

    // Convert node orbit array to vector
    std::vector<int> node_orbits(n), edge_orbits_vec(edge_data.num_edges);
    for (int i = 0; i < n; ++i) node_orbits[i] = orbits[i];
    for (int i = 0; i < edge_data.num_edges; ++i) edge_orbits_vec[i] = edge_data.edge_orbits[i];


    // Clear thread-local pointer
    edge_data_ptr = nullptr;

    return {node_orbits, edge_orbits_vec};
}

// Conversion methods
std::string GroupGraph::printGraph() const {
    std::ostringstream output;
    output << "Nodes:\n";
    for (const auto& entry : nodes) {
        output << "    Group " << entry.first << " (" << entry.second.ntype << ") (" << entry.second.pattern << ") ";
        output<< ": \n        Ports ";
        for (PortType port : entry.second.ports) {
            output << port << " ";
        }
        output << "\n        Hubs  ";
        for (NodeIDType hub : entry.second.hubs) {
            output << hub << " ";
        }
        output << "\n";
    }
    output << "Edges:\n";
    for (const auto& edge : edges) {
        output << "    Edge: " << std::get<0>(edge) << "(" << std::get<1>(edge) << ") -> "
               << std::get<2>(edge) << "(" << std::get<3>(edge) << ")  Order: " << std::get<4>(edge) << "\n";
    }
    return output.str();
}

std::string GroupGraph::toSmiles() const {
    std::unique_ptr<RDKit::RWMol> molecularGraph(new RDKit::RWMol());
    std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
    int atomCount = 0;

    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        std::string pattern = entry.second.pattern;
        std::unique_ptr<RDKit::ROMol> subGraph = createMol(pattern, node.isSmarts);
        nodePortToAtomIndex[std::to_string(nodeID)] = std::unordered_map<int, int>();
        for (size_t i = 0; i < node.ports.size(); ++i) {
            nodePortToAtomIndex[std::to_string(nodeID)][node.ports[i]] = atomCount + node.hubs[i];
        }
        atomCount += subGraph->getNumAtoms();
    }


    int atomId = -1;
    std::unordered_map<std::string, std::unordered_map<int, int>> nodeSubGraphIndicesToMolecularGraphIndices;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)] = std::unordered_map<int, int>();
        std::string pattern = node.pattern;
        std::unique_ptr<RDKit::ROMol> subGraph = createMol(pattern, node.isSmarts);
        for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*atom)->getIdx()] = atomId;
        }
    }


    atomId = -1;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        std::string pattern = node.pattern;
        std::unique_ptr<RDKit::ROMol> subGraph = createMol(pattern, node.isSmarts);
        for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            RDKit::Atom newAtom = **atom;
            molecularGraph->addAtom(&newAtom, true);
        }
        for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
            RDKit::Bond newBond = **bond;
            RDKit::Bond::BondType bondType = (*bond)->getBondType();
            int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getBeginAtomIdx()];
            int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getEndAtomIdx()];
            molecularGraph->addBond(atomIdx1, atomIdx2, bondType);
        }
    }


    for (const auto& edge : edges) {
        NodeIDType from = std::get<0>(edge);
        PortType fromPort = std::get<1>(edge);
        NodeIDType to = std::get<2>(edge);
        PortType toPort = std::get<3>(edge);
        unsigned int bondOrder = std::get<4>(edge);
        int fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
        int toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];

        molecularGraph->addBond(fromAtom, toAtom, static_cast<RDKit::Bond::BondType>(bondOrder));
    }

    return RDKit::MolToSmiles(*molecularGraph, true);
}

std::vector<std::vector<int>> GroupGraph::toEdgeGraph(const std::vector<std::pair<int, int>>& edge_list) const {
    int num_edges = edge_list.size();

    // Initialize edge graph as an adjacency matrix where each element is initially 0
    std::vector<std::vector<int>> edge_graph(num_edges, std::vector<int>(num_edges, 0));

    // Iterate through each pair of edges in the edge list
    for (int i = 0; i < num_edges; i++) {
        for (int j = i + 1; j < num_edges; j++) {
            const std::pair<int, int>& edge_i = edge_list[i];
            const std::pair<int, int>& edge_j = edge_list[j];

            // Check if the two edges share a common node
            if (edge_i.first == edge_j.first || edge_i.first == edge_j.second ||
                edge_i.second == edge_j.first || edge_i.second == edge_j.second) {
                // If they share a node, mark them as connected in the edge graph
                edge_graph[i][j] = 1;
                edge_graph[j][i] = 1;  // Symmetric adjacency matrix
            }
        }
    }

    return edge_graph;
}

std::unordered_map<std::string, int> GroupGraph::toVector() const {
    std::unordered_map<std::string, int> hist(nodetypes.size());
    for (const auto& entry : nodes) {
        const Group& node = entry.second;
        std::string ntype = node.ntype;
        hist[ntype] += 1;
    }
    return hist;
}

std::unique_ptr<AtomGraph> GroupGraph::toAtomicGraph() const {
    auto atomGraph = std::make_unique<AtomGraph>();
    const RDKit::PeriodicTable* pt = RDKit::PeriodicTable::getTable();
    std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
    int atomCount = 0;

    if (nodes.size() == 0) {
        throw std::invalid_argument("No nodes in the graph");
    }

    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        std::string pattern = entry.second.pattern;
        std::unique_ptr<RDKit::ROMol> subGraph = createMol(pattern, node.isSmarts);
        nodePortToAtomIndex[std::to_string(nodeID)] = std::unordered_map<int, int>();
        for (size_t i = 0; i < node.ports.size(); ++i) {
            nodePortToAtomIndex[std::to_string(nodeID)][node.ports[i]] = atomCount + node.hubs[i];
        }
        atomCount += subGraph->getNumAtoms();
    }



    int atomId = -1;
    std::unordered_map<std::string, std::unordered_map<int, int>> nodeSubGraphIndicesToMolecularGraphIndices;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)] = std::unordered_map<int, int>();
        std::string pattern = node.pattern;
        std::unique_ptr<RDKit::ROMol> subGraph = createMol(pattern, node.isSmarts);
        for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*atom)->getIdx()] = atomId;
        }
    }

    // Start building the atom graph from mappings defined above
    atomId = -1;
    for (const auto& entry : nodes) {
        NodeIDType nodeID = entry.first;
        const Group& node = entry.second;
        std::string pattern = node.pattern;
        std::unique_ptr<RDKit::ROMol> subGraph = createMol(pattern, node.isSmarts);
        for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
            atomId++;
            // RDKit::Atom newAtom = **atom;
            // molecularGraph->addAtom(&newAtom, true);
            int atomicNumber = (*atom)->getAtomicNum();
            int maxValence = pt->getDefaultValence(atomicNumber) + (*atom)->getFormalCharge();
            atomGraph->addNode((*atom)->getSymbol(), maxValence);

        }
        for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
            // RDKit::Bond newBond = **bond;
            int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getBeginAtomIdx()];
            int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeID)][(*bond)->getEndAtomIdx()];
            // molecularGraph->addBond(atomIdx1, atomIdx2, newBond.getBondType());
            unsigned int bondOrder = (*bond)->getBondTypeAsDouble();
            atomGraph->addEdge(atomIdx1, atomIdx2, bondOrder);
        }
    }


    for (const auto& edge : edges) {
        NodeIDType from = std::get<0>(edge);
        PortType fromPort = std::get<1>(edge);
        NodeIDType to = std::get<2>(edge);
        PortType toPort = std::get<3>(edge);
        NodeIDType fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
        NodeIDType toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];
        unsigned int bondOrder = std::get<4>(edge);
        atomGraph->addEdge(fromAtom, toAtom, bondOrder);
    }


    return atomGraph;
}

std::string GroupGraph::serialize() const {
    std::ostringstream oss;

    auto escapeString = [](const std::string& input) -> std::string {
        std::ostringstream oss;
        for (char c : input) {
            switch (c) {
                case '\n': oss << "\\n"; break;
                case '\t': oss << "\\t"; break;
                case '\r': oss << "\\r"; break;
                case '\\': oss << "\\\\"; break;
                case '\"': oss << "\\\""; break;
                default: oss << c; break;
            }
        }
        return oss.str();
    };

    oss << "{\n  \"nodes\": [\n";
    for (const auto& pair : nodes) {
        const Group& node = pair.second;
        oss << "    {\n      \"id\": " << pair.first
            << ",\n      \"ntype\": \"" << escapeString(node.ntype)
            << "\",\n      \"pattern\": \"" << escapeString(node.pattern)
            << "\",\n      \"isSmarts\": " << (node.isSmarts ? "true" : "false")
            << ",\n      \"ports\": [";
        for (size_t i = 0; i < node.ports.size(); ++i) {
            if (i > 0) oss << ",";
            oss << node.ports[i];
        }
        oss << "],\n      \"hubs\": [";
        for (size_t i = 0; i < node.hubs.size(); ++i) {
            if (i > 0) oss << ",";
            oss << node.hubs[i];
        }
        oss << "]\n    },\n";
    }
    if (!nodes.empty()) {
        oss.seekp(-2, oss.cur); // Remove last comma and newline
    }
    oss << "\n  ],\n  \"edges\": [\n";
    size_t i = 0;
    for (auto edge : edges) {
        oss << "    ["
            << std::get<0>(edge) << ","
            << std::get<1>(edge) << ","
            << std::get<2>(edge) << ","
            << std::get<3>(edge) << ","
            << std::get<4>(edge) << "]";
        if (i < edges.size() - 1) {
            oss << ",\n";
        }
        i++;
    }
    oss << "\n  ]\n}";
    return oss.str();
}

void GroupGraph::deserialize(const std::string& data) {
    using json = nlohmann::json;
    try {
        // Start with clean state
        nodes.clear();
        edges.clear();
        nodetypes.clear();


        json j = json::parse(data);

        // Pre-allocate space
        const auto& nodes_array = j["nodes"];
        nodes.reserve(nodes_array.size());

        // Process nodes
        for (const auto& node_data : nodes_array) {
            // Extract all data first
            NodeIDType id = node_data["id"].get<NodeIDType>();

            // Construct the group directly with its constructor
            std::string ntype = node_data["ntype"].get<std::string>();
            std::string pattern = node_data["pattern"].get<std::string>();
            std::vector<NodeIDType> hubs = node_data["hubs"].get<std::vector<NodeIDType>>();
            bool isSmarts = node_data["isSmarts"].get<bool>();

            // Create and insert the group
            Group group(ntype, pattern, hubs, isSmarts);

            // If ports were specified, override the default ports
            if (node_data.contains("ports")) {
                group.ports = node_data["ports"].get<std::vector<PortType>>();
            }

            // Use emplace with piecewise construction
            nodes.emplace(std::piecewise_construct,
                         std::forward_as_tuple(id),
                         std::forward_as_tuple(std::move(group)));
        }

        // Process edges
        const auto& edges_array = j["edges"];
        edges.reserve(edges_array.size());

        for (const auto& edge_data : edges_array) {
            edges.insert(
                std::make_tuple(
                    edge_data[0].get<NodeIDType>(),
                    edge_data[1].get<PortType>(),
                    edge_data[2].get<NodeIDType>(),
                    edge_data[3].get<PortType>(),
                    edge_data[4].get<unsigned int>()
                )
            );
        }

    } catch (const json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        throw;
    } catch (const std::exception& e) {
        std::cerr << "Error during deserialization: " << e.what() << std::endl;
        throw;
    }
}

void GroupGraph::toNautyGraph(int* n, int* m, graph** adj) const {
    std::unordered_map<NodeIDType, int> group_to_nauty;
    std::unordered_map<std::pair<NodeIDType, PortType>, int> port_to_nauty;
    std::unordered_map<std::tuple<NodeIDType, PortType, NodeIDType, PortType, unsigned int>, int> edge_to_nauty;

    int nodeIndex = 0;

    // Map GroupGraph nodes to nauty nodes
    for (const auto& [nodeID, group] : nodes) group_to_nauty[nodeID] = nodeIndex++;

    // Map port nodes
    for (const auto& [nodeID, group] : nodes) {
        for (PortType port : group.ports) {
            port_to_nauty[{nodeID, port}] = nodeIndex++;
        }
    }

    // Map edge nodes
    for (const auto& edge : edges) edge_to_nauty[edge] = nodeIndex++;

    *n = nodeIndex;  // Total number of nauty nodes
    *m = SETWORDSNEEDED(*n); // Compute `m` correctly

    // Allocate memory for adj matrix
    *adj = new graph[*n * (*m)]();

    std::fill(*adj, *adj + (*n * (*m)), 0); // Initialize adjacency matrix to 0

    // Build adjacency list
    for (const auto& [nodeID, group] : nodes) {
        int g_node = group_to_nauty[nodeID];
        for (PortType port : group.ports) {
            int p_node = port_to_nauty[{nodeID, port}];
            ADDONEEDGE(*adj, g_node, p_node, *m);
        }
    }

    for (const auto& edge : edges) {
        auto [src, srcPort, dst, dstPort, order] = edge;
        int e_node = edge_to_nauty[edge];
        int p1 = port_to_nauty[{src, srcPort}];
        int p2 = port_to_nauty[{dst, dstPort}];

        ADDONEEDGE(*adj, p1, e_node, *m);
        ADDONEEDGE(*adj, e_node, p2, *m);
    }
}

std::vector<setword> GroupGraph::canonize() const {
    int n, m;
    graph* adj = nullptr; // Initialize pointer

    toNautyGraph(&n, &m, &adj); // Now `adj` is allocated in toNautyGraph

    std::vector<int> lab(n), ptn(n), orbits(n);
    std::vector<setword> canong(n);

    DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    options.getcanon = TRUE;

    densenauty(adj, lab.data(), ptn.data(), orbits.data(), &options, &stats, m, n, canong.data());

    delete[] adj; // Free allocated memory

    return canong;
}


//#############################################################################################################
//#############################################################################################################
//#############################################################################################################

AtomGraph::AtomGraph()
    : nodes(), edges() {}

AtomGraph::AtomGraph(const AtomGraph& other)
    : nodes(other.nodes), edges(other.edges) {}

AtomGraph::Atom::Atom(const std::string& ntype){
    static std::unordered_map<std::string, int> standardElementValency = {
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"P", 3}, {"S", 2}, {"Cl", 1}, {"Br", 1}, {"I", 1}
    };
    this->ntype = ntype;
    if (standardElementValency.count(ntype)) {
        this->valency = standardElementValency[ntype];
    } else { // Error message if passing bad Atom Name
        std::stringstream err_msg;
        err_msg << "Element type '" << ntype << "' does not have a default valency. Valid element types are: ";
        for (const auto& pair : standardElementValency) {
            err_msg << pair.first << " ";
        }
        throw std::invalid_argument(err_msg.str());
    }
}

AtomGraph::Atom::Atom(const std::string& ntype, int valency){
    static std::unordered_map<std::string, int> standardElementValency = {
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"P", 3}, {"S", 2}, {"Cl", 1}, {"Br", 1}, {"I", 1}
    };
    this->ntype = ntype;
    if (valency == -1){
        if (standardElementValency.count(ntype)) {
            this->valency = standardElementValency[ntype];
        } else { // Error message if passing bad Atom Name
            std::stringstream err_msg;
            err_msg << "Element type '" << ntype << "' does not have a default valency. Valid element types are: ";
            for (const auto& pair : standardElementValency) {
                err_msg << pair.first << " ";
            }
            throw std::invalid_argument(err_msg.str());
        }
    } else {
        this->valency = valency;
    }
}

AtomGraph& AtomGraph::operator=(const AtomGraph& other) {
    if (this != &other) {
        nodes = other.nodes;
        edges = other.edges;
    }
    return *this;
}

std::string AtomGraph::Atom::toString() const {
    std::ostringstream output;
    output << "Atom " << " (" << ntype << ") Valency: " << valency;
    return output.str();
}

bool AtomGraph::operator==(const AtomGraph& other) const {
    // Check if the number of nodes and edges are the same
    if (this->nodes.size() != other.nodes.size()) {
        return false;
    }
    if (this->edges.size() != other.edges.size()) {
        return false;
    }

    // Convert AtomGraph to nauty graph
    int n = this->nodes.size(); // Assuming the number of nodes is the same for both graphs
    int m = SETWORDSNEEDED(n);

    // Use std::vector instead of DYNALLSTAT and DYNALLOC
    std::vector<setword> g1(m * n, 0); // Initialize graph 1
    std::vector<setword> g2(m * n, 0); // Initialize graph 2
    std::vector<int> lab1(n), ptn1(n), orbits1(n); // Label, partition, and orbits for graph 1
    std::vector<int> lab2(n), ptn2(n), orbits2(n); // Label, partition, and orbits for graph 2
    std::vector<setword> canong1(m * n, 0); // Canonical form for graph 1
    std::vector<setword> canong2(m * n, 0); // Canonical form for graph 2
    setword workspace[160]; // Workspace for nauty

    // Initialize nauty structures
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    // Convert AtomGraph to nauty graph for g1
    EMPTYGRAPH(g1.data(), m, n);
    for (const auto& [srcNode, dstNode, bondOrder] : this->edges) {
        ADDONEEDGE(g1.data(), srcNode, dstNode, m);
    }

    // Convert AtomGraph to nauty graph for g2
    EMPTYGRAPH(g2.data(), m, n);
    for (const auto& [srcNode, dstNode, bondOrder] : other.edges) {
        ADDONEEDGE(g2.data(), srcNode, dstNode, m);
    }

    // Call nauty to canonicalize the graphs
    options.getcanon = TRUE;
    nauty(g1.data(), lab1.data(), ptn1.data(), nullptr, orbits1.data(), &options, &stats, workspace, 160, m, n, canong1.data());
    nauty(g2.data(), lab2.data(), ptn2.data(), nullptr, orbits2.data(), &options, &stats, workspace, 160, m, n, canong2.data());

    // Compare the canonical forms to determine isomorphism
    if (memcmp(canong1.data(), canong2.data(), sizeof(setword) * m * n) != 0) {
        return false;
    }

    return true;
}

void AtomGraph::addNode(const std::string& type, const int valency) {
    int id = nodes.size();
    nodes[id] = Atom(type, valency);

}

void AtomGraph::addEdge(NodeIDType src, NodeIDType dst, unsigned int order) {
    if (nodes.find(src) == nodes.end() || nodes.find(dst) == nodes.end())
    {
        if (nodes.find(src) == nodes.end()) {
            throw std::invalid_argument("Atom " + std::to_string(src) + " does not exist");
        }
        else {
            throw std::invalid_argument("Atom " + std::to_string(dst) + " does not exist");
        }
    }
    if (getFreeValency(src) <= 0 && getFreeValency(dst) <= 0) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for both nodes");
    }
    if (getFreeValency(src) <= 0) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for the source node");
    }
    if (getFreeValency(dst) <= 0) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for the destination node");
    }
    if (edges.find(std::make_tuple(src, dst, order)) != edges.end() || edges.find(std::make_tuple(dst, src, order)) != edges.end()) {
        throw std::invalid_argument("Edge from " + std::to_string(src) + " to " + std::to_string(dst) + " already exists");
    }
    if (order > 4 || order < 1) {
        throw std::invalid_argument("Bond order of " + std::to_string(order) + " is invalid");
    }
    edges.insert(std::make_tuple(src, dst, order));
    edges.insert(std::make_tuple(dst, src, order));
}

std::vector<std::vector<std::pair<AtomGraph::NodeIDType, AtomGraph::NodeIDType>>> AtomGraph::substructureSearch(const AtomGraph& query, const std::vector<int>& hubs) const {
    /*
        Returns a list of all subgraph isomorphisms between the query graph and this graph
        format is a list of lists of pairs of node ids where (query_node_id, this_node_id) is a match

        @param query AtomGraph to look at isomorphism mapping
        @param hubs vector<ints> which indicate where the available hubs are

        @return vector of vectors of pairs of two NodeIds to map in query
    */
    std::vector<std::vector<std::pair<NodeIDType,NodeIDType>>> matches; // To store all matches
    std::unordered_map<NodeIDType, int> queryNeededFreeValency; // To store the number of hubs for each query node

    // Step 0: Pre-process query hubs
    for (const auto& node : query.nodes) {
        queryNeededFreeValency[node.first] = 0;
    }
    for (const auto& h : hubs) {
        queryNeededFreeValency[h]++;
    }

    // Step 1: Pre-filter nodes in the graph based on query node attributes
    std::unordered_map<NodeIDType, std::vector<NodeIDType>> candidateNodes; // Maps query nodes to possible candidates in the main graph
    for (const auto& queryNodePair : query.nodes) {
        const auto& queryNode = queryNodePair.second;
        const auto& queryID = queryNodePair.first;
        for (const auto &graphNodePair : nodes)
        {
            const auto& graphNode = graphNodePair.second;
            const auto& graphID = graphNodePair.first;

            // Match based on node type and valency
            if (queryNode.ntype == graphNode.ntype && queryNode.valency <= graphNode.valency) {
                candidateNodes[queryID].push_back(graphID);
            }
        }
    }

    // Step 2: Backtracking function to explore mappings
    std::function<void(std::unordered_map<NodeIDType, NodeIDType>&, std::unordered_set<NodeIDType>&)> backtrack =
        [&](std::unordered_map<NodeIDType, NodeIDType>& currentMapping, std::unordered_set<NodeIDType>& usedNodes) {
            // If all query nodes are mapped, validate hubs
            if (currentMapping.size() == query.nodes.size()) {
                // Check if the hubs specified match the query node hubs
                for (const auto& [id, count] : queryNeededFreeValency) {
                    NodeIDType graphNodeid = currentMapping[id];

                    if(this->getFreeValency(graphNodeid) != count) { // Check if number of bonds for query node matches the number of hubs
                        return;
                    }
                }
                for(const auto& [src, dst, order] : query.edges) {
                    NodeIDType graphNodeid = currentMapping[src];
                    auto it = edges.find(std::make_tuple(graphNodeid, currentMapping[dst], order));
                    if (it == edges.end()) { // Node has to be in the graph
                        return;
                    }
                }

                // Add the valid match
                std::vector<std::pair<NodeIDType, NodeIDType>> match;
                for (const auto& mapping : currentMapping) {
                    match.push_back(std::make_pair(mapping.first, mapping.second));
                }
                matches.push_back(match);
                return;
            }

            // Select the next unmapped query node
            NodeIDType nextQueryNode = -1;
            for (const auto& queryNodePair : query.nodes) {
                if (currentMapping.find(queryNodePair.first) == currentMapping.end()) {
                    nextQueryNode = queryNodePair.first;
                    break;
                }
            }

            if (nextQueryNode == -1) return; // No unmapped node found (shouldn't happen)

            // Try each candidate node for the selected query node
            for (NodeIDType candidate : candidateNodes[nextQueryNode]) {
                if (usedNodes.count(candidate)) continue;

                // Check edge consistency
                bool valid = true;
                for (const auto& [queryNeighbor, graphNeighbor] : currentMapping) {
                    // Check if the query graph has an edge between nextQueryNode and queryNeighbor
                    bool edgeExists = false;
                    // for (const auto& edgePair : query.edges.at(nextQueryNode)) {
                    //     if (edgePair.first == queryNeighbor) {
                    //         edgeExists = true;
                    //         break;
                    //     }
                    // }
                    for (const auto& [src, dst, order] : query.edges) {
                        if(src == nextQueryNode && dst == queryNeighbor) {
                            edgeExists = true;
                            break;
                        }
                    }

                    if (edgeExists) {
                        // Verify the corresponding edge exists in the main graph
                        bool graphEdgeExists = false;
                        // for (const auto& edgePair : edges.at(candidate)) {
                        //     if (edgePair.first == graphNeighbor) {
                        //         graphEdgeExists = true;
                        //         break;
                        //     }
                        // }
                        for (const auto& [src, dst, order] : edges) {
                            if(src == candidate && dst == graphNeighbor) {
                                graphEdgeExists = true;
                                break;
                            }
                        }
                        if (!graphEdgeExists) {
                            valid = false;
                            break;
                        }
                    }
                }

                if (!valid) continue;


                // Temporarily map the query node to the candidate
                currentMapping[nextQueryNode] = candidate;
                usedNodes.insert(candidate);

                // Recurse
                backtrack(currentMapping, usedNodes);

                // Backtrack
                currentMapping.erase(nextQueryNode);
                usedNodes.erase(candidate);
            }
        };

    // Step 3: Initialize and start the backtracking process
    std::unordered_map<NodeIDType, NodeIDType> currentMapping; // Maps query node IDs to graph node IDs
    std::unordered_set<NodeIDType> usedNodes; // Tracks already used graph nodes
    backtrack(currentMapping, usedNodes);

    return matches;
}

/**
 * Processing method for creating atom graphs from SMARTS strings.
 *
 *  Currently supported symbols
 *  `(,),[,],;,C,N,O,H,S,F,Br,Cl,I,-,=,#`
 *
 * TODO: Plenty more symbols to support.
 * `R,!,X,ints,*,@`
 * TODO: Add possibility to leverage RDKit parsing alternatively
 * TODO: Handle two letter elements, such as Br, Cl etc.
 *
 * @tparam pattern a string that will be processed into AtomGraph
 * @return void
 */
void AtomGraph::fromSmarts(const std::string& smarts) {
    nodes.clear();
    edges.clear();
    
    // Attempt to load via rdkit
    const RDKit::PeriodicTable* pt = RDKit::PeriodicTable::getTable();
    const auto& mol = createMol(smarts, true);
    if (mol) {
        try {
            for (size_t i=0; i<mol->getNumAtoms(); ++i) {
                const auto& atom = mol->getAtomWithIdx(i);
                int atomicNumber = atom->getAtomicNum();
                int maxValence = pt->getDefaultValence(atomicNumber) + atom->getFormalCharge();
                addNode(atom->getSymbol(), maxValence);
            }
            for (size_t i=0; i<mol->getNumBonds(); i++) {
                const auto& bond = mol->getBondWithIdx(i);
                double border = bond->getBondTypeAsDouble();
                int borderInt = (int)border;
                addEdge(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), borderInt);
            }
            return; // Return early if identified through rdkit and createMol cpp
        }
        catch (const std::exception& e) {
            std::cout << "RDKit parsing failed: " << e.what() << std::endl;
            // Continue to fallback parsing
        }
    };

    // After RDKit fallback fails...
    std::unordered_map<std::string, int> standardElementValency = {
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1},
        {"P", 3}, {"S", 2}, {"Cl", 1}, {"Br", 1}, {"I", 1},
    };

    std::vector<NodeIDType> centralNodeVec;
    std::unordered_map<int, NodeIDType> ringClosures;

    int prevDepth = 0;
    int currentDepth = 0;
    NodeIDType currentNode = 0;
    int bondOrder = 1;

    for (size_t i = 0; i < smarts.length(); ++i) {
        char c = smarts[i];

        if (c == '[') {
            // Bracketed atom, extract until ']'
            size_t end = smarts.find(']', i);
            if (end == std::string::npos) {
                throw std::invalid_argument("Unclosed bracket in SMARTS string: `" + smarts + "`");
            }

            std::string bracketContent = smarts.substr(i + 1, end - i - 1);
            i = end; // advance index

            // Extract atomic symbol and optional charge
            std::string symbol;
            int charge = 0;

            // Simple regex-free parser
            size_t j = 0;
            if (j + 1 < bracketContent.size() && islower(bracketContent[j + 1])) {
                symbol = bracketContent.substr(j, 2);
                j += 2;
            } else {
                symbol = bracketContent.substr(j, 1);
                j += 1;
            }

            // Look for '+' or '-'
            while (j < bracketContent.size()) {
                if (bracketContent[j] == '+') {
                    charge++;
                    j++;
                    while (j < bracketContent.size() && std::isdigit(bracketContent[j])) {
                        charge += bracketContent[j] - '0';
                        j++;
                    }
                } else if (bracketContent[j] == '-') {
                    charge--;
                    j++;
                    while (j < bracketContent.size() && std::isdigit(bracketContent[j])) {
                        charge -= bracketContent[j] - '0';
                        j++;
                    }
                } else {
                    j++;
                }
            }

            // Use standard valence if known, adjusted by charge
            int maxValence = 4;
            if (standardElementValency.count(symbol)) {
                maxValence = standardElementValency[symbol] + charge;
            }

            addNode(symbol, maxValence);
            currentNode = nodes.size() - 1;

            if (centralNodeVec.size() <= currentDepth) {
                centralNodeVec.resize(currentDepth + 1, currentNode);
            }

            if (centralNodeVec.empty()) {
                centralNodeVec.push_back(currentNode);
            } else if (currentDepth <= prevDepth) {
                addEdge(centralNodeVec[currentDepth], currentNode, bondOrder);
                centralNodeVec[currentDepth] = currentNode;
            } else {
                addEdge(centralNodeVec[currentDepth - 1], currentNode, bondOrder);
                centralNodeVec[currentDepth] = currentNode;
            }

            bondOrder = 1;
            prevDepth = currentDepth;
        }
        else if (std::isalpha(c)) {
            // Non-bracket atom, try 2-letter or 1-letter element
            std::string symbol;
            if (i + 1 < smarts.length() && islower(smarts[i + 1])) {
                symbol = smarts.substr(i, 2);
                i++;
            } else {
                symbol = std::string(1, c);
            }

            if (!standardElementValency.count(symbol)) {
                throw std::invalid_argument("Unknown atom type: " + symbol);
            }

            addNode(symbol, standardElementValency[symbol]);
            currentNode = nodes.size() - 1;

            if (centralNodeVec.size() <= currentDepth) {
                centralNodeVec.resize(currentDepth + 1, currentNode);
            }

            if (centralNodeVec.empty()) {
                centralNodeVec.push_back(currentNode);
            } else if (currentDepth <= prevDepth) {
                addEdge(centralNodeVec[currentDepth], currentNode, bondOrder);
                centralNodeVec[currentDepth] = currentNode;
            } else {
                addEdge(centralNodeVec[currentDepth - 1], currentNode, bondOrder);
                centralNodeVec[currentDepth] = currentNode;
            }

            bondOrder = 1;
            prevDepth = currentDepth;
        }
        else if (c == '(') {
            currentDepth++;
        }
        else if (c == ')') {
            currentDepth--;
        }
        else if (std::isdigit(c)) {
            int ringIndex = c - '0';
            if (ringClosures.find(ringIndex) != ringClosures.end()) {
                addEdge(currentNode, ringClosures[ringIndex], bondOrder);
                bondOrder = 1;
                ringClosures.erase(ringIndex);
            } else {
                ringClosures[ringIndex] = currentNode;
            }
        }
        else if (c == '-') {
            bondOrder = 1;
        }
        else if (c == '=') {
            bondOrder = 2;
        }
        else if (c == '#') {
            bondOrder = 3;
        }
        else {
            throw std::invalid_argument("Unsupported character in SMARTS: `" + std::string(1, c) + "` for SMILES `" + smarts + "`");
        }
    }

    if (!ringClosures.empty()) {
        std::cerr << "Unclosed rings detected: ";
        for (const auto& entry : ringClosures) {
            std::cerr << entry.first << " ";
        }
        std::cerr << std::endl;
        throw std::invalid_argument("Unclosed ring detected in SMILES string `" + smarts + "`");
    }
    if (currentDepth != 0) {
        throw std::invalid_argument("Unmatched parentheses in SMILES string `" + smarts + "`");
    }

}

void AtomGraph::fromSmiles(const std::string& smiles) {
    nodes.clear();
    edges.clear();

    std::unordered_map<std::string, int> standardElementValency = {
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"P", 3}, {"S", 2}, {"Cl", 1}, {"Br", 1}, {"I", 1}, {"c", 4}, {"n", 3}, {"o", 2}, {"s", 2}
    };

    std::stack<NodeIDType> nodeStack; // Stack to handle branching
    std::unordered_map<int, NodeIDType> ringClosures; // Map for ring closure indices
    NodeIDType lastNode = -1;
    int bondOrder = 1; // Default to single bond

    for (size_t i = 0; i < smiles.size(); ++i) {
        char c = smiles[i];

        if (std::isalpha(c)) {
            // Handle atom
            int valency = standardElementValency[std::string(1, c)];
            addNode(std::string(1, c), valency);
            NodeIDType currentNode = nodes.size() - 1;

            // If there's a previous node, add an edge with the current bond order
            if (lastNode != -1) {
                addEdge(lastNode, currentNode, bondOrder);
            }

            lastNode = currentNode;
            bondOrder = 1; // Reset bond order to single after use
        } else if (c == '(') {
            // Start a branch, push the current last node onto the stack
            nodeStack.push(lastNode);
        } else if (c == ')') {
            // End a branch, pop the last node from the stack
            if (!nodeStack.empty()) {
                lastNode = nodeStack.top();
                nodeStack.pop();
            } else {
                throw std::runtime_error("Unmatched closing parenthesis in SMILES string.");
            }
        } else if (std::isdigit(c)) {
            // Handle ring closure
            int ringIndex = c - '0';
            if (ringClosures.count(ringIndex)) {
                // Connect the current node to the ring closure with the current bond order
                addEdge(lastNode, ringClosures[ringIndex], bondOrder);
                ringClosures.erase(ringIndex);
            } else {
                // Store the current node as the ring closure point
                ringClosures[ringIndex] = lastNode;
            }
            bondOrder = 1; // Reset bond order to single after use
        } else if (c == '-') {
            // Set bond order to single, this may be useless
            bondOrder = 1;
        } else if (c == '=') {
            // Set bond order to double
            bondOrder = 2;
        } else if (c == '#') {
            // Set bond order to triple
            bondOrder = 3;
        } else {
            // Handle unsupported characters (e.g., invalid SMILES)
            throw std::invalid_argument(
                "Unsupported character in SMILES: `"
                + std::string(1, c)
                + "` for SMILES "
                + std::string(smiles)
            );
        }
    }

    // Basic error checking for unclosed rings
    if (!ringClosures.empty()) {
        std::cerr << "Unclosed rings detected: ";
        for (const auto& entry : ringClosures) {
            std::cerr << entry.first << " ";
        }
        std::cerr << std::endl;
        throw std::invalid_argument(
                "Unclosed ring detected in SMILES string `"
                + std::string(smiles)
                + "`"
            );
    }

}

int AtomGraph::getFreeValency(NodeIDType nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Cannot get free valency for non-existent node " + std::to_string(nodeID));
    }
    const Atom& node = nodes.at(nodeID);
    // if (edges.find(nodeID) == edges.end()) {
    //     return node.valency;
    // }
    // else{
    //     int occupied_electrons = 0;
    //     for (const auto& edge : edges.at(nodeID)) {
    //         occupied_electrons += std::get<1>(edge);
    //     }
    //     return node.valency - occupied_electrons;
    // }
    int totalOccupiedValency = 0;
    for (const auto& [src, dst, order] : edges) {
        if(src == nodeID) {
            totalOccupiedValency += order;
        }
    }
    return node.valency - totalOccupiedValency;
}

std::string AtomGraph::printGraph() const {
    std::ostringstream output;
    output << "Nodes:\n";
    for (const auto& entry : nodes) {
        output << "    Atom " << entry.first << " (" << entry.second.ntype << ")" << " Valency: " << entry.second.valency << "\n";
    }
    output << "Edges (without duplication):\n";
    std::unordered_set<std::tuple<NodeIDType, NodeIDType, unsigned int>> uniqueEdges;
    for (const auto& [src, dst, order] : edges) {
        if (uniqueEdges.find(std::make_tuple(src, dst, order)) == uniqueEdges.end()) {
            output << "    Edge: " << src << " <-> " << dst <<" Order: (" <<order<<")"<<"\n";
            uniqueEdges.insert(std::make_tuple(src, dst, order));
            uniqueEdges.insert(std::make_tuple(dst, src, order));
        }
    }
    return output.str();
}

std::vector<setword> AtomGraph::toNautyGraph() const {
    // Convert AtomGraph to a nauty graph representation
    int n = nodes.size(); // Number of nodes
    int m = SETWORDSNEEDED(n); // Size of one row of the adjacency matrix in setwords

    // Allocate storage for the graph
    std::vector<setword> g(m * n, 0); // Initialize nauty graph (adjacency matrix)

    // Initialize the nauty graph
    EMPTYGRAPH(g.data(), m, n);
    for (const auto& [src, dst, order] : edges) {
        // Add edges based on bond order
        for (int i = 0; i < order; ++i) {
            ADDONEEDGE(g.data(), src, dst, m); // Add edge from 'src' to 'dst'
        }
    }
    // Return the nauty graph representation
    return g;
}

std::vector<setword> AtomGraph::canonize() {
    // Convert AtomGraph to a Nauty graph representation
    std::vector<setword> g = this->toNautyGraph();

    // Prepare vectors and workspace
    int n = nodes.size(); // Number of nodes
    int m = SETWORDSNEEDED(n); // Size of one row of the adjacency matrix in setwords
    std::vector<int> lab(n), ptn(n), orbits(n); // Label, partition, and orbits
    std::vector<setword> canong(n);
    setword workspace[160]; // Workspace for nauty

    // Sort nodes by color and initialize `lab` and `ptn`
    std::vector<std::string> node_colors(n);
    for (int i = 0; i < n; ++i) node_colors[i] = nodes[i].ntype;

    std::unordered_map<std::string, int> color_to_index;
    int color_index = 0;
    for (const auto [id, atom] : nodes) {
        if (color_to_index.find(atom.ntype) == color_to_index.end()){
            color_to_index[atom.ntype] = color_index;
            color_index++;
        }
    }
    std::vector<std::pair<int, int>> color_sorted_nodes;
    for (int i = 0; i < n; ++i) color_sorted_nodes.emplace_back(color_to_index[node_colors[i]], i);
    std::sort(color_sorted_nodes.begin(), color_sorted_nodes.end());
    for (int i = 0; i < n; ++i) lab[i] = color_sorted_nodes[i].second;
    for (int i = 0; i < n - 1; ++i) ptn[i] = (color_sorted_nodes[i].first == color_sorted_nodes[i + 1].first) ? 1 : 0;
    ptn[n - 1] = 0;

    // Set up Nauty options
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = true; // Calculate the canonical labeling
    options.defaultptn = false; // Do not use default partition

    statsblk stats; // Statistics block

    // Run Nauty
    nauty(g.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options, &stats, workspace, 160, m, n, canong.data());

    return canong;
}

std::vector<std::vector<AtomGraph::NodeIDType>> AtomGraph::nodeAut() const {
    int n = nodes.size(); // Number of nodes
    int m = SETWORDSNEEDED(n); // Size of one row of the adjacency matrix in setwords

    // Prepare vectors and workspace
    std::vector<int> lab(n), ptn(n), orbits(n); // Label, partition, and orbits
    std::vector<setword> canong(m * n, 0);      // Canonical form
    setword workspace[160];                    // Workspace for nauty

    // Convert the AtomGraph to a Nauty graph representation
    std::vector<setword> g = this->toNautyGraph();

    // Set up Nauty options
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = true; // Calculate the canonical labeling

    statsblk stats; // Statistics block

    // Run Nauty
    nauty(g.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options, &stats, workspace, 160, m, n, canong.data());

    // Process the automorphism results
    // Orbits represent the equivalence classes of automorphism; map them to a vector of vectors
    std::vector<std::vector<AtomGraph::NodeIDType>> automorphisms(n);
    for (int i = 0; i < n; ++i) {
        automorphisms[orbits[i]].push_back(i);
    }

    // Remove empty entries from automorphisms
    automorphisms.erase(std::remove_if(automorphisms.begin(), automorphisms.end(),
                                       [](const std::vector<int>& v) { return v.empty(); }),
                        automorphisms.end());

    return automorphisms;
}

std::vector<AtomGraph::NodeIDType> AtomGraph::nodeOrbits() const {
    int n = nodes.size(); // Number of nodes
    int m = SETWORDSNEEDED(n); // Size of one row of the adjacency matrix in setwords

    // Prepare vectors and workspace
    std::vector<int> lab(n), ptn(n), orbits(n); // Label, partition, and orbits
    std::vector<setword> canong(m * n, 0);      // Canonical form
    setword workspace[160];                    // Workspace for nauty

    // Convert the AtomGraph to a Nauty graph representation
    std::vector<setword> g = this->toNautyGraph();

    // Set up Nauty options
    static DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = true; // Calculate the canonical labeling

    statsblk stats; // Statistics block

    // Run Nauty
    nauty(g.data(), lab.data(), ptn.data(), nullptr, orbits.data(), &options, &stats, workspace, 160, m, n, canong.data());

    return orbits;
}


// std::string AtomGraph::toSmiles(
// ) const {
// }
