#include <iostream>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <queue>
#include <stack>
#include <nlohmann/json.hpp>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <unordered_map>

#include "dataStructures.hpp"
#include "generate.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/MolOps.h>
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

struct NautyEdgeData {
    int num_edges;
    int edges[MAX_EDGES][2];
    int edge_orbits[MAX_EDGES];
    DisjointSet* uf = nullptr;
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
    std::unique_ptr<RDKit::ROMol> rdkitMol;

    // Create a molecule based on whether it's SMILES or SMARTS
    if (isSmarts) {
        rdkitMol.reset(RDKit::SmartsToMol(pattern));
    } else {
        rdkitMol.reset(RDKit::SmilesToMol(pattern));
    }

    if (!rdkitMol) {
        return nullptr; // Return nullptr if molecule creation failed
    }

    RDKit::RWMol processedMol;

    // Process atoms and add to the new molecule
    for (const auto& atom : rdkitMol->atoms()) {
        RDKit::Atom* new_atom = new RDKit::Atom(atom->getAtomicNum()); // Dynamically allocate new atom
        new_atom->setFormalCharge(atom->getFormalCharge()); // Keep charge for valency
        processedMol.addAtom(new_atom); // Add the new atom
        delete new_atom;
    }

    // Add bonds to processed molecule
    for (const auto& bond : rdkitMol->bonds()) {
        processedMol.addBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), bond->getBondType());
    }

    // Return a new ROMol constructed from processedMol
    return std::make_unique<RDKit::ROMol>(processedMol);
}

// =====================================================================
// Process-wide cache of parsed group patterns. toSmiles() and
// toAtomicGraph() are called once per generated graph in the inner loop
// of exhaustive_generate, and each call used to RDKit-parse every group's
// pattern from scratch — three times per call inside the original
// toAtomicGraph. The cache reduces that to one parse per unique
// (pattern, isSmarts) across the whole process. The cached payload is
// plain POD (atomic numbers, formal charges, element symbols, bond
// endpoints + types/orders) so we don't share RDKit ROMol objects across
// threads, only readable data.
// =====================================================================
namespace {

struct GroupMolCacheKey {
    std::string pattern;
    bool isSmarts;
    bool operator==(const GroupMolCacheKey& other) const {
        return isSmarts == other.isSmarts && pattern == other.pattern;
    }
};

struct GroupMolCacheKeyHash {
    std::size_t operator()(const GroupMolCacheKey& k) const {
        std::size_t h = std::hash<std::string>{}(k.pattern);
        h ^= (std::size_t)k.isSmarts * 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct CachedAtom {
    int atomicNumber;
    int formalCharge;
    std::string symbol;
};

struct CachedBond {
    int begin;
    int end;
    RDKit::Bond::BondType bondType;
    double bondOrder;
};

struct CachedMolData {
    std::vector<CachedAtom> atoms;
    std::vector<CachedBond> bonds;
};

std::unordered_map<GroupMolCacheKey, std::shared_ptr<const CachedMolData>, GroupMolCacheKeyHash>
    g_groupMolCache;
std::shared_mutex g_groupMolCacheMutex;

std::shared_ptr<const CachedMolData> getCachedMolData(
    const std::string& pattern, bool isSmarts
) {
    GroupMolCacheKey key{pattern, isSmarts};

    {
        std::shared_lock<std::shared_mutex> rlock(g_groupMolCacheMutex);
        auto it = g_groupMolCache.find(key);
        if (it != g_groupMolCache.end()) return it->second;
    }

    std::unique_lock<std::shared_mutex> wlock(g_groupMolCacheMutex);
    // Re-check after acquiring write lock — another thread may have
    // populated the entry while we waited.
    auto it = g_groupMolCache.find(key);
    if (it != g_groupMolCache.end()) return it->second;

    auto mol = createMol(pattern, isSmarts);
    if (!mol) {
        throw std::runtime_error("Failed to parse pattern: " + pattern);
    }

    auto data = std::make_shared<CachedMolData>();
    data->atoms.reserve(mol->getNumAtoms());
    for (const auto& atom : mol->atoms()) {
        data->atoms.push_back({
            atom->getAtomicNum(),
            atom->getFormalCharge(),
            atom->getSymbol()
        });
    }
    data->bonds.reserve(mol->getNumBonds());
    for (const auto& bond : mol->bonds()) {
        data->bonds.push_back({
            static_cast<int>(bond->getBeginAtomIdx()),
            static_cast<int>(bond->getEndAtomIdx()),
            bond->getBondType(),
            bond->getBondTypeAsDouble()
        });
    }

    g_groupMolCache.emplace(std::move(key), data);
    return data;
}

} // anonymous namespace

// Function to convert rdkit ROMol to AtomGraph
void createAtomGraphFromRDKit(const std::unique_ptr<RDKit::ROMol>& mol, AtomGraph &aG, bool validate=true) {
    const RDKit::PeriodicTable* pt = RDKit::PeriodicTable::getTable();
    for (size_t i=0; i<mol->getNumAtoms(); ++i) {
        const auto& atom = mol->getAtomWithIdx(i);
        int atomicNumber = atom->getAtomicNum();
        // Calculate explicit charge
        int charge = atom->getFormalCharge();
        int valence = pt->getDefaultValence(atomicNumber);
        if (atomicNumber == 16) { // unique handling for sulfur, which can be larger
            valence = 6;
        }
        // Iterate over the bonds and sum the bond orders
        // for (size_t i=0; i<mol->getNumBonds(); i++) {
        //     const auto& bond = mol->getBondWithIdx(i);
        //     double border = bond->getBondTypeAsDouble();
        //     int borderInt = (int)border;
        //     if (bond->getBeginAtomIdx() == atom->getIdx()) {
        //         valence -= borderInt;
        //     }
        //     else if (bond->getEndAtomIdx() == atom->getIdx()) {
        //         valence -= borderInt;
        //     }
        // }
        aG.addNode(atom->getSymbol(), valence+charge);
    }
    for (size_t i=0; i<mol->getNumBonds(); i++) {
        const auto& bond = mol->getBondWithIdx(i);
        double border = bond->getBondTypeAsDouble();
        aG.addEdge(bond->getBeginAtomIdx(), bond->getEndAtomIdx(), border, validate);
    }
}


// Core methods
GroupGraph::GroupGraph()
    : nodes(), edges(), nodetypes(), port_used_bits() {}

GroupGraph::GroupGraph(const GroupGraph& other)
    : nodes(other.nodes), edges(other.edges), nodetypes(other.nodetypes),
      port_used_bits(other.port_used_bits) {}

GroupGraph& GroupGraph::operator=(const GroupGraph& other) {
    if (this != &other) {
        nodes = other.nodes;
        edges = other.edges;
        nodetypes = other.nodetypes;
        port_used_bits = other.port_used_bits;
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

GroupGraph::Group::Group(const std::string& ntype, const std::string& pattern, const std::vector<int>& hubs, const std::string patternType){

    // Error handling
    if (ntype.empty() && pattern.empty()) {
        throw std::invalid_argument("Either SMARTS or type must be provided");
    }
    for (int hub : hubs) {
        if (hub < 0) {
            throw std::invalid_argument("Hub ID must be greater than or equal to 0");
        }
    }
    // Validate the Group inputs are possible
    AtomGraph atomGraph;
    if (patternType == "SMARTS") {
        atomGraph.fromSmarts(pattern);
        validateAtomisticAtomGraph(atomGraph, hubs, pattern, patternType);
    }
    else if (patternType == "SMILES") {
        atomGraph.fromSmiles(pattern);
        validateAtomisticAtomGraph(atomGraph, hubs, pattern, patternType);
    }
    else { // No validation for CG=Coarse Grained validation.
        atomGraph.fromNonAtomic(pattern); // Won't require validation, but parse as SMARTS
        for (int hub : hubs) { // Check that hubs can still be matched
            if (hub > static_cast<int>(atomGraph.nodes.size()) - 1) {
                std::stringstream hubsStr;
                for (const int hubnumber : hubs) {
                    hubsStr << hubnumber; // Append the integer followed by a space
                }
                throw std::invalid_argument("Hub Index ["+ hubsStr.str() + "] of " + std::to_string(hub) + " is greater than the number of atoms in the group");
            }
        }
    }
    this->ntype = ntype;
    this->pattern = pattern;
    this->hubs = hubs;
    this->patternType = patternType;
    std::vector<PortType> ports(hubs.size());
    std::iota(ports.begin(), ports.end(), 0);
    this->ports = ports;
}

void validateAtomisticAtomGraph (const AtomGraph atomGraph, const std::vector<int>& hubs, const std::string pattern, const std::string patternType) {
    // Validate pattern is a valid SMARTS string
    std::unique_ptr<RDKit::ROMol> mol = createMol(pattern, patternType == "SMARTS");
    if (!mol) {
        throw GrouperParseException("Invalid "+ patternType +": " + pattern + " provided");
    }
    // Validate connected molecules
    std::vector<std::vector<int>> moleculesVector;
    RDKit::MolOps::getMolFrags(*mol, moleculesVector);
    // std::vector<boost::shared_ptr<RDKit::ROMol>> moleculesVector = RDKit::MolOps::getMolFrags(*mol);
    if (moleculesVector.size()>1) {
        throw GrouperParseException("Invalid "+ patternType +": " + pattern + " with detached molecules.");
    }
    // Valency Checks
    std::unordered_map<int, int> atomFreeValency;
    for (const auto& [id, node] : atomGraph.nodes) {
        atomFreeValency[id] = node.valency;
        if (atomFreeValency[id] < 0) {
            throw GrouperParseException("Atom "+std::to_string(id)+" has a negative valency after adding the hub for pattern "+ pattern);
        }
    }
    // Validate hubs
    for (int hub : hubs) {
        if (hub > static_cast<int>(atomGraph.nodes.size()) - 1) {
            throw std::invalid_argument("Hub ID "+ std::to_string(hub) +" is greater than the number of atoms in the group");
        }
    }
    for (int hub : hubs) {
        atomFreeValency[hub] -= 1;
        if (atomFreeValency[hub] < 0) {
            throw GrouperParseException("Atom "+std::to_string(hub)+" has a negative valency after adding the hub for pattern "+ pattern);
        }
    }
}

std::vector<int> GroupGraph::Group::hubOrbits() const {

    AtomGraph atomGraph;
    if (patternType == "SMARTS") {
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
    // Cast value
    if (degree > static_cast<int>(hubs.size())) {
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
    if (patternType == "SMARTS") {
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
        int n_nodes = static_cast<int>(modifiedGraph.nodes.size());
        std::vector<std::string> node_colors(modifiedGraph.nodes.size());
        for (int i = 0; i < n_nodes; ++i) node_colors[i] = modifiedGraph.nodes[i].ntype;

        std::unordered_map<std::string, int> color_to_index;
        int color_index = 0;
        for (auto [id, atom] : modifiedGraph.nodes) {
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
    std::string patternType = "SMILES"
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
    if (!ntype.empty() && !pattern.empty()) {
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
        nodes[id] = Group(ntype, pattern, hubs, patternType);
        nodetypes[ntype] = hubs;
        port_used_bits.push_back(0);
    }
    // Case 1: Group type (ntype) is provided
    else if (!ntype.empty() && pattern.empty()) {
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
        nodes[id] = Group(ntype, pattern, hubs, patternType);
        port_used_bits.push_back(0);
    }
    else {
        std::string hubs_str = "[";
        for (int hub : hubs) {
            hubs_str += std::to_string(hub) + " ";
        }
        hubs_str += "]";
        throw std::invalid_argument("Invalid input for add_node -- ntype: " + ntype + " pattern: " + pattern + " hubs: " + hubs_str);
    }
}

void GroupGraph::addNode(Group group) {
    // Error handling
    if (group.ntype.empty()) {
        throw std::invalid_argument("Group type must be provided");
    }
    if (group.pattern.empty()) {
        throw std::invalid_argument("Group pattern must be provided");
    }
    for (int hub : group.hubs) {
        if (hub < 0) {
            throw std::invalid_argument("Hub ID must be greater than or equal to 0");
        }
    }
    for (const auto& entry : nodes) {
        if (entry.second.ntype == group.ntype && entry.second.pattern != group.pattern && entry.second.hubs != group.hubs) {
            throw std::invalid_argument("Group type already exists with different SMARTS/SMILES or hubs");
        }
    }
    int id = nodes.size();
    nodes[id] = group;
    nodetypes[group.ntype] = group.hubs;
    port_used_bits.push_back(0);
}

bool GroupGraph::addEdge(std::tuple<NodeIDType,PortType> fromNodePort, std::tuple<NodeIDType,PortType>toNodePort, double bondOrder, bool strict) {
    NodeIDType from = std::get<0>(fromNodePort);
    PortType fromPort = std::get<1>(fromNodePort);
    NodeIDType to = std::get<0>(toNodePort);
    PortType toPort = std::get<1>(toNodePort);

    if (from == to) {
        if (strict) throw std::invalid_argument("Source and destination nodes are the same");
        return false;
    }

    // Cache the nodes.find iterators so we don't hash the keys again
    // when we read the ports vector below.
    auto fromIt = nodes.find(from);
    if (fromIt == nodes.end()) {
        if (strict) throw std::invalid_argument("Source node does not exist");
        return false;
    }
    auto toIt = nodes.find(to);
    if (toIt == nodes.end()) {
        if (strict) throw std::invalid_argument("Destination node does not exist");
        return false;
    }

    const auto& fromPorts = fromIt->second.ports;
    const auto& toPorts = toIt->second.ports;
    if (std::find(fromPorts.begin(), fromPorts.end(), fromPort) == fromPorts.end()) {
        if (strict) throw std::invalid_argument("Source port does not exist");
        return false;
    }
    if (std::find(toPorts.begin(), toPorts.end(), toPort) == toPorts.end()) {
        if (strict) throw std::invalid_argument("Destination port does not exist");
        return false;
    }

    // Bit-test both endpoints against per-node port_used_bits before
    // committing. If either is already set we bail without mutating
    // anything, so no rollback is needed. port_used_bits is a vector
    // indexed directly by NodeID — no hash lookup, no allocation.
    if (fromPort < 0 || fromPort >= 64 || toPort < 0 || toPort >= 64) {
        if (strict) throw std::invalid_argument("Port index out of range for bitmask (0..63)");
        return false;
    }
    const std::uint64_t fromMask = std::uint64_t{1} << fromPort;
    if (port_used_bits[from] & fromMask) {
        if (strict) throw std::invalid_argument("Source port already in use");
        return false;
    }
    const std::uint64_t toMask = std::uint64_t{1} << toPort;
    if (port_used_bits[to] & toMask) {
        if (strict) throw std::invalid_argument("Destination port already in use");
        return false;
    }

    port_used_bits[from] |= fromMask;
    port_used_bits[to] |= toMask;
    edges.emplace(from, fromPort, to, toPort, bondOrder);
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
    // Preserve the vector's size and storage; just zero the bits.
    // No allocation, just a memset over n_nodes uint64_ts.
    std::fill(port_used_bits.begin(), port_used_bits.end(), 0);
}

void update_edge_orbits(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n) {

    if (!edge_data_ptr) return;

    int num_edges = edge_data_ptr->num_edges;
    int (*nauty_edges)[2] = edge_data_ptr->edges;
    int* edge_orbits = edge_data_ptr->edge_orbits;

    DisjointSet* uf = edge_data_ptr->uf;

    // Build edge index map from normalized edges to their indices
    std::map<std::pair<int, int>, int> edge_index;
    for (int i = 0; i < num_edges; ++i) {
        int u = std::min(nauty_edges[i][0], nauty_edges[i][1]);
        int v = std::max(nauty_edges[i][0], nauty_edges[i][1]);
        edge_index[{u, v}] = i;
    }

    // For each edge, map its endpoints under the current permutation and merge orbits
    for (int i = 0; i < num_edges; ++i) {
        int u = nauty_edges[i][0];
        int v = nauty_edges[i][1];
        int pu = perm[u];
        int pv = perm[v];

        auto permuted_edge = std::make_pair(std::min(pu, pv), std::max(pu, pv));

        auto it = edge_index.find(permuted_edge);
        if (it != edge_index.end()) {
            uf->unite(i, it->second);
        }
    }

    // Assign orbit representatives
    for (int i = 0; i < num_edges; ++i) {
        edge_orbits[i] = uf->find(i);
    }
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<std::vector<int>>> GroupGraph::computeOrbits(
    const std::vector<std::pair<int, int>>& edge_list,
    const std::vector<int>& node_colors,
    graph* g, int* lab, int* ptn, int* orbits, optionblk* options, statsblk* stats
) const {
    int n = nodes.size();
    int m = SETWORDSNEEDED(n);
    setword workspace[160]; // Nauty workspace

    // Allocate a local edge data struct for this thread
    NautyEdgeData edge_data;
    edge_data.num_edges = edge_list.size();
    DisjointSet uf(edge_list.size());
    edge_data.uf = &uf;
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


    // Compute hub orbits for each node
    std::vector<std::vector<int>> hub_orbits_vec(n);
    for (const auto& [id, node] : nodes) {
        hub_orbits_vec[id] = node.hubOrbits();
    }

    // Clear thread-local pointer
    edge_data_ptr = nullptr;

    return {node_orbits, edge_orbits_vec, hub_orbits_vec};
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<std::vector<int>>> GroupGraph::computeOrbits(
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
    DisjointSet uf(edge_list.size());
    edge_data.uf = &uf;
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


    // Compute hub orbits for each node
    std::vector<std::vector<int>> hub_orbits_vec(n);
    for (const auto& [id, node] : nodes) {
        hub_orbits_vec[id] = node.hubOrbits();
    }

    // Clear thread-local pointer
    edge_data_ptr = nullptr;

    return {node_orbits, edge_orbits_vec, hub_orbits_vec};
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
    using AtomIndexMap = std::unordered_map<int, int>;

    auto molecularGraph = std::make_unique<RDKit::RWMol>();
    std::unordered_map<NodeIDType, AtomIndexMap> nodePortToAtomIndex;

    int globalAtomIndex = 0;

    // Single pass: for each group, look up its parsed pattern in the
    // process-wide cache (parsing only happens on first sight of a
    // (pattern, isSmarts) pair), append atoms to the RWMol, record port
    // → global-atom mapping, then append intra-group bonds.
    for (const auto& [nodeID, node] : nodes) {
        bool isSmarts = node.patternType != "SMILES";
        const auto& cached = *getCachedMolData(node.pattern, isSmarts);

        AtomIndexMap& portToGlobal = nodePortToAtomIndex[nodeID];
        const int groupAtomBase = globalAtomIndex;

        for (size_t localIdx = 0; localIdx < cached.atoms.size(); ++localIdx) {
            const auto& a = cached.atoms[localIdx];
            RDKit::Atom newAtom(a.atomicNumber);
            newAtom.setFormalCharge(a.formalCharge);
            molecularGraph->addAtom(&newAtom, true);
            ++globalAtomIndex;
        }

        // Map (this group's port i) → global atom index of the hub atom
        // it sits on. node.hubs[i] is the local atom index inside the
        // pattern; node.ports[i] is the externally-visible port id.
        for (size_t i = 0; i < node.hubs.size(); ++i) {
            portToGlobal[node.ports[i]] = groupAtomBase + node.hubs[i];
        }

        for (const auto& b : cached.bonds) {
            molecularGraph->addBond(
                groupAtomBase + b.begin, groupAtomBase + b.end, b.bondType
            );
        }
    }

    // Inter-subgraph (cross-group) bonds.
    static const std::unordered_map<double, RDKit::Bond::BondType> bondOrderMap = {
        {1.0, RDKit::Bond::BondType::SINGLE},
        {2.0, RDKit::Bond::BondType::DOUBLE},
        {3.0, RDKit::Bond::BondType::TRIPLE},
        {1.5, RDKit::Bond::BondType::AROMATIC},
    };
    for (const auto &[from, fromPort, to, toPort, bondOrder] : edges) {
        int fromAtom = nodePortToAtomIndex.at(from).at(fromPort);
        int toAtom   = nodePortToAtomIndex.at(to).at(toPort);
        molecularGraph->addBond(fromAtom, toAtom, bondOrderMap.at(bondOrder));
    }

    return RDKit::MolToSmiles(*molecularGraph, true, false, -1, true, false);
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
    if (nodes.empty()) {
        throw std::invalid_argument("No nodes in the graph");
    }

    auto atomGraph = std::make_unique<AtomGraph>();
    const RDKit::PeriodicTable* pt = RDKit::PeriodicTable::getTable();

    // (group nodeID, port id) -> global atom index in the assembled
    // AtomGraph. Used to wire up cross-group edges below.
    std::unordered_map<NodeIDType, std::unordered_map<int, int>> nodePortToAtomIndex;

    int atomBase = 0;
    for (const auto& [nodeID, node] : nodes) {
        bool isSmarts = node.patternType != "SMILES";
        const auto& cached = *getCachedMolData(node.pattern, isSmarts);

        for (const auto& a : cached.atoms) {
            int maxValence = pt->getDefaultValence(a.atomicNumber) + a.formalCharge;
            atomGraph->addNode(a.symbol, maxValence);
        }
        for (const auto& b : cached.bonds) {
            atomGraph->addEdge(atomBase + b.begin, atomBase + b.end, b.bondOrder);
        }
        for (size_t i = 0; i < node.hubs.size(); ++i) {
            nodePortToAtomIndex[nodeID][node.ports[i]] = atomBase + node.hubs[i];
        }
        atomBase += static_cast<int>(cached.atoms.size());
    }

    for (const auto& edge : edges) {
        auto [from, fromPort, to, toPort, bondOrder] = edge;
        atomGraph->addEdge(
            nodePortToAtomIndex.at(from).at(fromPort),
            nodePortToAtomIndex.at(to).at(toPort),
            bondOrder
        );
    }

    return atomGraph;
}

std::vector<setword> GroupGraph::canonizeAtomic() const {
    // Fused toAtomicGraph()->canonize(). Builds the same colored
    // atom-level nauty encoding that AtomGraph::canonize produces
    // (one vertex per atom with element-symbol color, plus one
    // auxiliary vertex per bond colored by bond order), but directly
    // from cached group patterns + GroupGraph edges — without
    // allocating an AtomGraph or running its add{Node,Edge} code per
    // call. Profile-driven: AtomGraph::canonize was 37% of total
    // runtime and toAtomicGraph another 7%; fusing them removes the
    // AtomGraph round-trip while preserving the exact same dedup key.
    if (nodes.empty()) {
        return {};
    }

    struct AtomInfo {
        std::string symbol;
    };
    struct BondInfo {
        int begin;
        int end;
        double order;
    };

    // Collect atoms (per-group, in `nodes` iteration order) and bonds
    // (intra-group from cached patterns, then cross-group from edges).
    std::vector<AtomInfo> atomList;
    std::vector<BondInfo> bondList;
    // (group nodeID, port id) -> global atom index, used to wire up
    // cross-group bonds below.
    std::unordered_map<NodeIDType, std::unordered_map<int, int>> nodePortToAtom;

    int atomBase = 0;
    for (const auto& [nodeID, node] : nodes) {
        bool isSmarts = node.patternType != "SMILES";
        const auto& cached = *getCachedMolData(node.pattern, isSmarts);

        for (const auto& a : cached.atoms) {
            atomList.push_back({a.symbol});
        }
        for (const auto& b : cached.bonds) {
            bondList.push_back({atomBase + b.begin, atomBase + b.end, b.bondOrder});
        }
        for (size_t i = 0; i < node.hubs.size(); ++i) {
            nodePortToAtom[nodeID][node.ports[i]] = atomBase + node.hubs[i];
        }
        atomBase += static_cast<int>(cached.atoms.size());
    }

    for (const auto& edge : edges) {
        auto [from, fromPort, to, toPort, bondOrder] = edge;
        bondList.push_back({
            nodePortToAtom.at(from).at(fromPort),
            nodePortToAtom.at(to).at(toPort),
            bondOrder
        });
    }

    // From here this mirrors AtomGraph::canonize: one nauty vertex per
    // atom, one per bond (colored by order), color-string partition for
    // stable order, densenauty + appended canonical color sequence.
    const int numAtoms = static_cast<int>(atomList.size());
    const int numBonds = static_cast<int>(bondList.size());
    const int n = numAtoms + numBonds;
    const int m = SETWORDSNEEDED(n);

    std::vector<setword> g(static_cast<size_t>(m) * n, 0);
    EMPTYGRAPH(g.data(), m, n);

    int bondVertex = numAtoms;
    for (const auto& b : bondList) {
        ADDONEEDGE(g.data(), b.begin, bondVertex, m);
        ADDONEEDGE(g.data(), bondVertex, b.end, m);
        ++bondVertex;
    }

    std::vector<std::string> color_str(n);
    for (int i = 0; i < numAtoms; ++i) {
        color_str[i] = "a:" + atomList[i].symbol;
    }
    for (int i = 0; i < numBonds; ++i) {
        color_str[numAtoms + i] = "b:" + std::to_string(bondList[i].order);
    }

    std::vector<std::pair<std::string, int>> color_sorted;
    color_sorted.reserve(n);
    for (int i = 0; i < n; ++i) color_sorted.emplace_back(color_str[i], i);
    std::sort(color_sorted.begin(), color_sorted.end());

    std::vector<int> lab(n), ptn(n), orbits(n);
    for (int i = 0; i < n; ++i) lab[i] = color_sorted[i].second;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = (color_sorted[i].first == color_sorted[i + 1].first) ? 1 : 0;
    ptn[n - 1] = 0;

    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;

    std::vector<setword> canong(static_cast<size_t>(m) * n, 0);
    densenauty(g.data(), lab.data(), ptn.data(), orbits.data(),
               &options, &stats, m, n, canong.data());

    canong.reserve(canong.size() + n);
    std::hash<std::string> hasher;
    for (int i = 0; i < n; ++i) {
        canong.push_back(static_cast<setword>(hasher(color_str[lab[i]])));
    }
    return canong;
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
            << "\",\n      \"patternType\": \"" << escapeString(node.patternType)
            << "\",\n      \"ports\": [";
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
            std::string patternType = node_data["patternType"].get<std::string>();

            // Create and insert the group
            Group group(ntype, pattern, hubs, patternType);

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
    std::unordered_map<std::tuple<NodeIDType, PortType, NodeIDType, PortType, double>, int> edge_to_nauty;

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
    if (nodes.empty()) {
        return {};
    }

    // Encode the GroupGraph as a colored simple graph for nauty. Three vertex
    // types are interleaved: one nauty vertex per Group, one per (group, port),
    // and one per edge. A per-vertex color (group ntype, port hub-label, edge
    // bond order) is fed to nauty via lab/ptn so the canonical form respects
    // the labels — without that, two graphs differing only in chemistry would
    // canonicalize to the same value.
    std::unordered_map<NodeIDType, int> group_to_nauty;
    std::unordered_map<std::pair<NodeIDType, PortType>, int> port_to_nauty;
    std::unordered_map<std::tuple<NodeIDType, PortType, NodeIDType, PortType, double>, int> edge_to_nauty;

    int nodeIndex = 0;
    for (const auto& [nodeID, group] : nodes) group_to_nauty[nodeID] = nodeIndex++;
    for (const auto& [nodeID, group] : nodes) {
        for (PortType port : group.ports) {
            port_to_nauty[{nodeID, port}] = nodeIndex++;
        }
    }
    for (const auto& edge : edges) edge_to_nauty[edge] = nodeIndex++;

    int n = nodeIndex;
    int m = SETWORDSNEEDED(n);

    std::vector<setword> g(static_cast<size_t>(m) * n, 0);
    EMPTYGRAPH(g.data(), m, n);

    for (const auto& [nodeID, group] : nodes) {
        int gnode = group_to_nauty[nodeID];
        for (PortType port : group.ports) {
            int pnode = port_to_nauty[{nodeID, port}];
            ADDONEEDGE(g.data(), gnode, pnode, m);
        }
    }
    for (const auto& edge : edges) {
        auto [src, srcPort, dst, dstPort, order] = edge;
        int enode = edge_to_nauty[edge];
        int p1 = port_to_nauty[{src, srcPort}];
        int p2 = port_to_nauty[{dst, dstPort}];
        ADDONEEDGE(g.data(), p1, enode, m);
        ADDONEEDGE(g.data(), enode, p2, m);
    }

    // Compute a color *string* for each nauty vertex. Group nodes use the
    // ntype, port nodes use the hub label, edge nodes use the bond order.
    std::vector<std::string> color_str(n);
    for (const auto& [nodeID, group] : nodes) {
        color_str[group_to_nauty[nodeID]] = "g:" + group.ntype;
    }
    for (const auto& [nodeID, group] : nodes) {
        for (size_t i = 0; i < group.ports.size(); ++i) {
            int hub_label = (i < group.hubs.size()) ? group.hubs[i] : -1;
            color_str[port_to_nauty[{nodeID, group.ports[i]}]] =
                "p:" + std::to_string(hub_label);
        }
    }
    for (const auto& edge : edges) {
        color_str[edge_to_nauty[edge]] =
            "e:" + std::to_string(std::get<4>(edge));
    }

    // Sort vertices by the string color so the partition class order is a
    // function of chemistry only — not of the insertion order in `nodes`.
    // This is what makes canonize() reproducible across reorderings.
    std::vector<std::pair<std::string, int>> color_sorted;
    color_sorted.reserve(n);
    for (int i = 0; i < n; ++i) color_sorted.emplace_back(color_str[i], i);
    std::sort(color_sorted.begin(), color_sorted.end());

    std::vector<int> lab(n), ptn(n), orbits(n);
    for (int i = 0; i < n; ++i) lab[i] = color_sorted[i].second;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = (color_sorted[i].first == color_sorted[i + 1].first) ? 1 : 0;
    ptn[n - 1] = 0;

    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;

    std::vector<setword> canong(static_cast<size_t>(m) * n, 0);
    densenauty(g.data(), lab.data(), ptn.data(), orbits.data(),
               &options, &stats, m, n, canong.data());

    // Append the canonical-order color sequence so the returned key encodes
    // both topology and labels. We hash the color *string* (not a per-call
    // integer code) so two graphs with disjoint color sets do not alias.
    canong.reserve(canong.size() + n);
    std::hash<std::string> hasher;
    for (int i = 0; i < n; ++i) {
        canong.push_back(static_cast<setword>(hasher(color_str[lab[i]])));
    }
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
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"P", 3}, {"S", 2}, {"Cl", 1}, {"Br", 1}, {"I", 1}, {"*", 12}
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
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"P", 3}, {"S", 2}, {"Cl", 1}, {"Br", 1}, {"I", 1}, {"*", 12}
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

void AtomGraph::addNode(Atom atom) {
    int id = nodes.size();
    nodes[id] = atom;
}

void AtomGraph::addEdge(NodeIDType src, NodeIDType dst, double order, bool validate) {
    if (nodes.find(src) == nodes.end() || nodes.find(dst) == nodes.end())
    {
        if (nodes.find(src) == nodes.end()) {
            throw std::invalid_argument("Atom " + std::to_string(src) + " does not exist");
        }
        else {
            throw std::invalid_argument("Atom " + std::to_string(dst) + " does not exist");
        }
    }
    if (getFreeValency(src) <= 0 && getFreeValency(dst) <= 0 && validate) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for both nodes");
    }
    if ((getFreeValency(src) <= 0) && validate) {
        throw std::invalid_argument("Adding edge from " + std::to_string(src) + " to " + std::to_string(dst) + " would exceed the valency for the source node");
    }
    if ((getFreeValency(dst) <= 0)  && validate) {
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
                std::unordered_set<NodeIDType> mappedGraphNodes;
                for(const auto& pair : currentMapping) {
                    mappedGraphNodes.insert(pair.second);
                }

                for (const auto& [queryNodeId, count] : queryNeededFreeValency) {
                    NodeIDType graphNodeId = currentMapping[queryNodeId];
                    int external_bonds = 0;
                    for (const auto& edge : this->edges) {
                        if (std::get<0>(edge) == graphNodeId && mappedGraphNodes.find(std::get<1>(edge)) == mappedGraphNodes.end()) {
                            external_bonds++;
                        }
                    }

                    if (external_bonds != count) {
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
 *
 * @tparam pattern a string that will be processed into AtomGraph
 * @return void
 */
void AtomGraph::fromSmarts(const std::string& smarts) {
    nodes.clear();
    edges.clear();

    // Attempt to load via rdkit
    const auto& mol = createMol(smarts, true);
    if (mol) {
        createAtomGraphFromRDKit(mol, *this);
        return;
    };

    // If RDKit fails...
    std::unordered_map<std::string, int> standardElementValency = {
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1},
        {"P", 3}, {"S", 6}, {"Cl", 1}, {"Br", 1}, {"I", 1},
    };

    std::vector<NodeIDType> centralNodeVec;
    std::unordered_map<int, NodeIDType> ringClosures;

    int prevDepth = 0;
    int currentDepth = 0;
    NodeIDType currentNode = 0;
    double bondOrder = 1;

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

            if (static_cast<int>(centralNodeVec.size()) <= currentDepth) {
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

            if (static_cast<int>(centralNodeVec.size()) <= currentDepth) {
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
            throw GrouperParseException("Unsupported character in SMARTS: `" + std::string(1, c) + "` for SMILES `" + smarts + "`");
        }
    }

    if (!ringClosures.empty()) {
        std::cerr << "Unclosed rings detected: ";
        for (const auto& entry : ringClosures) {
            std::cerr << entry.first << " ";
        }
        std::cerr << std::endl;
        throw GrouperParseException("Unclosed ring detected in SMILES string `" + smarts + "`");
    }
    if (currentDepth != 0) {
        throw GrouperParseException("Unmatched parentheses in SMILES string `" + smarts + "`");
    }

}

void AtomGraph::fromSmiles(const std::string& smiles) {
    nodes.clear();
    edges.clear();

    // Attempt to load via rdkit
    const auto& mol = createMol(smiles, false);
    if (mol) {
        createAtomGraphFromRDKit(mol, *this);
        return;
    };

    std::unordered_map<std::string, int> standardElementValency = {
        {"H", 1}, {"B", 3}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1}, {"P", 3}, {"S", 6}, {"Cl", 1}, {"Br", 1}, {"I", 1}, {"c", 4}, {"n", 3}, {"o", 2}, {"s", 2}
    };

    std::stack<NodeIDType> nodeStack; // Stack to handle branching
    std::unordered_map<int, NodeIDType> ringClosures; // Map for ring closure indices
    NodeIDType lastNode = -1;
    double bondOrder = 1; // Default to single bond

    for (size_t i = 0; i < smiles.size(); ++i) {
        char c = smiles[i];

        if (std::isalpha(c)) {
            // Handle atom
            int valency = standardElementValency[std::string(1, c)];
            if (!valency) {
                throw GrouperParseException(
                    "SMILES character `" + std::string(1, c)
                    + "` not in standard element map. Try to add brackets around each element "
                    "for clarity: i.e. 'Li'->'[Li]', or try setting `Grouper.Group` patternType argument to `NONATOMIC`."
                );
            }
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
        } else if (c == '[') {
            // store characters within brackets for next element
            std::string next_elem = "";
            for (size_t j = i+1; j < smiles.size(); ++j) {// iter through next i elements
                char bracketChar = smiles[j] ;
                if (bracketChar == ']') {
                    i = j; //reset to end of parsed brackets
                    break;
                }
                else if (std::isalpha(bracketChar)){
                    next_elem += bracketChar;
                }
            if (next_elem.length() < 1) {
                throw GrouperParseException(
                    "Failed to Parse SMILES of: " + smiles +
                    " at element of " + smiles.substr(i, j)
                );
            }
            }
            int valency = standardElementValency[next_elem];
            addNode(next_elem, valency);
            NodeIDType currentNode = nodes.size() - 1;
            // If there's a previous node, add an edge with the current bond order
            if (lastNode != -1) {
                addEdge(lastNode, currentNode, bondOrder);
            }
            lastNode = currentNode;
            bondOrder = 1; // Reset bond order to single after use
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

void AtomGraph::fromNonAtomic(const std::string& smarts) {
    nodes.clear();
    edges.clear();

    // Attempt to load via rdkit
    const auto& mol = createMol(smarts, true);
    if (mol) {
        createAtomGraphFromRDKit(mol, *this, false);
        return;
    }
    else {
        throw GrouperNotYetImplementedException("NONATOMIC SMARTS must still be parsable by RDKit.");
    }
}


int AtomGraph::getFreeValency(NodeIDType nodeID) const {
    if (nodes.find(nodeID) == nodes.end()) {
        throw std::invalid_argument("Cannot get free valency for non-existent node " + std::to_string(nodeID));
    }
    const Atom& node = nodes.at(nodeID);
    int totalOccupiedValency = 0;
    for (const auto& [src, dst, order] : edges) {
        if (src == nodeID) {
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
    std::unordered_set<std::tuple<NodeIDType, NodeIDType, double>> uniqueEdges;
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
    int n = nodes.size();
    int m = SETWORDSNEEDED(n);

    // Allocate storage for the graph
    std::vector<setword> g(m * n, 0);

    // Initialize the nauty graph
    EMPTYGRAPH(g.data(), m, n);

    // Add all edges (just once per edge, regardless of bond order)
    for (const auto& [src, dst, order] : edges) {
        ADDONEEDGE(g.data(), src, dst, m);
    }

    return g;
}

int AtomGraph::getNodeIndex(int node_id) const {
    int index = 0;
    for (const auto& [id, node] : nodes) {
        if (id == node_id) return index;
        index++;
    }
    return -1; // Not found
}

std::vector<setword> AtomGraph::canonize() const {
    if (nodes.empty()) {
        return {};
    }

    // One nauty vertex per atom; one auxiliary vertex per bond colored
    // by bond order, so single/double/triple/aromatic bonds canonicalize
    // distinctly even though densenauty doesn't natively support edge
    // weights. This is the same trick used in GroupGraph::canonize().
    const int numAtoms = static_cast<int>(nodes.size());
    const int numBonds = static_cast<int>(edges.size());
    const int n = numAtoms + numBonds;
    const int m = SETWORDSNEEDED(n);

    std::vector<setword> g(static_cast<size_t>(m) * n, 0);
    EMPTYGRAPH(g.data(), m, n);

    // Walk edges in unordered_set iteration order. The order is
    // implementation-defined, but we'll feed a colored partition to
    // nauty so the canonical form is invariant to it anyway.
    std::vector<double> bondOrders;
    bondOrders.reserve(numBonds);
    int bondVertex = numAtoms;
    for (const auto& [src, dst, order] : edges) {
        ADDONEEDGE(g.data(), src, bondVertex, m);
        ADDONEEDGE(g.data(), bondVertex, dst, m);
        bondOrders.push_back(order);
        ++bondVertex;
    }

    // Color each nauty vertex with a string. Atoms get "a:<element>",
    // bond vertices get "b:<order>". AtomGraph::addNode assigns IDs in
    // 0..numAtoms-1 by construction (id = nodes.size() at insert time),
    // so atom i in the partition is nodes.at(i).
    std::vector<std::string> color_str(n);
    for (int i = 0; i < numAtoms; ++i) {
        color_str[i] = "a:" + nodes.at(i).ntype;
    }
    for (int i = 0; i < numBonds; ++i) {
        color_str[numAtoms + i] = "b:" + std::to_string(bondOrders[i]);
    }

    // Sort vertices by color string so the partition order is a function
    // of chemistry only — not of insertion order. Without this the
    // canonical form would shift if the same molecule were built atom-
    // by-atom in a different sequence.
    std::vector<std::pair<std::string, int>> color_sorted;
    color_sorted.reserve(n);
    for (int i = 0; i < n; ++i) color_sorted.emplace_back(color_str[i], i);
    std::sort(color_sorted.begin(), color_sorted.end());

    std::vector<int> lab(n), ptn(n), orbits(n);
    for (int i = 0; i < n; ++i) lab[i] = color_sorted[i].second;
    for (int i = 0; i < n - 1; ++i)
        ptn[i] = (color_sorted[i].first == color_sorted[i + 1].first) ? 1 : 0;
    ptn[n - 1] = 0;

    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;

    std::vector<setword> canong(static_cast<size_t>(m) * n, 0);
    densenauty(g.data(), lab.data(), ptn.data(), orbits.data(),
               &options, &stats, m, n, canong.data());

    // Append the canonical-order color sequence so two graphs with
    // disjoint element sets do not alias. We hash the color *string*
    // (not a per-call integer code) so the encoding is globally stable.
    canong.reserve(canong.size() + n);
    std::hash<std::string> hasher;
    for (int i = 0; i < n; ++i) {
        canong.push_back(static_cast<setword>(hasher(color_str[lab[i]])));
    }
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
