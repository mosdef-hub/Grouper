#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>


struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ (h2 << 1); // Combine hashes
    }
};

// Helper function to create adjacency matrix
void createAdjMatrix(const std::vector<std::pair<int, int>>& edge_list, std::vector<std::vector<int>>& adj_matrix) {
    for (const auto& edge : edge_list) {
        int u = edge.first;
        int v = edge.second;
        adj_matrix[u][v] = 1;
        adj_matrix[v][u] = 1; // Undirected graph
    }
}

// Helper function to check if two edge permutations are equivalent
bool areEdgePermutationsEquivalent(const std::vector<std::pair<int, int>>& perm1,
                                   const std::vector<std::pair<int, int>>& perm2) {
    if (perm1.size() != perm2.size()) {
        return false;
    }

    std::unordered_set<std::pair<int, int>, pair_hash> set1(perm1.begin(), perm1.end());
    std::unordered_set<std::pair<int, int>, pair_hash> set2(perm2.begin(), perm2.end());

    return set1 == set2;
}