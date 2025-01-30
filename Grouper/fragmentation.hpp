#include <unordered_map>
#include <string>
#include "dataStructures.hpp"


std::vector<GroupGraph> fragment(
    const std::string& smiles, 
    const std::unordered_set<GroupGraph::Group>& nodeDefs
);

std::vector<int> _findMatchingIndices(const int* orbits, int targetOrbit, int numAtoms);

// std::unordered_map<GroupGraph::Group, std::unordered_set<GroupGraph::Group>> determineNodeComposition(const std::unordered_map<std::string, GroupGraph::Group>& nodeDefs);


struct TupleHash {
    template <typename T1, typename T2, typename T3>
    std::size_t operator ()(const std::tuple<T1, T2, T3>& t) const {
        auto h1 = std::hash<T1>{}(std::get<0>(t)); 
        auto h2 = std::hash<T2>{}(std::get<1>(t)); 
        auto h3 = std::hash<T3>{}(std::get<2>(t)); 
        // Combine the hashes for each element of the tuple
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};