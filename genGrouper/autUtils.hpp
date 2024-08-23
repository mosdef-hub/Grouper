#ifndef GENERATE_H
#define GENERATE_H

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <sstream>
#include <iterator>
#include <tuple>


void createAdjMatrix(const std::vector<std::pair<int, int>>& edge_list, std::vector<std::vector<int>>& adj_matrix);

bool areEdgePermutationsEquivalent(const std::vector<std::pair<int, int>>& perm1,
                                   const std::vector<std::pair<int, int>>& perm2);

#endif // GENERATE_H
