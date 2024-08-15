#include <unordered_map>
#include <string>
#include "dataStructures.hpp"


GroupGraph fragment(
    const std::string& smiles, 
    const std::unordered_map<std::string, GroupGraph::Node>& nodeDefs
);