#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <iostream>
#include "process_colored_graphs.cpp"

int main() {

    std::unordered_map<std::string, std::vector<int>> node_types = {
        {"N", {0}},
        {"CO", {0, 1}},
        {"CC", {0, 1, 2, 3}}
    };
    std::unordered_map<int, std::string> int_to_node_type = {
        {0, "N"},
        {1, "CO"},
        {2, "CC"},
        {3, "CC"},
        {4, "CC"}
    };
    const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex = {
        {"N", {{0, 0}}},
        {"CO", {{0, 0}, {1, 0}}},
        {"CC", {{0, 0}, {1, 0}, {2, 1}, {3, 1}}}
    };
    const std::unordered_map<std::string, std::string>& nodeTypeToSmiles = {
        {"N", "N"},
        {"CO", "CO"},
        {"CC", "C=C"}
    };


    std::ifstream input_file("/Users/kieran/projects/molGrouper/molGrouper/cpp_code/vcolg_out.txt");
    if (!input_file.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }
    std::string line;
    while (std::getline(input_file, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }
        auto result = process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, nodeTypeToSmiles, nodeTypePortToIndex, false);
        for (auto it : result.second){
            std::cout << it << std::endl;
        }
    }
    input_file.close();

    return 0;
}

