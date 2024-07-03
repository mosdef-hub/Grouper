#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <iostream>

#include "GroupGraph.cpp"
#include "process_colored_graphs.cpp"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

void removeDuplicatesGPU(thrust::device_vector<unsigned long>& hash_vector);

// Function to update and display progress bar
void update_progress(int current, int total) {
    float progress = static_cast<float>(current) / total;
    int bar_width = 100;
    
    std::cout << "[";
    int pos = bar_width * progress;
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << current<< "/"<<total << " (" << int(progress * 100.0) << "%)\r";
    std::cout.flush();
}

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


    std::ifstream input_file("/raid6/homes/kierannp/foo/gpufoo/molGrouper/molGrouper/cpp_code/vcolg_out.txt");
    if (!input_file.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

    std::string line;
    int line_count = 0;
    int total_lines = 0;
    while (std::getline(input_file, line)) {
        if (!line.empty()) {  // Optionally skip empty lines
            ++total_lines;
        }
    }

    std::cout << "Processing " << total_lines << " lines..." << std::endl;
    // Reset file position to start from the beginning
    input_file.clear(); // clear any error flags
    input_file.seekg(0, std::ios::beg);

    thrust::host_vector<unsigned long> smiles_basis;
    std::unordered_map<unsigned long, std::string> hash_to_smiles;
    while (std::getline(input_file, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }
        auto result = process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, nodeTypeToSmiles, nodeTypePortToIndex, false);
        
        for (auto it : result.first){
            hash_to_smiles[it.first] = it.second;
        }
        for (auto it : result.second){
            smiles_basis.push_back(it);
        }
        ++line_count;
        update_progress(line_count, total_lines);
    }
    std::cout << std::endl;

    thrust::device_vector<unsigned long> device_hash_vector = smiles_basis;
    removeDuplicatesGPU(device_hash_vector);
    cudaDeviceSynchronize();
    smiles_basis.resize(device_hash_vector.size());
    thrust::copy(device_hash_vector.begin(), device_hash_vector.end(), smiles_basis.begin());



    for (auto it : smiles_basis){
        std::cout << hash_to_smiles[it] << std::endl;
    }

    input_file.close();

    return 0;
}

