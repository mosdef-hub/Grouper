#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <omp.h>
#include <stdio.h>

#include "dataStructures.hpp"
#include "processColoredGraphs.hpp"

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
    std::cout << "] " << current << "/" << total << " (" << int(progress * 100.0) << "%)\r";
    std::cout.flush();
}

std::unordered_set<std::string> exhaustiveGenerate(
    int n_nodes, 
    std::unordered_set<GroupGraph::Node> node_defs, 
    std::string input_file_path,
    int num_procs = 32,
    bool verbose = false
) {
    std::ifstream input_file(input_file_path);
    if (!input_file.is_open()) {
        throw std::runtime_error("Error opening input file...");
    }

    std::string line;
    std::vector<std::string> lines;
    int total_lines = 0;

    // Read all lines into a vector
    while (std::getline(input_file, line)) {
        if (!line.empty()) {  // Optionally skip empty lines
            lines.push_back(line);
            ++total_lines;
        }
    }

    std::cout << "Processing " << total_lines << " lines..." << std::endl;

    input_file.close(); // Close the input file as it's no longer needed

    std::unordered_set<std::string> smiles_basis;

    omp_set_num_threads(num_procs);      // Set the number of threads to match

    std::cout<< "Using "<<num_procs << " processors" << std::endl;

    #pragma omp parallel
    {

        // int thread_id = omp_get_thread_num();
        // std::cout<< "Thread " << thread_id << " started" << std::endl;
        std::cout.flush();

        // Thread-local data structures
        std::unordered_set<std::string> local_smiles_basis;

        #pragma omp for schedule(dynamic) nowait
        for (int i = 0; i < total_lines; ++i) {
            // std::cout<< "Thread " << thread_id << " processing line " << i << std::endl;
            std::unordered_set<std::string> result = process_nauty_output(lines[i], node_defs, verbose);

            for (auto it : result) {
                local_smiles_basis.insert(it);
            }

            #pragma omp critical
            {
                update_progress(i + 1, total_lines);
            }
        }
        // std::cout << "Thread " << thread_id << " finished" << std::endl;
        // std::cout.flush();

        // Merge thread-local results into global results
        #pragma omp critical
        {
            smiles_basis.insert(local_smiles_basis.begin(), local_smiles_basis.end());
        }
    }

    std::cout << std::endl;

    return smiles_basis;
}
