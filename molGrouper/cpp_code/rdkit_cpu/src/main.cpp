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

#include "GroupGraph.cpp"
#include "process_colored_graphs.cpp"

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

    std::ifstream input_file("/Users/kieran/foo/molGrouper/molGrouper/cpp_code/rdkit_cpu/vcolg_out.txt");
    if (!input_file.is_open()) {
        std::cerr << "Error opening input file." << std::endl;
        return 1;
    }

    std::cout
    << omp_get_num_procs() << " "
    << omp_get_max_threads() << " "
    << omp_get_thread_limit() << std::endl;

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

    int num_threads = 4; // Set this to your desired number of threads
    omp_set_num_threads(num_threads);

    #pragma omp parallel
    {
        // Thread-local data structures
        std::unordered_set<std::string> local_smiles_basis;

        #pragma omp for schedule(dynamic) nowait
        for (int i = 0; i < total_lines; ++i) {
            std::unordered_set<std::string> result = process_nauty_graph_vcolg_output(lines[i], node_types, int_to_node_type, nodeTypeToSmiles, nodeTypePortToIndex, false);

            for (auto it : result) {
                local_smiles_basis.insert(it);
            }

            #pragma omp critical
            {
                update_progress(i + 1, total_lines);
            }
        }

        // Merge thread-local results into global results
        #pragma omp critical
        {
            smiles_basis.insert(local_smiles_basis.begin(), local_smiles_basis.end());
        }
    }

    std::cout << std::endl;

    for (auto it : smiles_basis) {
        std::cout << it << std::endl;
    }

    return 0;
}
