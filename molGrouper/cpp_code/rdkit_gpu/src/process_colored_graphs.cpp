#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>
#include <numeric>
#include <exception>
#include <iterator>
#include <set>
#include <tuple>
#include <functional>
#include <unordered_map>
#include <sys/resource.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

void removeDuplicatesGPU(thrust::device_vector<unsigned long>& hash_vector);

void logMemoryUsage(const std::string& message) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << message << " - Memory usage: " << usage.ru_maxrss << " kilobytes" << std::endl;
}

template<typename T>
size_t getHostVectorMemoryUsage(const thrust::host_vector<T>& vec) {
    return vec.capacity() * sizeof(T);
}

template<typename T>
size_t getVectorMemoryUsage(const std::vector<T>& vec) {
    return vec.capacity() * sizeof(T);
}

// Helper function to estimate memory usage of unordered maps
template<typename K, typename V>
size_t getUnorderedMapMemoryUsage(const std::unordered_map<K, V>& umap) {
    return umap.size() * (sizeof(K) + sizeof(V)) + sizeof(umap);
}

// Helper function to estimate memory usage of unordered sets
template<typename T>
size_t getUnorderedSetMemoryUsage(const std::unordered_set<T>& uset) {
    return uset.size() * sizeof(T) + sizeof(uset);
}

class Logger {
public:
    template <typename T>
    void logSize(const std::string& name, const T& var) {
        sizes[name] = sizeof(var);
        types[name] = typeid(var).name();
    }

    void printSizes() const {
        for (const auto& entry : sizes) {
            std::cout << "Size of " << entry.first << " (" << types.at(entry.first) << "): " << entry.second << " bytes" << std::endl;
        }
    }

private:
    std::unordered_map<std::string, size_t> sizes;
    std::unordered_map<std::string, std::string> types;
};


void get_sensible_port_combos(
    const std::vector<std::vector<std::pair<int, int>>>& v, 
    const GroupGraph& gG, 
    const std::vector<std::pair<int, int>>& edge_list, 
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    const std::unordered_map<std::string, std::string>& nodeTypeToSmiles,
    const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex,
    thrust::host_vector<unsigned long>& host_hash_vector, 
    std::unordered_map<unsigned long, std::string>& hash_to_smiles,
    const bool verbose = false) {
        
    if (v.empty()) {
        return;
    }

    GroupGraph temp_gG(node_types);
    std::hash<std::string> hasher;
    unsigned long hash_value;
    std::string smiles;

    // Calculate the total size of the Cartesian product
    auto size = std::accumulate(v.begin(), v.end(), 1,
                                [](size_t s, const std::vector<std::pair<int, int>>& sub) { return s * sub.size(); });

    std::cout << "Size of Cartesian product: " << size << std::endl;

    // Initialize counters for each vector of pairs
    std::vector<size_t> counters(v.size(), 0);

    for (size_t i = 0; i < size; ++i) {
        std::vector<std::pair<int, int>> combination;

        for (size_t j = 0; j < v.size(); ++j) {
            combination.push_back(v[j][counters[j]]);
        }

        try {
            temp_gG = gG;
            for (size_t j = 0; j < edge_list.size(); ++j) {
                temp_gG.addEdge(
                    edge_list[j].first, combination[j].first,
                    edge_list[j].second, combination[j].second
                    );
            }

            // convert to molecular graph
            std::string smiles = temp_gG.toSmiles(nodeTypeToSmiles, nodeTypePortToIndex);
            // smiles = "CO"; // For testing purposes

            // Convert to hash
            hash_value = static_cast<unsigned long>(hasher(smiles));
            host_hash_vector.push_back(hash_value);
        } 
        catch (const std::exception& e) {
            if (verbose) {
                std::cerr << "Couldn't produce graph from vcolg output: " << e.what() << std::endl;
            }
        }

        // Update counters for the next combination
        for (size_t j = v.size(); j-- > 0;) {
            if (++counters[j] < v[j].size()) {
                break;
            }
            counters[j] = 0;
        }
    }
}


std::pair<std::unordered_map<unsigned long, std::string>, thrust::host_vector<unsigned long>> process_nauty_graph_vcolg_output(
        const std::string& line, 
        const std::unordered_map<std::string, std::vector<int>>& node_types,
        const std::unordered_map<int, std::string>& int_to_node_type,
        const std::unordered_map<std::string, std::string>& nodeTypeToSmiles,
        const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex,
        bool verbose = false) {

    // logMemoryUsage("Before processing");
    std::cout<<"Before processing"<<std::endl;

    std::vector<GroupGraph> group_graphs_list;
    std::unordered_set<std::string> graph_basis;

    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        return {};
    }

    std::string node_description_str = line.substr(0, split_pos);
    std::string edge_description_str = line.substr(split_pos + 2);

    std::istringstream node_desc_iss(node_description_str);
    std::vector<std::string> node_description(std::istream_iterator<std::string>{node_desc_iss},
                                              std::istream_iterator<std::string>());

    std::istringstream edge_desc_iss(edge_description_str);
    std::vector<std::string> edge_description(std::istream_iterator<std::string>{edge_desc_iss},
                                              std::istream_iterator<std::string>());

    // logMemoryUsage("After splitting input");

    std::vector<std::pair<int, int>> edge_list;
    for (size_t i = 0; i < edge_description.size(); i += 2) {
        edge_list.emplace_back(std::stoi(edge_description[i]), std::stoi(edge_description[i + 1]));
    }

    int n_vertices = std::stoi(node_description[0]);
    int n_edges = std::stoi(node_description[1]);
    std::vector<int> colors;
    for (size_t i = 2; i < node_description.size(); ++i) {
        colors.push_back(std::stoi(node_description[i]));
    }

    std::vector<std::pair<int, int>> non_colored_edge_list;

    GroupGraph gG(node_types);
    for (int i = 0; i < n_vertices; ++i) {
        gG.addNode(i, int_to_node_type.at(colors[i]));
    }

    std::map<int, std::vector<std::pair<int, int>>> edge_index_to_edge_color;
    for (size_t i = 0; i < edge_list.size(); ++i) {
        int src = edge_list[i].first;
        int dst = edge_list[i].second;
        std::vector<int> src_ports = node_types.at(int_to_node_type.at(colors[src]));
        std::vector<int> dst_ports = node_types.at(int_to_node_type.at(colors[dst]));
        non_colored_edge_list.push_back({src, dst});
        std::vector<std::pair<int, int>> port_combinations;
        for (int sp : src_ports) {
            for (int dp : dst_ports) {
                port_combinations.push_back({sp, dp});
            }
        }
        edge_index_to_edge_color[i] = port_combinations;
    }

    // logMemoryUsage("After processing edges");

    // Use a graph library to create a graph from non_colored_edge_list and check degrees
    std::map<int, int> node_degrees; // Placeholder for node degrees
    for (const auto& edge : non_colored_edge_list) {
        node_degrees[edge.first]++;
        node_degrees[edge.second]++;
    }
    for (const auto& [node, degree] : node_degrees) {
        if (degree > node_types.at(int_to_node_type.at(colors[node])).size()) {
            return {};
        }
    }

    std::vector<std::vector<std::pair<int, int>>> edge_colors;
    for (const auto& [index, port_combinations] : edge_index_to_edge_color) {
        edge_colors.push_back(port_combinations);
    }

    // logMemoryUsage("After checking degrees");
    thrust::host_vector<unsigned long> host_hash_vector; //unsigned long because hash is 64 bit on our system
    std::unordered_map<unsigned long, std::string> hash_to_smiles;


    // Generate and filter port combinations
    get_sensible_port_combos(
        edge_colors, 
        gG, 
        non_colored_edge_list, 
        node_types, 
        nodeTypeToSmiles,
        nodeTypePortToIndex,
        host_hash_vector, 
        hash_to_smiles,
        verbose
    );

    logMemoryUsage("After processing port combinations");


    thrust::device_vector<unsigned long> device_hash_vector = host_hash_vector;
    removeDuplicatesGPU(device_hash_vector);
    cudaDeviceSynchronize();
    host_hash_vector.resize(device_hash_vector.size());
    thrust::copy(device_hash_vector.begin(), device_hash_vector.end(), host_hash_vector.begin());




    return {hash_to_smiles, host_hash_vector};
}



// int main() {

//     std::unordered_map<std::string, std::vector<int>> node_types = {
//         {"N", {0}},
//         {"CO", {0, 1}},
//         {"CC", {0, 1, 2, 3}}
//     };
//     std::unordered_map<int, std::string> int_to_node_type = {
//         {0, "N"},
//         {1, "CO"},
//         {2, "CC"},
//         {3, "CC"},
//         {4, "CC"}
//     };
//     std::ifstream input_file("/Users/kieran/projects/molGrouper/molGrouper/cpp_code/vcolg_out.txt");
//     if (!input_file.is_open()) {
//         std::cerr << "Error opening input file." << std::endl;
//         return 1;
//     }
//     std::string line;
//     while (std::getline(input_file, line)) {
//         if (line.empty()) {
//             continue; // Skip empty lines
//         }
//         auto result = process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, true);
//         std::cout << "Processed line. Result size: " << result.size() << std::endl;
//     }
//     input_file.close();

//     return 0;
// }
