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
#include <filesystem>
#include <random>
#include <algorithm>

#include <libpq-fe.h>
#include <omp.h>
#include <stdio.h>

#include <nauty/nauty.h>

#include "dataStructures.hpp"
#include "processColoredGraphs.hpp"
#include "generate.hpp"

#define MAX_EDGES 100

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

struct hash_vector {
    std::size_t operator()(const std::vector<setword>& v) const {
        std::size_t seed = 0;
        for (const auto& i : v) {
            boost::hash_combine(seed, i);
        }
        return seed;
    }
};

std::unordered_map<std::string, std::string> parseConfig(const std::string& configFile) {
    std::unordered_map<std::string, std::string> configParams;
    std::ifstream config(configFile);
    std::string line;
    while (std::getline(config, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            configParams[key] = value;
        }
    }
    return configParams;
}

void executeQuery(PGconn* conn, const std::string& query, const std::vector<std::string>& params = {}) {
    const char* paramValues[params.size()];
    for (size_t i = 0; i < params.size(); ++i) {
        paramValues[i] = params[i].c_str();
    }
    PGresult* res = PQexecParams(conn, query.c_str(), params.size(), nullptr, paramValues, nullptr, nullptr, 0);
    if (PQresultStatus(res) != PGRES_COMMAND_OK) {
        std::cerr << "Query execution failed: " << PQerrorMessage(conn) << std::endl;
    }
    PQclear(res);
}

std::tuple< int, std::vector<int>, std::vector<std::pair<int, int>> > parse_nauty_graph_line_tmp(const std::string& line, const std::unordered_set<GroupGraph::Group>& node_defs) {
    // Split the line into node_description and edge_description
    size_t split_pos = line.find("  ");
    if (split_pos == std::string::npos) {
        throw std::runtime_error("Invalid nauty output line...");
    }

    std::string node_description_str = line.substr(0, split_pos);
    std::string edge_description_str = line.substr(split_pos + 2);

    std::istringstream node_desc_iss(node_description_str);
    std::vector<std::string> node_description(std::istream_iterator<std::string>{node_desc_iss},
                                              std::istream_iterator<std::string>());

    std::istringstream edge_desc_iss(edge_description_str);
    std::vector<std::string> edge_description(std::istream_iterator<std::string>{edge_desc_iss},
                                              std::istream_iterator<std::string>());

    std::vector<std::pair<int, int>> edge_list;
    for (size_t i = 0; i < edge_description.size(); i += 2) {
        edge_list.emplace_back(std::stoi(edge_description[i]), std::stoi(edge_description[i + 1]));
    }



    int n_vertices = std::stoi(node_description[0]);
    // int n_edges = std::stoi(node_description[1]);
    std::vector<int> colors;
    for (size_t i = 2; i < node_description.size(); ++i) {
        colors.push_back(std::stoi(node_description[i]));
    }

    if (node_defs.size() < static_cast<size_t>(*std::max_element(colors.begin(), colors.end()) + 1)) {
        throw std::runtime_error("Number of nodes in node_defs does not match the number of nodes in the nauty_output_file...");
    }

    return {n_vertices, colors, edge_list};
}

bool check_max_bond_not_exceeded_tmp(
    const std::vector<std::pair<int,int>>& edge_list,
    const std::vector<int>& colors,
    const std::unordered_map<std::string, std::vector<int>>& node_types,
    const std::unordered_map<int, std::string>& int_to_node_type) {

    std::vector<int> node_bond_count;
    for (std::size_t i = 0; i < colors.size(); ++i) {
        std::string node_type = int_to_node_type.at(colors[i]);
        int n_ports = node_types.at(node_type).size();
        node_bond_count.push_back(n_ports);
    }
    for (const auto& edge : edge_list) {
        int src = edge.first;
        int dst = edge.second;
        node_bond_count[src] -= 1;
        node_bond_count[dst] -= 1;
    }
    for (const auto& count : node_bond_count) {
        if (count < 0) {
            return false;
        }
    }
    return true;
}




void insertGraph(PGconn* conn, const GroupGraph& graph, const std::string& table_name) {
    std::string smiles = graph.toSmiles();
    std::string graph_data = graph.serialize();
    int n_nodes = graph.nodes.size();

    const char* paramValues[3];
    paramValues[0] = smiles.c_str();
    paramValues[1] = graph_data.c_str();
    paramValues[2] = std::to_string(n_nodes).c_str();

    std::string insert_query = "INSERT INTO " + table_name + " (smiles, graph_data, n_nodes) VALUES ($1, $2, $3) ON CONFLICT (smiles) DO NOTHING";
    PGresult* insertRes = PQexecParams(conn, insert_query.c_str(), 3, nullptr, paramValues, nullptr, nullptr, 0);

    if (PQresultStatus(insertRes) != PGRES_COMMAND_OK && PQresultStatus(insertRes) != PGRES_TUPLES_OK) {
    }
    PQclear(insertRes);
}

thread_local size_t isomorphism_checks = 0;

std::pair<std::unordered_set<GroupGraph>, size_t> exhaustiveGenerate(
    int n_nodes,
    std::unordered_set<GroupGraph::Group> node_defs,
    int num_procs = -1,
    std::string vcolg_output_file = "",
    std::unordered_map<std::string, int> positiveConstraints = {},
    std::unordered_set<std::string> negativeConstraints = {},
    std::string config_path = ""
) {


    // Error handling
    if (n_nodes < 1) {
        throw std::invalid_argument("Number of nodes must be greater than 0...");
    }
    if (node_defs.size() < 1) {
        throw std::invalid_argument("Group definitions must not be empty...");
    }
    if (num_procs <= -1){
        num_procs = omp_get_max_threads();
    }
    if (!positiveConstraints.empty()){
        for (const auto& constraint : positiveConstraints) {
            if (constraint.second < 0) {
                throw std::invalid_argument("Positive constraint value must be greater than or equal to 0...");
            }
        }
    }


    if (vcolg_output_file.empty()) {
        // Call nauty
        std::string geng_command = "geng " + std::to_string(n_nodes) + " -c > geng_out.txt";
        std::string vcolg_command = "vcolg geng_out.txt -T -m" + std::to_string(node_defs.size()) + " > vcolg_out.txt";
        system(geng_command.c_str());
        system(vcolg_command.c_str());
    }

    vcolg_output_file = vcolg_output_file.empty() ? "vcolg_out.txt" : vcolg_output_file;

    // Read the input file
    std::ifstream input_file(vcolg_output_file);
    if (!input_file.is_open()) {
        throw std::runtime_error("Error opening vcolg output file, check vcolg_output_file argument and make sure nauty is installed...");
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

    std::cout << "Processing " << total_lines << " lines from vcolg output..." << std::endl;

    input_file.close(); // Close the input file as it's no longer needed

    if (total_lines == 0) {
        throw std::runtime_error("No lines found in vcolg output file...");
    }

    std::unordered_set<GroupGraph> global_basis;
    omp_set_num_threads(num_procs);      // Set the number of threads to match
    int n_finished = 0;
    std::cout<< "Using "<<num_procs << " processors" << std::endl;

    if (!config_path.empty()) {
        std::unordered_map<std::string, std::string> configParams = parseConfig(config_path);
        std::string table_name = configParams["table_name"];
        std::string conninfo = "dbname=" + configParams["dbname"] +
                       " user=" + configParams["user"] +
                       " password=" + configParams["password"] +
                       " hostaddr=" + configParams["hostaddr"] +
                       " port=" + configParams["port"];

        PGconn* main_conn = PQconnectdb(conninfo.c_str());
        if (PQstatus(main_conn) != CONNECTION_OK) {
            std::cerr << "Connection to database failed: " << PQerrorMessage(main_conn) << std::endl;
            PQfinish(main_conn);
            throw std::runtime_error("Database connection failed.");
        }
        std::string create_table_query =
            "CREATE TABLE IF NOT EXISTS " + table_name + " (" +
            "smiles TEXT PRIMARY KEY, " +
            "graph_data TEXT, " +
            "n_nodes INT" +
            ");";
        executeQuery(main_conn, create_table_query);
        PQfinish(main_conn);

        #pragma omp parallel
        {
            PGconn* conn = PQconnectdb(conninfo.c_str());
            if (PQstatus(conn) != CONNECTION_OK) {
                #pragma omp critical
                {
                    std::cerr << "Thread " << omp_get_thread_num() << " failed to connect to DB: " << PQerrorMessage(conn) << std::endl;
                }
            } else {
                std::unordered_set<GroupGraph> local_basis;
                int n = 20; // Max number of nodes (adjustable)
                int m = SETWORDSNEEDED(n);
                std::vector<setword> g(m * n, 0);
                std::vector<int> lab(n, 0), ptn(n, 0), orbits(n, 0);
                DEFAULTOPTIONS_GRAPH(options);
                statsblk stats;

                #pragma omp for schedule(dynamic) nowait
                for (int i = 0; i < total_lines; ++i) {
                    local_basis.clear();
                    process_nauty_output(
                        lines[i], node_defs, &local_basis,
                        positiveConstraints, negativeConstraints,
                        g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats
                    );

                    for (const auto& graph : local_basis) {
                        insertGraph(conn, graph, table_name);
                    }

                    if ((++n_finished % 1) == 0 || n_finished == total_lines) {
                        #pragma omp critical
                        {
                            update_progress(n_finished, total_lines);
                        }
                    }
                }
                PQfinish(conn);
            }
        }
        std::cout << std::endl;
        std::cout << "Graphs saved to database." << std::endl;
        return {};
    } else {
        size_t total_isomorphism_checks = 0;
        size_t previous_n_colorings = 0;
        std::unordered_set<std::string> canon_basis;
        std::unordered_map<std::string, std::unordered_set<std::string>> vcolg_line_to_unique_smiles;
        int max_threads = num_procs;
        std::vector<std::unordered_set<GroupGraph>> all_local_bases(max_threads);

    #pragma omp parallel
    {
            isomorphism_checks = 0; // Reset for each thread
            int tid = omp_get_thread_num();
            std::unordered_set<GroupGraph>& local_basis = all_local_bases[tid];
            int n = 20; // Max number of nodes (adjustable)
            int m = SETWORDSNEEDED(n);
            std::vector<setword> g(m * n, 0);
            std::vector<int> lab(n, 0), ptn(n, 0), orbits(n, 0);
            DEFAULTOPTIONS_GRAPH(options);
            statsblk stats;

            #pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < total_lines; ++i) {
                process_nauty_output(
                    lines[i], node_defs, &local_basis,
                    positiveConstraints, negativeConstraints,
                    g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats
                );

                if ((++n_finished % 1) == 0 || n_finished == total_lines) {
                    #pragma omp critical
                    {
                        update_progress(n_finished, total_lines);
                    }
                }
            }
            #pragma omp critical
            {
                total_isomorphism_checks += isomorphism_checks;
            }
        }

        for (const auto& local_basis : all_local_bases) {
            for (const auto& graph : local_basis) {
                std::string smiles = graph.toSmiles();
                if (canon_basis.find(smiles) == canon_basis.end()) {
                    canon_basis.insert(smiles);
                    global_basis.insert(graph);
                }
            }
        }
        std::cout << std::endl;
        std::cout<< "Number of unique graphs: " << global_basis.size() << std::endl;

        // This part is not parallelized, but it should be fast enough
        for (int i = 0; i < total_lines; ++i) {
            std::unordered_set<GroupGraph> local_basis;
            int n = 20; // Max number of nodes (adjustable)
            int m = SETWORDSNEEDED(n);
            std::vector<setword> g(m * n, 0);
            std::vector<int> lab(n, 0), ptn(n, 0), orbits(n, 0);
            DEFAULTOPTIONS_GRAPH(options);
            statsblk stats;
            process_nauty_output(
                lines[i], node_defs, &local_basis,
                positiveConstraints, negativeConstraints,
                g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats
            );
            for (const auto& graph : local_basis) {
                vcolg_line_to_unique_smiles[lines[i]].insert(graph.toSmiles());
            }
        }

        std::cout << "\nUnique graphs per vcolg line:" << std::endl;
        for (const auto& pair : vcolg_line_to_unique_smiles) {
            std::cout << "  Line: \"" << pair.first << "\", Unique Graphs: " << pair.second.size();
            for (const auto& smi : pair.second) {
                std::cout << " " << smi;
            }
            std::cout << std::endl;
        }

        return {global_basis, total_isomorphism_checks};
    }
}

// This is the randomGenerate that utilizes the nauty library
std::pair<std::unordered_set<GroupGraph>, size_t> randomGenerate(
    int n_nodes,
    const std::unordered_set<GroupGraph::Group>& node_defs,
    int num_graphs = 100,
    int num_procs = -1,
    const std::unordered_map<std::string, int>& positiveConstraints = {},
    const std::unordered_set<std::string>& negativeConstraints = {}
) {
    if (n_nodes < 1) {
        throw std::invalid_argument("Number of nodes must be greater than 0...");
    }
    if (node_defs.empty()) {
        throw std::invalid_argument("Group definitions must not be empty...");
    }
    if (num_graphs < 1) {
        throw std::invalid_argument("Number of graphs must be greater than 0...");
    }
    for (const auto& constraint : positiveConstraints) {
        if (constraint.second < 0) {
            throw std::invalid_argument("Positive constraint value must be >= 0...");
        }
    }
    if (num_procs <= -1) {
        num_procs = omp_get_max_threads();
    }

    std::string geng_command = "geng " + std::to_string(n_nodes) + " -c > geng_out.txt";
    std::string vcolg_command = "vcolg geng_out.txt -T -m" + std::to_string(node_defs.size()) + " > vcolg_out.txt";
    system(geng_command.c_str());
    system(vcolg_command.c_str());

    std::ifstream input_file("vcolg_out.txt");
    if (!input_file.is_open()) {
        throw std::runtime_error("Error opening input file...");
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input_file, line)) {
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    input_file.close();
    if (lines.empty()) {
        throw std::runtime_error("No valid graphs found...");
    }

    std::vector<std::tuple<int, std::vector<int>, std::vector<std::pair<int, int>>>> possible_node_colored_graphs;

    std::unordered_set<GroupGraph> global_basis;
    std::unordered_set<std::string> global_smiles_basis;
    // std::random_device rd;

    // Create necessary maps
    std::unordered_map<int, std::string> int_to_node_type;
    std::unordered_map<int, GroupGraph::Group> color_to_group;
    std::unordered_map<std::string, std::vector<int>> node_types;
    std::unordered_map<std::string, std::vector<int>> node_type_to_hub;
    std::unordered_map<int, std::string> int_to_pattern;
    std::unordered_map<std::string, std::string> type_to_pattern;
    for (const auto& node : node_defs) {
        node_types[node.ntype] = node.ports;
        type_to_pattern[node.ntype] = node.pattern;
    }
    for (const auto& node: node_defs) {
        for (const auto& h : node.hubs) {
            node_type_to_hub[node.ntype].push_back(h);
        }
    }
    int i = 0;
    for (const auto& [node_type, ports] : node_types) {
        int_to_node_type[i] = node_type;
        int_to_pattern[i] = type_to_pattern[node_type];
        for (auto& g : node_defs) {
            if (g.ntype == node_type) {
                color_to_group[i] = g;
                break;  // Exit the loop once a match is found
            }
        }

        i++;
    }

    // Process lines and filter out incompatible graphs
    std::vector<std::tuple<int, std::vector<int>, std::vector<std::pair<int, int>>>> valid_node_colored_graphs;
    for (int line_idx = 0; line_idx < lines.size(); ++line_idx) {
        const auto& l = lines[line_idx];
        auto [n_vertices, colors, edge_list] = parse_nauty_graph_line_tmp(l, node_defs);
        if (!check_max_bond_not_exceeded_tmp(edge_list, colors, node_types, int_to_node_type)) continue;

        bool is_valid_graph = true;
        std::unordered_map<int, int> node_degrees;
        for (const auto& [src, dst] : edge_list) {
            node_degrees[src]++;
            node_degrees[dst]++;
        }

        for (int i = 0; i < n_vertices; ++i) {
            int color = colors[i];
            int degree = node_degrees[i];
            if (color_to_group.at(color).hubs.size() < degree) {
                is_valid_graph = false;
                break;
            }
        }

        if (is_valid_graph) {
            valid_node_colored_graphs.push_back(std::make_tuple(line_idx, colors, edge_list));
        }
    }
    possible_node_colored_graphs = valid_node_colored_graphs;

    if (possible_node_colored_graphs.empty()) {
        throw std::runtime_error("No valid node-colored graphs found after filtering based on hub counts.");
    }

    std::uniform_int_distribution<> dist_graph(0, int(possible_node_colored_graphs.size()) - 1);


    // Find maximum degree of any node in the node colored graphs
    std::unordered_map<int, int> max_degree; // color -> max_degree
    for (const auto& [line_idx, colors, edge_list] : possible_node_colored_graphs) {
        std::unordered_map<int, int> index_to_degree;
        for (const auto& [src, dst]: edge_list) {
            index_to_degree[src]++;
            index_to_degree[dst]++;
        }
        for (const auto& [index, degree] : index_to_degree) {
            if (max_degree.find(colors[index]) == max_degree.end()) {
                max_degree[colors[index]] = degree;
            } else {
                max_degree[colors[index]] = std::max(max_degree[colors[index]], degree);
            }
        }
    }

    // Step 1: Precompute all possible attachments for each group
    std::unordered_map<GroupGraph::Group, std::unordered_map<int, std::vector<std::vector<int>>>> group_attachments;
    for(const auto& [color, degree] : max_degree) {
        for( int i = 1; i <= degree; i++) {
            group_attachments[color_to_group[color]][i] = color_to_group[color].getPossibleAttachments(i);
        }
    }

    bool reachedMaxAttempts = true;
    size_t total_isomorphism_checks = 0;

    omp_set_num_threads(num_procs);
    #pragma omp parallel
    {
        isomorphism_checks = 0; // Reset for each thread
        std::unordered_set<GroupGraph> local_basis;
        // std::random_device rd;
        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution<> dist_graph(0, possible_node_colored_graphs.size() - 1);

        #pragma omp for schedule(dynamic)
        for (int attempt = 0; attempt < num_graphs * 1000; ++attempt) {
            if (global_basis.size() >= num_graphs) { // Check if another thread has already reached the target
                continue; // Continue to next iteration, but don't generate new graphs
            }
            GroupGraph candidate_graph;
            auto [vcolg_line_idx, colors, edge_list] = possible_node_colored_graphs[dist_graph(gen)];
            // printf("trying to generate graph from vcolg line %d: %s\n", vcolg_line_idx, lines[vcolg_line_idx].c_str());
            for (int color : colors) {
                GroupGraph::Group group = color_to_group[color];
                candidate_graph.addNode(group.ntype, group.pattern, group.hubs, group.patternType);
            }
            std::unordered_map<int, int> node_degrees;
            for (const auto& [src, dst] : edge_list) {
                node_degrees[src]++;
                node_degrees[dst]++;
            }
            std::vector<std::vector<int>> choosen_attachments;
            int c = 0;
            for (int color : colors) {
                int node_degree = node_degrees[c];
                std::vector<std::vector<int>> possible_attachments = group_attachments[color_to_group[color]][node_degree];
                if (possible_attachments.empty()) {
                    continue; // Skip this graph if no valid attachments for this degree
                }
                std::uniform_int_distribution<> dist_attachment(0, possible_attachments.size() - 1);
                choosen_attachments.push_back(possible_attachments[dist_attachment(gen)]);
                c++;
            }

            // Shuffle attachments for each node once to ensure random port allocation
            for (auto& attachments : choosen_attachments) {
                std::shuffle(attachments.begin(), attachments.end(), gen);
            }

            std::unordered_map<int, int> current_degree;
            for (int index = 0; index < n_nodes; index++) {
                current_degree[index] = 0;
            }
            for (const auto& [src, dst] : edge_list) {
                candidate_graph.addEdge(
                    {src, choosen_attachments[src][current_degree[src]]},
                    {dst, choosen_attachments[dst][current_degree[dst]]}
                );
                current_degree[src]++;
                current_degree[dst]++;
            }
            std::string smiles = candidate_graph.toSmiles();
            #pragma omp critical
            {
                if (global_smiles_basis.insert(smiles).second) { // If SMILES is unique
                    global_basis.insert(candidate_graph); // Add graph to the result set
                    update_progress(global_basis.size(), num_graphs);
                    if (int(global_basis.size()) >= num_graphs) {
                        reachedMaxAttempts = false; // Target reached, so max attempts not reached
                    }
                }
            }
        }
        #pragma omp critical
        {
            total_isomorphism_checks += isomorphism_checks;
        }
    }

    if (reachedMaxAttempts) {
        std::cout << "Warning: Maximum number of attempts reached without generating desired number of graphs..." << std::endl;
    }

    return {global_basis, total_isomorphism_checks};
}
