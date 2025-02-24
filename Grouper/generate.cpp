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

#include <libpq-fe.h>
#include <omp.h>
#include <stdio.h>

#include <nauty/nauty.h>

#include "dataStructures.hpp"
#include "processColoredGraphs.hpp"

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

std::unordered_set<GroupGraph> exhaustiveGenerate(
    int n_nodes, 
    std::unordered_set<GroupGraph::Group> node_defs, 
    std::string nauty_path,
    std::string input_file_path = "",
    int num_procs = -1,
    std::unordered_map<std::string, int> positiveConstraints = {},
    std::unordered_set<std::string> negativeConstraints = {},
    std::string config_path = "",
    bool verbose = false
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
    if (nauty_path.empty()) {
        throw std::invalid_argument("Nauty path must not be empty...");
    }
    if (!positiveConstraints.empty()){
        for (const auto& constraint : positiveConstraints) {
            if (constraint.second < 0) {
                throw std::invalid_argument("Positive constraint value must be greater than or equal to 0...");
            }
        }
    }

    if (verbose) {
        std::cout << "Number of nodes: " << n_nodes << std::endl;
        std::cout << "Number of node definitions: " << node_defs.size() << std::endl;
        std::cout << "Input file path: " << input_file_path << std::endl;
        std::cout << "Number of processors: " << num_procs << std::endl;
    }

    if (input_file_path.empty()) {
        // Call nauty
        std::string geng_command = nauty_path + "/geng " + std::to_string(n_nodes) + " -ctf > geng_out.txt";
        std::string vcolg_command = nauty_path + "/vcolg geng_out.txt -T -m" + std::to_string(node_defs.size()) + " > vcolg_out.txt";

        if (verbose) {
            std::cout << "Calling geng..." << std::endl;
        }
        system(geng_command.c_str());
        system(vcolg_command.c_str());
    }

    input_file_path = input_file_path.empty() ? "vcolg_out.txt" : input_file_path;

    // Read the input file
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

    if (total_lines == 0) {
        throw std::runtime_error("No lines found in input file...");
    }

    std::unordered_set<GroupGraph> global_basis;
    std::unordered_set<std::vector<setword>, hash_vector> canon_basis;

    omp_set_num_threads(num_procs);      // Set the number of threads to match

    std::cout<< "Using "<<num_procs << " processors" << std::endl;

    
    #pragma omp parallel
    {
        // Thread-local nauty structures using std::vector
        int n = 20; // Max number of nodes (adjustable)
        int m = SETWORDSNEEDED(n);
        std::vector<setword> g(m * n, 0);
        std::vector<int> lab(n, 0), ptn(n, 0), orbits(n, 0);
        DEFAULTOPTIONS_GRAPH(options);
        statsblk stats;

        // Thread-local basis set
        std::unordered_set<GroupGraph> local_basis;

        #pragma omp for schedule(dynamic) nowait
        for (int i = 0; i < total_lines; ++i) {
            process_nauty_output(
                lines[i], 
                node_defs,
                &local_basis,
                positiveConstraints, 
                negativeConstraints, 
                verbose,
                g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats
            );

            #pragma omp critical
            {
                update_progress(i + 1, total_lines);
            }
        }

        #pragma omp critical
        {
            for (const auto& graph : local_basis) {
                // if (canon_basis.find(graph.toSmiles()) == canon_basis.end()) {
                //     canon_basis.insert(graph.toSmiles());
                //     global_basis.insert(graph);
                // }
                if (canon_basis.find(graph.canonize()) == canon_basis.end()) {
                    canon_basis.insert(graph.canonize());
                    global_basis.insert(graph);
                }
            }
        }
    }

    std::cout << std::endl;


    if (!config_path.empty()) {

        std::unordered_map<std::string, std::string> configParams = parseConfig(config_path);

        std::string table_name = configParams["table_name"]; // Get the table name from config

        std::string conninfo = "dbname=" + configParams["dbname"] +
                       " user=" + configParams["user"] +
                       " password=" + configParams["password"] +
                       " hostaddr=" + configParams["hostaddr"] +
                       " port=" + configParams["port"];

        PGconn* conn = PQconnectdb(conninfo.c_str());
        if (PQstatus(conn) != CONNECTION_OK) {
                std::cerr << "Connection to database failed: " << PQerrorMessage(conn) << std::endl;
                PQfinish(conn);
            }

        try {
                // Sanity check: Check if connection is in a valid state
                if (PQstatus(conn) != CONNECTION_OK) {
                    std::cerr << "Database connection is not in a valid state: " << PQerrorMessage(conn) << std::endl;
                    PQfinish(conn);
                }

                // Sanity check: Check if the table exists and create it if it doesn't
                std::string create_table_query = 
                    "CREATE TABLE IF NOT EXISTS " + table_name + " ("
                    "smiles TEXT PRIMARY KEY, "
                    "graph_data TEXT, "   
                    "n_nodes INT"         
                    ");";

                executeQuery(conn, create_table_query);

                // Prepare SQL statement for insertion
                std::string insert_query = "INSERT INTO " + table_name +" (smiles, graph_data, n_nodes) VALUES ($1, $2, $3)";

                // Start a transaction block
                PGresult* res = PQexec(conn, "BEGIN");
                if (PQresultStatus(res) != PGRES_COMMAND_OK) {
                    std::cerr << "BEGIN command failed: " << PQerrorMessage(conn) << std::endl;
                    PQclear(res);
                    PQfinish(conn);
                }
                PQclear(res);

                // Loop through the global_basis set and insert each GroupGraph into the database
                for (const auto& graph : global_basis) {
                    std::string smiles = graph.toSmiles();
                    std::string graph_data = graph.serialize();
                    int n_nodes = graph.nodes.size();

                    // Prepare parameters for the insert query
                    const char* paramValues[3];
                    paramValues[0] = smiles.c_str();
                    paramValues[1] = graph_data.c_str();
                    paramValues[2] = std::to_string(n_nodes).c_str();

                    // Execute the insert query
                    PGresult* insertRes = PQexecParams(conn, insert_query.c_str(), 3, nullptr, paramValues, nullptr, nullptr, 0);
                    if (PQresultStatus(insertRes) != PGRES_COMMAND_OK) {
                        std::cerr << "Insert command failed for SMILES: " << smiles << ", Error: " << PQerrorMessage(conn) << std::endl;
                        PQclear(insertRes);
                        // Optionally, you could choose to rollback the transaction here if any insert fails
                        continue;  // Skip to the next graph if insertion fails
                    }
                    PQclear(insertRes);
                }

                // Commit the transaction
                res = PQexec(conn, "COMMIT");
                if (PQresultStatus(res) != PGRES_COMMAND_OK) {
                    std::cerr << "COMMIT command failed: " << PQerrorMessage(conn) << std::endl;
                    PQclear(res);
                    PQfinish(conn);
                }
                PQclear(res);

        } catch (const std::exception& e) {
                std::cerr << "Database error: " << e.what() << std::endl;
                PQfinish(conn);
            }

            PQfinish(conn);
            std::cout << "Connection closed." << std::endl;
    }
    
    std::cout<< "Number of unique graphs: " << global_basis.size() << std::endl;

    return global_basis;
}
