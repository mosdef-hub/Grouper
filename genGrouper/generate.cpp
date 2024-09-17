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

#include <libpq-fe.h>
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

// Function to execute an SQL query with parameters
void executeQuery(PGconn* conn, const std::string& query, const std::vector<std::string>& params) {
    std::vector<const char*> param_values(params.size());
    for (size_t i = 0; i < params.size(); ++i) {
        param_values[i] = params[i].c_str();
    }

    // Execute the query with parameters
    PGresult* res = PQexecParams(
        conn,                // Connection
        query.c_str(),       // SQL query
        params.size(),       // Number of parameters
        NULL,                // Parameter types (NULL means let libpq infer types)
        param_values.data(), // Parameter values
        NULL,                // Parameter lengths (NULL means use default lengths)
        NULL,                // Parameter formats (NULL means use text format)
        0                    // Result format (0 for text, 1 for binary)
    );

    // Check the result status
    if (PQresultStatus(res) != PGRES_COMMAND_OK) {
        std::cerr << "Query failed: " << PQerrorMessage(conn) << std::endl;
    }
    PQclear(res);
}

std::unordered_set<GroupGraph> exhaustiveGenerate(
    int n_nodes, 
    std::unordered_set<GroupGraph::Node> node_defs, 
    std::string nauty_path,
    std::string input_file_path = "",
    int num_procs = -1,
    std::unordered_map<std::string, int> positiveConstraints = {},
    std::unordered_set<std::string> negativeConstraints = {},
    bool writeToDatabase = false,
    bool verbose = false
) {
    // Error handling

    if (n_nodes < 1) {
        throw std::invalid_argument("Number of nodes must be greater than 0...");
    }
    if (node_defs.size() < 1) {
        throw std::invalid_argument("Node definitions must not be empty...");
    }
    if (num_procs <= -1){
        num_procs = omp_get_max_threads();
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
    std::unordered_set<std::string> canon_basis;

    omp_set_num_threads(num_procs);      // Set the number of threads to match

    std::cout<< "Using "<<num_procs << " processors" << std::endl;

    #pragma omp parallel
    {

        // int thread_id = omp_get_thread_num();
        // std::cout<< "Thread " << thread_id << " started" << std::endl;
        std::cout.flush();

        // Thread-local data structures
        std::unordered_set<GroupGraph> local_basis;

        #pragma omp for schedule(dynamic) nowait
        for (int i = 0; i < total_lines; ++i) {

            process_nauty_output(
                lines[i], 
                node_defs,
                &local_basis,
                positiveConstraints, 
                negativeConstraints, 
                verbose
            );

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
            for (const auto& graph : local_basis) {
                if (canon_basis.find(graph.toSmiles()) == canon_basis.end()) {
                    canon_basis.insert(graph.toSmiles());
                    global_basis.insert(graph);
                }
            }
        }
    }

    std::cout << std::endl;

    std::cout<< "Generated " << global_basis.size()<< " unique graphs" << std::endl;

    if (writeToDatabase) {
        // Database connection parameters
        // dbname, user, password, hostaddr, port, 
        // Connection parameters
        const char* conninfo = "dbname=foo user=tmpuser password=password hostaddr=127.0.0.1 port=5432";
        
        // Establish a connection to the PostgreSQL database
        PGconn* conn = PQconnectdb(conninfo);
        if (PQstatus(conn) != CONNECTION_OK) {
            std::cerr << "Connection to database failed: " << PQerrorMessage(conn) << std::endl;
            PQfinish(conn);
        }
        try {
            // Prepare SQL statement for insertion
            std::string insert_query = "INSERT INTO graph_table (smiles, graph_data) VALUES ($1, $2)";

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
                // Convert the GroupGraph to a SMILES string
                std::string smiles = graph.toSmiles();

                // Serialize or otherwise convert GroupGraph to a format for storage (e.g., JSON, binary)
                std::string graph_data = graph.serialize();  // Assuming you have a serialize method

                // Execute the prepared statement
                std::vector<std::string> params = { smiles, graph_data };
                executeQuery(conn, insert_query, params);
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

        // Close the connection
        PQfinish(conn);
        std::cout << "Connection closed." << std::endl;
    }

    return global_basis;
}
