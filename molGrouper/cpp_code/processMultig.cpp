#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
// #include "groupGraph.cpp"

/*
run geng, vcolg, multig

How can you account for symmetries in the group graph structure
*/



class MultigConverter {
public:
    static std::pair<int, int> multi_to_pair(int multi, int max_multi) {
        if (1 <= multi && multi <= max_multi) {
            int x = (multi - 1) % static_cast<int>(std::sqrt(max_multi)) + 1;
            int y = (multi - 1) / static_cast<int>(std::sqrt(max_multi)) + 1;

            x -= 1; // Convert to 0-indexed
            y -= 1; // Convert to 0-indexed
            return std::make_pair(x, y);
        } else {
            throw std::invalid_argument("Input multi must be in the range 1 to max_multi.");
        }
    }

    static GroupGraph multig_line_to_graph(
        const std::string& line,
        const std::unordered_map<int, std::string>& int_to_node_type,
        const std::unordered_map<std::string, std::vector<int>>& node_int_to_port,
        const std::unordered_map<std::string, std::vector<int>>& node_types,
        bool verbose = false
    ) {
        std::vector<std::string> edge_description;
        size_t pos = line.find("  ");
        if (pos != std::string::npos) {
            edge_description = split(line.substr(pos + 2), ' ');
        }

        std::vector<std::tuple<int, int, int>> edge_list;
        for (size_t i = 0; i < edge_description.size(); i += 3) {
            int source = std::stoi(edge_description[i]);
            int destination = std::stoi(edge_description[i + 1]);
            int edge_type = std::stoi(edge_description[i + 2]);
            edge_list.emplace_back(source, destination, edge_type);
        }

        std::vector<std::string> node_description = split(line.substr(0, pos), ' ');
        int n_vertices = std::stoi(node_description[0]);
        int n_edges = std::stoi(node_description[1]);
        std::vector<std::string> colors(node_description.begin() + 2, node_description.end());

        int max_n_attachments = 0;
        for (const auto& entry : node_types) {
            max_n_attachments = std::max(max_n_attachments, static_cast<int>(entry.second.size()));
        }

        GroupGraph groupG(node_types);
        const GroupGraph emptyG({{"placeholder1", {1,}}, {"placeholder2", {4,}}});
        // Add nodes
        std::string node_type;
        std::unordered_map<int, std::string>  index_to_node_type;
        for (int i = 0; i < n_vertices; ++i) {
            node_type = int_to_node_type.find(std::stoi(colors[i]))->second;
            index_to_node_type[i] = node_type;
            groupG.addNode(i, node_type);
        }

        // Add edges
        bool result;
        int port1, port2, source_node, dest_node, port_int_1, port_int_2;
        std::string source_node_type, dest_node_type;
        
        for (const auto e : edge_list) {
            port_int_1, port_int_2;
            std::tie(port_int_1, port_int_2) = multi_to_pair(std::get<2>(e), max_n_attachments * max_n_attachments);
            source_node = std::get<0>(e);
            dest_node = std::get<1>(e);
            source_node_type = index_to_node_type[source_node];
            dest_node_type = index_to_node_type[dest_node];
            port1 = node_int_to_port.find(source_node_type)->second[port_int_1];
            port2 = node_int_to_port.find(dest_node_type)->second[port_int_2];
            result = groupG.addEdge(std::get<0>(e), port1, std::get<1>(e), port2);
            if (!result) {
                break;
            }
        }
        if (result){
            return groupG;
        }
        else {
            return emptyG;
        }
        
    }

    static std::vector<GroupGraph> parse_multig_file() {
        // Specify the file path
        std::vector<GroupGraph> converted_graphs;
        std::string filePath = "multig_out.txt";

        std::unordered_map<std::string, std::vector<GroupGraph::PortType>> node_types;
        node_types["NH2"] = {0};
        node_types["CO"] = {0, 1};
        node_types["CC"] = {0, 1, 2, 3};

        std::unordered_map<int, std::string> int_to_node_type = {{0, "NH2"}, {1, "CO"}, {2, "CC"}};  // Maps node type integers to node type strings
        std::unordered_map<std::string, std::vector<int>> node_int_to_port = {{"NH2", {0}}, {"CO", {0, 1}}, {"CC", {0, 1, 2, 3}}};  // Maps node type strings to port integers


        // Open the file
        std::ifstream inputFile(filePath);

        // Check if the file is opened successfully
        if (!inputFile.is_open()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            return converted_graphs; 
        }

        // Read and process the contents of the file
        std::string line;
        GroupGraph output_graph(node_types);
        while (std::getline(inputFile, line)) {
            // Process each line as needed
            output_graph = multig_line_to_graph(
                line,
                int_to_node_type,
                node_int_to_port,
                node_int_to_port
            );
            converted_graphs.push_back(output_graph);
        }

        // Close the file
        inputFile.close();

        return converted_graphs;
    }
private:
    static std::vector<std::string> split(const std::string& s, char delimiter) {
        std::vector<std::string> tokens;
        size_t start = 0, end = 0;
        while ((end = s.find(delimiter, start)) != std::string::npos) {
            tokens.push_back(s.substr(start, end - start));
            start = end + 1;
        }
        tokens.push_back(s.substr(start));
        return tokens;
    }
};



// int main() {
//     auto output = MultigConverter::parse_multig_file();
    

//     return 0;
// }