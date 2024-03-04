#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>

/**
 * @brief Represents a graph with nodes, ports, and edges.
 */

class GroupGraph {
private:
    using NodeID = int; // Change this type according to your needs
    using PortType = int; // Change this type according to your needs

    /**
     * @brief Represents a node in the graph.
     */
    struct Node {
        std::string type;
        std::vector<PortType> ports;
    };

    std::unordered_map<NodeID, Node> nodes; ///< Map of node IDs to their respective nodes.
    std::vector<std::tuple<NodeID, PortType, NodeID, PortType>> edges; ///< List of edges connecting nodes.
    std::unordered_map< std::string, std::vector<PortType> > nodeTypes; ///< Map of node types to their respective ports.

public:
    /**
     * @brief Constructs a GroupGraph with specified node types.
     * @param types A map of node types to their respective ports.
     */
    GroupGraph(const std::unordered_map<std::string, std::vector<PortType>>& types) {
        // Constructor to initialize the node types
        nodeTypes = types;
    }

    /**
     * @brief Adds a node to the graph.
     * @param id The identifier for the node.
     * @param type The type of the node.
     * @throws std::invalid_argument if the node already exists.
     */
    void addNode(NodeID id, const std::string& type) {
        // Add a node to the graph
        if (nodes.find(id) != nodes.end()) {
            throw std::invalid_argument("Node already exists");
        }
        nodes[id] = {type, nodeTypes[type]};
    }

    /**
     * @brief Adds an edge between two nodes in the graph.
     * @param from The ID of the first node.
     * @param fromPort The port on the first node.
     * @param to The ID of the second node.
     * @param toPort The port on the second node.
     * @throws std::invalid_argument if the nodes or ports do not exist.
     */
    void addEdge(NodeID from, PortType fromPort, NodeID to, PortType toPort) {
        // Add an edge to the graph
        if (nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) {
            throw std::invalid_argument("Node does not exist");
        }
        if (std::find(nodes[from].ports.begin(), nodes[from].ports.end(), fromPort) == nodes[from].ports.end()) {
            throw std::invalid_argument("Port does not exist");
        }
        edges.push_back(std::make_tuple(from, fromPort, to, toPort));
    }
    
    /**
     * @brief Gets the number of free ports on a node.
     * @param nodeID The identifier for the node.
     * @return The number of free ports.
     * @throws std::invalid_argument if the node does not exist.
     */
    int n_free_ports(NodeID nodeID) const {
        // Get the number of free ports on a node
        if (nodes.find(nodeID) == nodes.end()) {
            throw std::invalid_argument("Node does not exist");
        }

        const Node& node = nodes.at(nodeID);
        int occupied_ports = 0;

        for (const auto& edge : edges) {
            if (std::get<0>(edge) == nodeID) {
                occupied_ports++;
            }
        }

        return node.ports.size() - occupied_ports;
    }

    /**
     * @brief Converts the GroupGraph to an undirected graph.
     */
    void make_undirected() {
        std::vector<std::tuple<NodeID, PortType, NodeID, PortType>> undirectedEdges;
        for (const auto& edge : edges) {
            NodeID node1 = std::get<0>(edge);
            PortType port1 = std::get<1>(edge);
            NodeID node2 = std::get<2>(edge);
            PortType port2 = std::get<3>(edge);

            // Add an undirected edge
            undirectedEdges.push_back(std::make_tuple(node2, port2, node1, port1));
        }
        // Add the undirected edges to the main edges container
        edges.insert(edges.begin(), undirectedEdges.begin(), undirectedEdges.end());
    }

    /**
     * @brief Prints the nodes and edges of the GroupGraph.
     */
    void printGraph() const {
        // Print the nodes and edges of the GroupGraph
        std::cout << "Nodes:\n";
        for (const auto& entry : nodes) {
            std::cout << "Node " << entry.first << " ("<<entry.second.type<<") " << ": Ports ";
            for (PortType port : entry.second.ports) {
                std::cout << port << " ";
            }
            std::cout << "\n";
        }

        std::cout << "\nEdges:\n";
        for (const auto& edge : edges) {
            std::cout << "Edge: " 
            << std::get<0>(edge) << "(" << std::get<1>(edge)  << ")"
            << " -> "
            << std::get<2>(edge) << "(" << std::get<3>(edge)<< ")"<<"\n";
        }
    }
};

int main() {
    // Example usage
    GroupGraph myGraph({{"type1", {1, 2, 3}}, {"type2", {4, 5, 6}}});

    myGraph.addNode(1, "type1");
    myGraph.addNode(2, "type2");
    myGraph.addEdge(1,2, 2,5);
    myGraph.addEdge(1,3, 2,4);

    myGraph.make_undirected();

    myGraph.printGraph();

    std::cout << "Free Ports for Node 1: " << myGraph.n_free_ports(1) << "\n";
    std::cout << "Free Ports for Node 2: " << myGraph.n_free_ports(2) << "\n\n";


    return 0;
}
