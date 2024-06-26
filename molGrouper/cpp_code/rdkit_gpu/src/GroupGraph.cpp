#include <iostream>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>

/**
 * @brief Represents a graph with nodes, ports, and edges.
 */

class GroupGraph {
public:
    /**
     * @brief Constructs a GroupGraph with specified node types.
     * @param types A map of node types to their respective ports.
     */
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


    GroupGraph(const std::unordered_map<std::string, std::vector<PortType>>& types) {
        // Constructor to initialize the node types
        nodeTypes = types;
    }

    /**
     * @brief Copy constructor for GroupGraph.
     */
    GroupGraph(const GroupGraph& other) {
        // Copy nodes
        for (const auto& pair : other.nodes) {
            nodes[pair.first] = pair.second;
        }
        // Copy edges
        edges = other.edges;
        // Copy node types
        nodeTypes = other.nodeTypes;
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
    bool addEdge(NodeID from, PortType fromPort, NodeID to, PortType toPort, bool verbose = false ) { // return true if sucessful false otherwise
        if (n_free_ports(from) <= 0){
            if (verbose) {
                std::cout << "Source node doesn't have enough ports!"<<std::endl;
            }
            return false;
            // throw std::invalid_argument("Source node doesn't have enough ports!");
        }
        if (n_free_ports(to) <= 0){
            if (verbose) {
                std::cout << "Destination node doesn't have enough ports!"<<std::endl;
            }
            return false;
            // throw std::invalid_argument("Source node doesn't have enough ports!");
        }
        // Add an edge to the graph
        if (nodes.find(from) == nodes.end() || nodes.find(to) == nodes.end()) {
            throw std::invalid_argument("Node does not exist");
        }
        if (std::find(nodes[from].ports.begin(), nodes[from].ports.end(), fromPort) == nodes[from].ports.end()) {
            throw std::invalid_argument("Port does not exist");
        }
        edges.push_back(std::make_tuple(from, fromPort, to, toPort));
        // edges.push_back(std::make_tuple(to, toPort, from, fromPort)); // Uncomment this line if you want to add an edge in both directions
        return true;
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

    int numNodes() const {
        return nodes.size();
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

        std::cout << "Edges:\n";
        for (const auto& edge : edges) {
            std::cout << "Edge: " 
            << std::get<0>(edge) << "(" << std::get<1>(edge)  << ")"
            << " -> "
            << std::get<2>(edge) << "(" << std::get<3>(edge)<< ")"<<"\n";
        }
    }

    RDKit::ROMol* toMolecularGraph(
        const std::unordered_map<std::string, std::string>& nodeTypeToSmiles,
        const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex
    ) const {
        RDKit::RWMol* molecularGraph = new RDKit::RWMol();

        // Mapping from node and port in the group graph to atom index in the molecular graph
        std::unordered_map<std::string, std::unordered_map<int, int>> nodePortToAtomIndex;
        int atomCount = 0;

        // Iterate over nodes
        for (const auto& entry : nodes) {
            NodeID nodeId = entry.first;
            const Node& node = entry.second;

            std::string nodeType = node.type;
            std::vector<PortType> ports = node.ports;

            std::string smiles = nodeTypeToSmiles.at(nodeType);
            RDKit::ROMol* subGraph = RDKit::SmilesToMol(smiles);

            nodePortToAtomIndex[std::to_string(nodeId)] = std::unordered_map<int, int>();

            // Assign atom indices for ports
            for (size_t i = 0; i < ports.size(); ++i) {
                nodePortToAtomIndex[std::to_string(nodeId)][ports[i]] = atomCount + nodeTypePortToIndex.at(nodeType).at(ports[i]);
            }

            atomCount += subGraph->getNumAtoms();
            delete subGraph;
        }

        // Mapping from node and subgraph indices to molecular graph indices
        int atomId = -1;
        std::unordered_map<std::string, std::unordered_map<int, int>> nodeSubGraphIndicesToMolecularGraphIndices;
        for (const auto& entry : nodes) {
            NodeID nodeId = entry.first;
            const Node& node = entry.second;

            nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)] = std::unordered_map<int, int>();
            std::string smiles = nodeTypeToSmiles.at(node.type);
            RDKit::ROMol* subGraph = RDKit::SmilesToMol(smiles);

            for (auto atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
                atomId++;
                nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)][(*atom)->getIdx()] = atomId;
            }

            delete subGraph;
        }

        // Add atoms and bonds from the subgraphs to the molecular graph
        atomId = -1;
        for (const auto& entry : nodes) {
            NodeID nodeId = entry.first;
            const Node& node = entry.second;

            std::string smiles = nodeTypeToSmiles.at(node.type);
            RDKit::ROMol* subGraph = RDKit::SmilesToMol(smiles);

            for (RDKit::ROMol::AtomIterator atom = subGraph->beginAtoms(); atom != subGraph->endAtoms(); ++atom) {
                atomId++;
                RDKit::Atom* newAtom = new RDKit::Atom(**atom);
                molecularGraph->addAtom(newAtom, true);
            }

            // Add bonds from the subgraph
            for (RDKit::ROMol::BondIterator bond = subGraph->beginBonds(); bond != subGraph->endBonds(); ++bond) {
                RDKit::Bond* newBond = new RDKit::Bond(**bond);
                int atomIdx1 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)][(*bond)->getBeginAtomIdx()];
                int atomIdx2 = nodeSubGraphIndicesToMolecularGraphIndices[std::to_string(nodeId)][(*bond)->getEndAtomIdx()];
                (*molecularGraph).addBond(atomIdx1, atomIdx2, newBond->getBondType());
            }

            delete subGraph;
        }

        // Add bonds from the group graph to the molecular graph
        for (const auto& edge : edges) {
            NodeID from = std::get<0>(edge);
            PortType fromPort = std::get<1>(edge);
            NodeID to = std::get<2>(edge);
            PortType toPort = std::get<3>(edge);

            int fromAtom = nodePortToAtomIndex[std::to_string(from)][fromPort];
            int toAtom = nodePortToAtomIndex[std::to_string(to)][toPort];

            (*molecularGraph).addBond(fromAtom, toAtom, RDKit::Bond::SINGLE);
        }

        return molecularGraph;
    }
};

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
//     const std::unordered_map<std::string, std::unordered_map<int, int>>& nodeTypePortToIndex = {
//         {"N", {{0, 0}}},
//         {"CO", {{0, 0}, {1, 0}}},
//         {"CC", {{0, 0}, {1, 0}, {2, 1}, {3, 1}}}
//     };
//     const std::unordered_map<std::string, std::string>& nodeTypeToSmiles = {
//         {"N", "N"},
//         {"CO", "CO"},
//         {"CC", "C=C"}
//     };

//     GroupGraph gG(node_types);
//     gG.addNode(0, "N");
//     gG.addNode(1, "CC");

//     gG.addEdge(0, 0, 1, 0);

//     // gG.printGraph();

//     RDKit::ROMol* moleGraph = gG.toMolecularGraph(nodeTypeToSmiles, nodeTypePortToIndex);

//     std::string smiles = RDKit::MolToSmiles(*moleGraph);
//     std::cout << "SMILES: " << smiles << std::endl;

//     return 0;
// }
