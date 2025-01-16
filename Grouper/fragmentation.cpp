#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <nauty/nauty.h>

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <iostream>  // For debugging output


#include "dataStructures.hpp"


struct TupleHash {
    template <typename T1, typename T2, typename T3>
    std::size_t operator ()(const std::tuple<T1, T2, T3>& t) const {
        auto h1 = std::hash<T1>{}(std::get<0>(t)); 
        auto h2 = std::hash<T2>{}(std::get<1>(t)); 
        auto h3 = std::hash<T3>{}(std::get<2>(t)); 

        // Combine the hashes using a better mix (commonly used strategy)
        h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h2
        h1 ^= h3 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h3
        return h1;
    }
};

namespace std {
    template <>
    struct hash<std::tuple<char, char, int>> {
        std::size_t operator()(const std::tuple<char, char, int>& t) const {
            auto h1 = std::hash<char>{}(std::get<0>(t));
            auto h2 = std::hash<char>{}(std::get<1>(t));
            auto h3 = std::hash<int>{}(std::get<2>(t));

            // Combine the hashes using a better mix (commonly used strategy)
            h1 ^= h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h2
            h1 ^= h3 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2);  // Mixing h1 and h3
            return h1;
        }
    };
}

 std::vector<int> _findMatchingIndices(const int* orbits, int targetOrbit, int numAtoms) {
        std::vector<int> matchingIndices;

        for (int i = 0; i < numAtoms; ++i) {
            if (orbits[i] == targetOrbit) {
                matchingIndices.push_back(i);
            }
        }

        return matchingIndices;
}

std::unordered_map<GroupGraph::Node, std::unordered_set<GroupGraph::Node>> determineNodeComposition(const std::unordered_map<std::string, GroupGraph::Node>& nodeDefs){
    std::unordered_map<GroupGraph::Node, std::unordered_set<GroupGraph::Node>> nodeComposition;
    for (const auto& nodeDef1 : nodeDefs) {
        const GroupGraph::Node& nodeData = nodeDef1.second;
        const std::string& smiles = nodeData.smiles; // Need to use smiles because rdkit can't actually do substructure matching with smarts
        std::unique_ptr<RDKit::ROMol> original(RDKit::SmilesToMol(smiles));

        for (const auto& nodeDef2 : nodeDefs) {
            const std::string& smarts2 = nodeDef2.first;
            const GroupGraph::Node& nodeData2 = nodeDef2.second;
            std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts2));
            std::vector<RDKit::MatchVectType> matches;
            RDKit::SubstructMatch(*original, *query, matches);
            if (!matches.empty()) {
                nodeComposition[nodeData].insert(nodeData2);
            }
        }
        if (nodeComposition[nodeData].empty()) {
            nodeComposition[nodeData].insert(nodeData);
        }
    }
    return nodeComposition;
}

GroupGraph fragment(
    const std::string& smiles, 
    const std::unordered_map<std::string, GroupGraph::Node>& nodeDefs
) {

    GroupGraph groupGraph;

    // Parse the SMILES string to an RDKit molecule
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));

    // Error handling
    if (!mol) {
        throw std::invalid_argument("Invalid SMILES string");
    }
    for (const auto& nodeDef: nodeDefs) {
        std::unique_ptr<RDKit::ROMol> subgraph(RDKit::SmartsToMol(nodeDef.first));
        if (!subgraph) {
            throw std::invalid_argument("Invalid SMARTS pattern: " + nodeDef.first);
        }
    }

    // Create a vector of node types sorted by the number of atoms in descending order
    std::vector<std::pair<std::string, GroupGraph::Node>> sortedNodeDefs(nodeDefs.begin(), nodeDefs.end());
    std::sort(sortedNodeDefs.begin(), sortedNodeDefs.end(), 
        [](const auto& a, const auto& b) { 
            std::unique_ptr<RDKit::ROMol> molA(RDKit::SmartsToMol(a.first));
            std::unique_ptr<RDKit::ROMol> molB(RDKit::SmartsToMol(b.first));
            return molA->getNumAtoms() > molB->getNumAtoms(); 
        });

    // Maps to keep track of atom indices and corresponding information
    std::unordered_map<int, GroupGraph::NodeIDType> atomToNodeid; // Maps atom indices to node IDs
    std::unordered_map<int, std::vector<int>> atomToPorts;  // atom to matching port indices
    std::unordered_map<int, std::string> atomToSmarts;  // Maps atom indices to SMARTS patterns
    std::unordered_map<std::string, std::vector<int>> smartsToOrbits;  // Maps Smarts to orbit indices
    std::unordered_map<int, std::vector<int>> atomToOrbit;  // Maps atom indices to other atom indices in the same orbit
    std::unordered_map<int, bool> atomIndividuation;  // Maps atom indices to whether differentiated from other atoms in the same orbit via a bond
    GroupGraph::NodeIDType nodeId = 0;
    for (std::size_t i = 0; i < mol->getNumAtoms(); ++i) {
        atomToOrbit[i] = std::vector<int>();
    }
    // Nauty varibles
    int n = mol->getNumAtoms(); //maximal number of nodes
    int m = SETWORDSNEEDED(n);    // Convert AtomGraph to nauty graph
    DYNALLSTAT(graph, g, g_sz);
    DYNALLOC2(graph, g, g_sz, n, m, "malloc");
    int lab[n], ptn[n], orbits[n];

    // Determine if nodes can be made by composing other nodes
    // std::unordered_map<GroupGraph::Node, std::unordered_set<GroupGraph::Node>> nodeComposition = determineNodeComposition(nodeDefs);

    for (const auto& nodeDef : sortedNodeDefs) {
        const std::string& smarts = nodeDef.first;
        const GroupGraph::Node& nodeData = nodeDef.second;
        
        std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));

        std::vector<RDKit::MatchVectType> matches;
        RDKit::SubstructMatch(*mol, *query, matches);

        if (matches.empty()) {
            continue;
        }

        if (query->getNumAtoms() != 1) {  // Need to differentitate isomorphic atoms in the molecule to enable fragmentation that is independent of atom order when assembled in the mol
            // Initialize nauty structures
            EMPTYGRAPH(g, m, n);
            static DEFAULTOPTIONS_GRAPH(options);
            statsblk stats;
            setword workspace[160];
            graph canong[g_sz];  // Canonical form
            // Convert rdkit graph to nauty graph
            for (const auto& bond : query->bonds()) { // Convert rdkit graph to nauty graph
                int beginAtomIdx = bond->getBeginAtomIdx();
                int endAtomIdx = bond->getEndAtomIdx();
                ADDONEEDGE(g, beginAtomIdx, endAtomIdx, m);
            }
            // Call nauty to compute automorphisms
            options.getcanon = TRUE;
            nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 160, m, n, canong);
            // Save the orbits of each node type
            smartsToOrbits[smarts] = std::vector<int>(orbits, orbits + query->getNumAtoms());
        }


        for (const auto& match : matches) {
            // match is a vector of pairs (query atom index, mol atom index)
            bool alreadyMatched = false;
            for (const auto& atomMatch : match) {
                int atomIdx = atomMatch.second;

                if (atomToNodeid.count(atomIdx)) {
                    alreadyMatched = true;
                    break;
                }
            }
            if (alreadyMatched) {
                continue;
            }

            // Add the node to the graph
            groupGraph.addNode(
                nodeData.ntype, 
                nodeData.smiles, 
                nodeData.hubs
            );
            GroupGraph::NodeIDType currentId = nodeId++;

            // Fill mapping of atom to other atoms in the same orbit
            if (query->getNumAtoms() != 1) {
                for (const auto& atomMatch : match) {
                    int atomIdx = atomMatch.second;
                    int orbit = smartsToOrbits[smarts][atomMatch.first];
                    std::vector<int> orbit_no_bond_types = _findMatchingIndices(orbits, orbit, n);
                    std::unordered_map<int, std::unordered_set<std::tuple<char, char, int>, TupleHash>> orbitToBondTypes;

                    for (const auto& orbitidx: orbit_no_bond_types){
                        const RDKit::Atom* atom = query->getAtomWithIdx(orbitidx);
                        auto bonds = query->atomBonds(atom);
                        std::unordered_set<std::tuple<char, char, int>> bondTypes;
                        for (const auto& bond : bonds) {
                            int beginAtomIdx = bond->getBeginAtomIdx();
                            int endAtomIdx = bond->getEndAtomIdx();
                            char beginSymbol = mol->getAtomWithIdx(beginAtomIdx)->getSymbol()[0];
                            char endSymbol = mol->getAtomWithIdx(endAtomIdx)->getSymbol()[0];
                            bondTypes.insert(std::make_tuple(beginSymbol, endSymbol, bond->getBondType()));
                            bondTypes.insert(std::make_tuple(endSymbol, beginSymbol, bond->getBondType()));
                        }
                        orbitToBondTypes[orbitidx].insert(bondTypes.begin(), bondTypes.end());
                    }
                    std::vector<int> finalOrbit;
                    // If bondtype is different between two atoms in the same orbit, individuate and change orbit
                    finalOrbit.push_back(atomIdx);
                    atomIndividuation[atomIdx] = true;

                    for (const auto& orbitidx : orbit_no_bond_types) {
                        if (orbitidx == atomIdx) {
                            continue;
                        }
                        if (orbitToBondTypes[atomIdx] == orbitToBondTypes[orbitidx]) {
                            finalOrbit.push_back(orbitidx);
                            atomIndividuation[orbitidx] = false;
                        }
                    }
                    // Update atomToOrbit with the finalOrbit
                    for (const auto& idx : finalOrbit) {
                        atomToOrbit[idx] = finalOrbit;
                        if (finalOrbit.size() > 1) {
                            atomIndividuation[idx] = false;
                        }
                    }
                }
            }
            else{
                int atomIdx = match[0].second;
                atomIndividuation[atomIdx] = true;
                atomToOrbit[atomIdx].push_back(atomIdx);
            }
            for (const auto& atomMatch : match) {
                int atomIdx = atomMatch.second;
                atomToNodeid[atomIdx] = currentId;
                atomToSmarts[atomIdx] = smarts;  // Track which SMARTS pattern the atom belongs to
                // check if atom is associated with multiple ports
                std::vector<int> hubIndices;
                for (size_t idx = 0; idx < nodeData.hubs.size(); ++idx) {
                    const auto& hub = nodeData.hubs[idx];
                    if (hub == atomMatch.first) {
                        hubIndices.push_back(idx);
                    }
                }
                atomToPorts[atomIdx] = hubIndices;
            }
        }
    }

    // If any atoms are not matched to a node, throw an error
    std::vector<int> unmatchedAtoms;
    std::string unmatchedAtomsStr = "";
    for (const auto& atom : mol->atoms()) {
        int atomIdx = atom->getIdx();
        if (!atomToNodeid.count(atomIdx)) {
            unmatchedAtoms.push_back(atomIdx);
        }
    }
    if (!unmatchedAtoms.empty()) {
        for (const auto& atomIdx : unmatchedAtoms) {
            unmatchedAtomsStr += std::to_string(atomIdx) + " ";
        }
        std::cout<<groupGraph.printGraph();
        throw std::invalid_argument("Atoms not matched to any node: " + unmatchedAtomsStr);
    }
    
    // std::cout<<"atomToOrbit: ";
    // for (const auto& [atomIdx, orbit] : atomToOrbit) {
    //     std::cout<<"atomIdx: "<<atomIdx<<" orbit: ";
    //     for (const auto& idx : orbit) {
    //         std::cout<<idx<<" ";
    //     }
    //     std::cout<<std::endl;
    // }
    // std::cout<<std::endl;

    // print all maps assoicated with atom index
    // for (int atomIdx = 0; atomIdx < mol->getNumAtoms(); atomIdx++) {
    //     printf("atomIdx: %d ports: ", atomIdx);
    //     for (const auto& port : atomToPorts[atomIdx]) {
    //         printf("%d ", port);
    //     }
    //     printf(" nodeid: %d", atomToNodeid[atomIdx]);
    //     printf(" element: %s", mol->getAtomWithIdx(atomIdx)->getSymbol().c_str());
    //     printf(" individuation: %d", atomIndividuation[atomIdx]);
    //     printf(" smarts: %s", atomToSmarts[atomIdx].c_str());
    //     printf(" atomOrbits: ");
    //     for (const auto& orbit : atomToOrbit[atomIdx]) {
    //         printf("%d ", orbit);
    //     }
    //     printf("\n");
    // }



    // Add edges based on connectivity in the original molecule
    for (const auto& bond : mol->bonds()) {
        int beginAtomIdx = bond->getBeginAtomIdx();
        int endAtomIdx = bond->getEndAtomIdx();


        GroupGraph::NodeIDType fromNode = atomToNodeid[beginAtomIdx];
        GroupGraph::NodeIDType toNode = atomToNodeid[endAtomIdx];
        std::string fromSmarts = atomToSmarts[beginAtomIdx];
        std::string toSmarts = atomToSmarts[endAtomIdx];


        // Check if both atoms are part of the same SMARTS subgraph
        if (fromNode != toNode) {
            // see if atoms are individuated
            if (atomIndividuation[beginAtomIdx] && atomIndividuation[endAtomIdx]) {
                // Get the port numbers for the atoms
                std::vector<int> fromPorts = atomToPorts[beginAtomIdx];
                std::vector<int> toPorts = atomToPorts[endAtomIdx];
                if (!fromPorts.empty() && !toPorts.empty()) {
                    int fromPort = fromPorts.back();
                    int toPort = toPorts.back();
                    atomToPorts[beginAtomIdx].pop_back();
                    atomToPorts[endAtomIdx].pop_back();
                    groupGraph.addEdge({fromNode, fromPort}, {toNode, toPort}, true);
                }
                else {
                    std::cout<<groupGraph.printGraph();
                    std::cout<<"beginAtomIdx: "<<beginAtomIdx<<" endAtomIdx: "<<endAtomIdx<<std::endl;
                        // print all maps assoicated with atom index
                    for (std::size_t atomIdx = 0; atomIdx < mol->getNumAtoms(); atomIdx++) {
                        printf("atomIdx: %zu ports: ", atomIdx);
                        for (const auto& port : atomToPorts[atomIdx]) {
                            printf("%d ", port);
                        }
                        printf(" nodeid: %d", atomToNodeid[atomIdx]);
                        printf(" element: %s", mol->getAtomWithIdx(atomIdx)->getSymbol().c_str());
                        printf(" individuation: %d", atomIndividuation[atomIdx]);
                        printf(" smarts: %s", atomToSmarts[atomIdx].c_str());
                        printf(" atomOrbits: ");
                        for (const auto& orbit : atomToOrbit[atomIdx]) {
                            printf("%d ", orbit);
                        }
                        printf("\n");
                    }

                    throw std::invalid_argument("Not enough ports for bond");
                }
            
            }
            else if (atomIndividuation[beginAtomIdx] && !atomIndividuation[endAtomIdx]){
                for (const auto& atomid : atomToOrbit[endAtomIdx]){
                    if (atomToPorts[atomid].empty()){
                        continue;
                    }
                    int toPort = atomToPorts[atomid].back();
                    int fromPort = atomToPorts[beginAtomIdx].back();
                    atomToPorts[beginAtomIdx].pop_back();
                    atomToPorts[atomid].pop_back();
                    atomToPorts[endAtomIdx] = atomToPorts[atomid];
                    atomToPorts[atomid].clear();
                    atomIndividuation[atomid] = true;
                    groupGraph.addEdge({fromNode, fromPort}, {toNode, toPort}, true);
                    break;
                }
            }
            else if (!atomIndividuation[beginAtomIdx] && atomIndividuation[endAtomIdx]){
                for (const auto& atomid : atomToOrbit[beginAtomIdx]){
                    if (atomToPorts[atomid].empty()){
                        continue;
                    }
                    int fromPort = atomToPorts[atomid].back();
                    int toPort = atomToPorts[endAtomIdx].back();
                    // Change the ports associated with beginAtomIdx
                    atomToPorts[atomid].pop_back();
                    atomToPorts[beginAtomIdx] = atomToPorts[atomid];
                    atomToPorts[atomid].clear();
                    atomToPorts[endAtomIdx].pop_back();
                    atomIndividuation[atomid] = true;
                    groupGraph.addEdge({fromNode, fromPort}, {toNode, toPort}, true);
                    break;
                }
            }
            else{
                for (const auto& atomid : atomToOrbit[beginAtomIdx]){
                    if (atomToPorts[atomid].empty()){
                        continue;
                    }
                    int fromPort = atomToPorts[atomid].back();
                    for (const auto& atomid2 : atomToOrbit[endAtomIdx]){
                        if (atomToPorts[atomid2].empty()){
                            continue;
                        }
                        int toPort = atomToPorts[atomid2].back();
                        atomToPorts[atomid].pop_back();
                        atomToPorts[atomid2].pop_back();
                        atomToPorts[beginAtomIdx] = atomToPorts[atomid];
                        atomToPorts[endAtomIdx] = atomToPorts[atomid2];
                        atomToPorts[atomid].clear();
                        atomToPorts[atomid2].clear();
                        atomIndividuation[atomid] = true;
                        atomIndividuation[atomid2] = true;
                        groupGraph.addEdge({fromNode, fromPort}, {toNode, toPort}, true);
                        break;
                    }
                }
            }
        }
    }
    return groupGraph;
}