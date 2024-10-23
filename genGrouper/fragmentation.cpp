#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <iostream>  // For debugging output
#include "dataStructures.hpp"

GroupGraph fragment(
    const std::string& smiles, 
    const std::unordered_map<std::string, GroupGraph::Node>& nodeDefs
) {
    GroupGraph groupGraph;

    // Parse the SMILES string to an RDKit molecule
    std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) {
        throw std::invalid_argument("Invalid SMILES string");
    }

    // Create a vector of node types sorted by the number of atoms in descending order
    std::vector<std::pair<std::string, GroupGraph::Node>> sortedNodeDefs(nodeDefs.begin(), nodeDefs.end());
    std::sort(sortedNodeDefs.begin(), sortedNodeDefs.end(), 
        [](const auto& a, const auto& b) { 
            std::unique_ptr<RDKit::ROMol> molA(RDKit::SmartsToMol(a.first));
            std::unique_ptr<RDKit::ROMol> molB(RDKit::SmartsToMol(b.first));
            return molA->getNumAtoms() > molB->getNumAtoms(); 
        });

    // Map to keep track of atom indices and corresponding node IDs
    std::unordered_map<int, GroupGraph::NodeIDType> atomMapping;
    std::unordered_map<int, std::string> atomToSmarts;  // Maps atom indices to SMARTS patterns
    GroupGraph::NodeIDType nodeId = 0;

    for (const auto& nodeDef : sortedNodeDefs) {
        const std::string& smarts = nodeDef.first;
        const GroupGraph::Node& nodeData = nodeDef.second;
        
        std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
        if (!query) {
            throw std::invalid_argument("Invalid SMARTS pattern: " + smarts);
        }

        std::vector<RDKit::MatchVectType> matches;
        RDKit::SubstructMatch(*mol, *query, matches);

        for (const auto& match : matches) {
            bool alreadyMatched = false;
            for (const auto& atomMatch : match) {
                int atomIdx = atomMatch.second;
                if (atomMapping.count(atomIdx)) {
                    alreadyMatched = true;
                    break;
                }
            }
            if (alreadyMatched) {
                continue;
            }

            groupGraph.addNode(
                nodeData.ntype, 
                RDKit::MolToSmiles(*query), 
                nodeData.hubs
            );
            GroupGraph::NodeIDType currentId = nodeId++;
            for (const auto& atomMatch : match) {
                int atomIdx = atomMatch.second;
                atomMapping[atomIdx] = currentId;
                atomToSmarts[atomIdx] = smarts;  // Track which SMARTS pattern the atom belongs to
            }
        }
    }

    // Add edges based on connectivity in the original molecule
    for (const auto& bond : mol->bonds()) {
        int beginAtomIdx = bond->getBeginAtomIdx();
        int endAtomIdx = bond->getEndAtomIdx();

        if (atomMapping.count(beginAtomIdx) && atomMapping.count(endAtomIdx)) {
            GroupGraph::NodeIDType fromNode = atomMapping[beginAtomIdx];
            GroupGraph::NodeIDType toNode = atomMapping[endAtomIdx];
            std::string fromSmarts = atomToSmarts[beginAtomIdx];
            std::string toSmarts = atomToSmarts[endAtomIdx];

            // std::cout<< "Checking edge between atoms: " << beginAtomIdx << " and " << endAtomIdx << std::endl;
            // std::cout << "From SMARTS: " << fromSmarts << " To SMARTS: " << toSmarts << std::endl;
            // std::cout << "From Node: " << fromNode << " To Node: " << toNode << std::endl;
            
            // Check if both atoms are part of the same SMARTS subgraph
            if (fromSmarts != toSmarts && fromNode != toNode) {
                std::cout << "Adding edge between nodes: " << fromNode << " and " << toNode << std::endl;
                groupGraph.addEdge({fromNode, 0}, {toNode, 0}, false);  // Assuming port 0 for simplicity
            }
        }
    }

    return groupGraph;
}
