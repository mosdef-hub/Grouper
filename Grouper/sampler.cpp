#include "sampler.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <random>
#include <GraphMol/Descriptors/Crippen.h>

// Constructor
GroupGraphSampler::GroupGraphSampler(double temperature, int maxIterations)
    : temperature(temperature), maxIterations(maxIterations), rng(std::random_device{}()) {}

// Perform Monte Carlo sampling
void GroupGraphSampler::sample(GroupGraph& initialGraph) {
    GroupGraph currentGraph = initialGraph;
    double currentEnergy = calculateEnergy(currentGraph);

    for (int i = 0; i < maxIterations; ++i) {
        GroupGraph newGraph = perturbGraph(currentGraph);
        double newEnergy = calculateEnergy(newGraph);

        if (acceptMove(currentEnergy, newEnergy)) {
            currentGraph = newGraph;
            currentEnergy = newEnergy;
        }

        std::cout << "Iteration " << i + 1 << ": Energy = " << currentEnergy << "\n";
    }

    initialGraph = currentGraph;
}

// Set temperature
void GroupGraphSampler::setTemperature(double temperature) {
    if (temperature <= 0) {
        throw std::invalid_argument("Temperature must be positive.");
    }
    this->temperature = temperature;
}

// Set maximum iterations
void GroupGraphSampler::setMaxIterations(int maxIterations) {
    if (maxIterations <= 0) {
        throw std::invalid_argument("Maximum iterations must be positive.");
    }
    this->maxIterations = maxIterations;
}

// Calculate energy using RDKit logP
double GroupGraphSampler::calculateEnergy(const GroupGraph& graph) {
    std::string smiles = graph.toSmiles();
    RDKit::ROMol* mol = RDKit::SmilesToMol(smiles);
    if (!mol) {
        throw std::runtime_error("Failed to parse SMILES: " + smiles);
    }

    double logP = 0.0;
    double molarRefractivity = 0.0;
    RDKit::Descriptors::calcCrippenDescriptors(*mol, logP, molarRefractivity);
    delete mol;

    return -logP;
}

// Perturb the graph by adding nodes and modifying ports
GroupGraph GroupGraphSampler::perturbGraph(const GroupGraph& graph) {
    GroupGraph newGraph = graph;

    if (getRandomInt(0, 1) == 0) {
        // Add a new node
        std::string newType = "C";
        std::string newSmarts = "C";
        std::vector<GroupGraph::NodeIDType> newHubs;
        if (!newGraph.nodes.empty()) {
            // Connect the new node to an existing random node
            int randomNodeIndex = getRandomInt(0, newGraph.nodes.size() - 1);
            auto it = std::next(newGraph.nodes.begin(), randomNodeIndex);
            newHubs.push_back(it->first); // Add the selected node ID as a hub
        }
        newGraph.addNode(newType, newSmarts, newHubs);
    } else {
        // Modify ports of an existing node
        if (!newGraph.nodes.empty()) {
            int randomNodeIndex = getRandomInt(0, newGraph.nodes.size() - 1);
            auto it = std::next(newGraph.nodes.begin(), randomNodeIndex);
            GroupGraph::NodeIDType nodeId = it->first;
            GroupGraph::Node& node = it->second;

            if (!node.ports.empty()) {
                // Modify a random port
                int randomPortIndex = getRandomInt(0, node.ports.size() - 1);
                node.ports[randomPortIndex] = getRandomInt(0, 100); // Change to a random new port value
            }
        }
    }

    return newGraph;
}

// Accept or reject the new graph based on Metropolis criterion
bool GroupGraphSampler::acceptMove(double energyOld, double energyNew) {
    if (energyNew < energyOld) {
        return true;
    }

    double acceptanceProbability = std::exp((energyOld - energyNew) / temperature);
    double randomValue = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
    return randomValue < acceptanceProbability;
}

// Generate a random integer in the range [min, max]
int GroupGraphSampler::getRandomInt(int min, int max) {
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng);
}

// Generate a random port from the available ports
int GroupGraphSampler::getRandomPort(const std::vector<int>& ports) {
    if (ports.empty()) {
        throw std::runtime_error("No ports available to choose from.");
    }
    return ports[getRandomInt(0, ports.size() - 1)];
}
