#ifndef GROUP_GRAPH_SAMPLER_HPP
#define GROUP_GRAPH_SAMPLER_HPP

#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <unordered_set>
#include <random>
#include "dataStructures.hpp"

// Monte Carlo Sampler for GroupGraph
class GroupGraphSampler {
public:
    GroupGraphSampler(double temperature, int maxIterations);

    // Perform Monte Carlo sampling
    void sample(GroupGraph& initialGraph);

    // Set temperature
    void setTemperature(double temperature);

    // Set maximum iterations
    void setMaxIterations(int maxIterations);

private:
    double temperature;
    int maxIterations;
    std::mt19937 rng; // Random number generator

    // Calculate energy using RDKit logP
    double calculateEnergy(const GroupGraph& graph);

    // Perturb the graph (e.g., add/remove edges or nodes)
    GroupGraph perturbGraph(const GroupGraph& graph);

    // Accept or reject the new graph based on Metropolis criterion
    bool acceptMove(double energyOld, double energyNew);

    int getRandomInt(int min, int max);

    int getRandomPort(const std::vector<int>& ports);


};

#endif // GROUP_GRAPH_SAMPLER_HPP
