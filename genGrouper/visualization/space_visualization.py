import typing as t
import matplotlib.pyplot as plt
from numpy.typing import ArrayLike
from grakel.kernels import Kernel
from networkx import Graph
import networkx as nx
from rdkit import Chem
from rdkit.Chem import rdmolops


def chemical_space_network(space_set : t.Set[str], kernel: Kernel, cutoff: float, known_values: dict = None, options: dict = None) -> Graph:
    # Convert SMILES to graphs using RDKit
    graphs = []
    for smile in space_set:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            adj_matrix = rdmolops.GetAdjacencyMatrix(mol)
            atom_labels = {i:atom.GetSymbol() for i, atom in enumerate(mol.GetAtoms())}
            graphs.append((adj_matrix, atom_labels))
    # Initialize Weisfeiler-Lehman graph kernel
    graph_kernel = kernel(n_iter=5, normalize=True)
    # Calculate the kernel matrix
    kernel_matrix = graph_kernel.fit_transform(graphs)
    # Create a network graph
    G = nx.Graph()
    # Add nodes to the graph
    for i, smile in enumerate(space_set):
        G.add_node(i, label=smile, known_value=known_values.get(smile, None) if known_values is not None else None)
    # Add edges based on kernel cutoff
    for i in range(len(kernel_matrix)):
        for j in range(i+1, len(kernel_matrix)):
            if kernel_matrix[i][j] > cutoff:
                G.add_edge(i, j, weight=kernel_matrix[i][j])
    return G

def calculate_edge_density(space_set, kernel, cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, .7, .8, .9]):
    edge_density = []
    for cutoff in cutoffs:
        G = chemical_space_network(space_set, kernel, cutoff)
        edge_density.append(nx.density(G))

    return list(zip(cutoffs, edge_density))

from typing import Dict, List
import networkx as nx
import numpy as np
from numpy.typing import ArrayLike

def calculate_modularity(space_set, kernel, community: str = None, cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, .7, .8, .9], options: dict = None) -> Dict[str, ArrayLike]:
    modularity = []
    for cutoff in cutoffs:
        G = chemical_space_network(space_set, kernel, cutoff)
        if community is None:
            modularity.append(nx.community.modularity(G, nx.community.greedy_modularity_communities(G)))
        else:
            modularity.append(nx.community.modularity(G, community))

    return list(zip(cutoffs, modularity))



def calculate_network_properties(space_network, options: dict = None) -> t.Dict[str, ArrayLike]:
    # Calculate network properties
    properties = {}
    properties['average_clustering'] = []
    properties['average_shortest_path'] = []
    properties['diameter'] = []
    properties['radius'] = []
    properties['density'] = []


    for i, node in enumerate(space_network.nodes()):
        subgraph = space_network.subgraph([n for n in space_network.neighbors(node)] + [node])
        properties['average_clustering'].append(nx.average_clustering(subgraph))
        properties['average_shortest_path'].append(nx.average_shortest_path_length(subgraph))
        properties['diameter'].append(nx.diameter(subgraph))
        properties['radius'].append(nx.radius(subgraph))
        properties['density'].append(nx.density(subgraph))

    return properties

    

def plot_feature_similarity(feature_set: t.Mapping[str, ArrayLike], kernel: t.Callable[[ArrayLike, ArrayLike], float], cutoff: float, options: dict) -> plt.Figure:
    pass
