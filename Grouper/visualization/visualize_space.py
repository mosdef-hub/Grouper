import typing as t
from typing import Dict

import matplotlib.pyplot as plt
import networkx as nx
# from grakel.kernels import Kernel
from networkx import Graph
from numpy.typing import ArrayLike
from rdkit import Chem
from rdkit.Chem import rdmolops
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np
from collections import Counter

from Grouper import GroupGraph


def chemical_space_network(
    space_set: t.Set[GroupGraph],
    kernel,
    cutoff: float,
    known_values: dict = None,
    options: dict = None,
) -> Graph:
    graphs = []
    for g in space_set:
        smile = g.to_smiles()
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            adj_matrix = rdmolops.GetAdjacencyMatrix(mol)
            atom_labels = {i: atom.GetSymbol() for i, atom in enumerate(mol.GetAtoms())}
            graphs.append((adj_matrix, atom_labels))
    graph_kernel = kernel(n_iter=5, normalize=True)
    kernel_matrix = graph_kernel.fit_transform(graphs)
    G = nx.Graph()
    for i, g in enumerate(space_set):
        smile = g.to_smiles()
        if known_values is None:
            color = 0
        else:
            color = 0 if known_values[smile] is None else known_values[smile]
        G.add_node(i, label=smile, property=color)
    for i in range(len(kernel_matrix)):
        for j in range(i + 1, len(kernel_matrix)):
            if kernel_matrix[i][j] > cutoff:
                G.add_edge(i, j, weight=kernel_matrix[i][j])
    return G


def calculate_edge_density(
    space_set, kernel, cutoffs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
):
    edge_density = []
    for cutoff in cutoffs:
        G = chemical_space_network(space_set, kernel, cutoff)
        edge_density.append(nx.density(G))
    return list(zip(cutoffs, edge_density))


def calculate_modularity(
    space_set,
    kernel,
    community: str = None,
    cutoffs=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
    options: dict = None,
) -> Dict[str, ArrayLike]:
    modularity = []
    for cutoff in cutoffs:
        G = chemical_space_network(space_set, kernel, cutoff)
        if community is None:
            modularity.append(
                nx.community.modularity(
                    G, nx.community.greedy_modularity_communities(G)
                )
            )
        else:
            modularity.append(nx.community.modularity(G, community))

    return list(zip(cutoffs, modularity))


def calculate_network_properties(
    space_network, options: dict = None
) -> t.Dict[str, ArrayLike]:
    # Calculate network properties
    properties = {}
    properties["average_clustering"] = []
    properties["average_shortest_path"] = []
    properties["diameter"] = []
    properties["radius"] = []
    properties["density"] = []

    for i, node in enumerate(space_network.nodes()):
        subgraph = space_network.subgraph(
            [n for n in space_network.neighbors(node)] + [node]
        )
        properties["average_clustering"].append(nx.average_clustering(subgraph))
        properties["average_shortest_path"].append(
            nx.average_shortest_path_length(subgraph)
        )
        properties["diameter"].append(nx.diameter(subgraph))
        properties["radius"].append(nx.radius(subgraph))
        properties["density"].append(nx.density(subgraph))

    return properties


def plot_feature_similarity(
    feature_set: t.Mapping[str, ArrayLike],
    kernel: t.Callable[[ArrayLike, ArrayLike], float],
    cutoff: float,
    options: dict,
) -> plt.Figure:
    pass

def calculate_centrality(space_network: Graph) -> Dict[str, Dict]:
    """Compute various centrality measures."""
    return {
        "degree_centrality": nx.degree_centrality(space_network),
        "betweenness_centrality": nx.betweenness_centrality(space_network),
        "closeness_centrality": nx.closeness_centrality(space_network),
        "eigenvector_centrality": nx.eigenvector_centrality(space_network),
    }


def visualize_chemical_space(kernel_matrix: ArrayLike, method: str = "PCA") -> plt.Figure:
    """Project the chemical space into a 2D representation using PCA or t-SNE."""
    if method == "PCA":
        reduced = PCA(n_components=2).fit_transform(kernel_matrix)
    elif method == "tSNE":
        reduced = TSNE(n_components=2, perplexity=30, metric="precomputed").fit_transform(1 - kernel_matrix)
    else:
        raise ValueError("Unsupported dimensionality reduction method")

    fig, ax = plt.subplots()
    ax.scatter(reduced[:, 0], reduced[:, 1])
    ax.set_title(f"{method} Projection of Chemical Space")
    return fig


def calculate_diversity_index(kernel_matrix: ArrayLike) -> float:
    """Compute molecular diversity index based on pairwise similarity."""
    return np.mean(1 - kernel_matrix[np.triu_indices_from(kernel_matrix, k=1)])


def find_shortest_paths(space_network: Graph, source: str, target: str) -> t.List[t.List[str]]:
    """Find shortest paths between two molecules in the chemical space."""
    return list(nx.all_shortest_paths(space_network, source=source, target=target, weight="weight"))


def calculate_substructure_frequencies(space_set: t.Set[GroupGraph]) -> Dict[str, int]:
    """Identify frequently occurring substructures in the chemical space."""
    substructure_counts = Counter()
    for g in space_set:
        for i in g.nodes.keys():
            substructure_counts[g.nodes[i].type] += 1
    return substructure_counts

def visualize_substructure_correlation(space_set: t.Set[GroupGraph], graph_property: ArrayLike, options: dict = None) -> plt.Figure:
    """Visualize the correlation between substructures in the chemical space."""
    from scipy.stats import pearsonr
    def to_ordered_vector(g, order):
        node_hist = g.to_vector()
        ordered_vector = []
        for o in order:
            if o not in node_hist:
                ordered_vector.append(0)
            else:
                ordered_vector.append(node_hist[o])
        return ordered_vector
    
    nodes_types = set()
    for g in space_set:
        for i in g.nodes.keys():
            nodes_types.add(g.nodes[i].type)
    order = list(nodes_types)
    # Calculate Pearson correlation coefficient
    data = []
    for g in space_set:
        ordered_vector = to_ordered_vector(g, order)
        data.append(ordered_vector)

    data = np.array(data)
    properties = np.array(graph_property)
    correlations = []
    for i in range(data.shape[1]):
        corr, _ = pearsonr(data[:, i], properties)
        correlations.append(corr)
    fig, ax = plt.subplots()
    ax.bar(range(len(correlations)), correlations)
    ax.set_xticks(range(len(correlations)))
    ax.set_xticklabels(order, rotation=45)
    ax.set_title("Substructure Correlation")
    return fig

