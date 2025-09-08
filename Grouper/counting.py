# We took inspiration from https://github.com/killiansheriff/Polya


import itertools

import igraph
import matplotlib.pyplot as plt
import numpy as np
import sympy
import sympy as sp
from sympy.combinatorics.perm_groups import PermutationGroup


class GraphBase:
    def generate_permutations(self, values, coords):
        signs = list(itertools.product([1, -1], repeat=len(values)))
        permutations = []
        for sign in signs:
            value = tuple(val * s for val, s in zip(values, sign))
            permutations.extend(list((itertools.permutations(value, 3))))
        coords = np.vstack((coords, np.array(list(set(permutations)))))
        return coords

    def get_edges(self, vertexpositions, nn_dst, atol=0.1):
        # Subtract each point from all the other points
        diff = vertexpositions[:, np.newaxis, :] - vertexpositions[np.newaxis, :, :]
        # Compute the norm of the differences along the last axis
        distances = np.linalg.norm(diff, axis=-1)
        # Add edges to the graph
        mask = np.isclose(distances, nn_dst, atol=atol)
        edge_index_0, edge_index_1 = np.where(mask)
        edges = [(edge_index_0[i], edge_index_1[i]) for i in range(len(edge_index_0))]
        edges = np.unique(np.sort(edges), axis=0)
        return edges, distances

    def get_cycle_index(self, permutation_group):
        cycle_types = [p.cycle_structure for p in permutation_group]
        monomials = [
            np.prod(
                [sp.symbols(f"s_{ctype}") ** cycle[ctype] for ctype in cycle.keys()]
            )
            for cycle in cycle_types
        ]
        nnodes = np.sum([key * value for key, value in cycle_types[0].items()])
        group_size = len(permutation_group) + 1  # add identity
        cycle_index = np.sum(monomials) + sp.symbols("s_1") ** nnodes
        return cycle_index / group_size  # need divided size of group

    def apply_enumeration_theorem(self, cycle_index, ntypes):
        n_nodes = self.graph.vcount()
        # replace s_i with the sum of the powers of 1 for each variable and factorize
        num_orbits = sp.factor(
            cycle_index.subs(
                [
                    (sp.symbols(f"s_{i}"), sum([1**i for _ in range(ntypes)]))
                    for i in range(1, n_nodes + 1)
                ]
            )
        )
        return num_orbits

    def generate_inventory(self, pattern_types, cycle_index=None, nnodes=None):
        """
        Prints the pattern inventory using the cycle index polynomial
        for the specified number of atom types (ntypes), with optional atom type substitutions.

        Params:
            ntypes (int): Number of atom types.
            atom_types (dict, optional): A dictionary mapping atom type symbols to their names,
                                            e.g., {"t1": "C", "t2": "O"}.
        """

        # Get the number of atom types
        ntypes = len(pattern_types)
        # Get the number of nodes
        if nnodes is None:
            nnodes = self.graph.vcount()
        # Get the cycle index polynomial
        if cycle_index is None:
            cycle_index = self.get_cycle_index(
                PermutationGroup(np.array(self.graph.get_automorphisms_vf2()))
            )
        # define symbolic variables for d1 to d10
        types = sp.symbols(f"t1:{ntypes + 1}")
        # replace s_i with the sum of the powers of the d variables and factorize
        p_g = sp.factor(
            cycle_index.subs(
                [
                    (
                        sp.symbols(f"s_{i}"),
                        np.sum([types[j] ** i for j in range(ntypes)]),
                    )
                    for i in range(1, nnodes + 1)
                ]
            )
        )
        # Substitute atom type symbols from the dictionary into the polynomial
        substitutions = [
            (
                sp.symbols(f"t{i + 1}"),
                sp.symbols(pattern_types.get(f"t{i + 1}", f"t{i + 1}")),
            )
            for i in range(ntypes)
        ]
        p_g_substituted = p_g.subs(substitutions)
        return p_g_substituted

    def filter_inventory(self, inventory, edge_types, colors, edge_list):
        connections = []
        for e in edge_list:
            connection = (colors[e[0]], colors[e[1]])
            connections.append(connection)

        expanded_inventory = sympy.expand(inventory)
        filtered_inventory = []
        for expression in expanded_inventory.args:
            base, exponent = expression.as_base_exp()
            factors = base.as_ordered_factors()
            valid = True
            for factor in factors:
                connection_type = edge_types[str(factor)]
                if (
                    connection_type not in connections
                    or connection_type[::-1] not in connections
                ):
                    valid = False
                    break
            if valid:
                filtered_inventory.append(expression)
        return filtered_inventory

    def plot(self):
        # Quick viz plotting
        # import matplotlib.pyplot as plt
        fig, ax = plt.subplots(nrows=1, ncols=1)
        layout = self.graph.layout_kamada_kawai()
        igraph.plot(self.graph, layout=layout, target=ax)
        ax.set_title(self.graph_name)
        plt.axis("off")
        return fig, ax


class NodeColoredGraph(GraphBase):
    def __init__(self, edge_list, node_colors):
        super().__init__()
        self.edge_list = edge_list
        self.node_colors = node_colors
        self.graph_name = "Node Colored Graph"

        g = igraph.Graph(directed=False)
        # Add vertices based on the number of unique nodes in the edge list
        nodes = set()
        for edge in self.edge_list:
            nodes.update(edge)
        g.add_vertices(max(nodes) + 1)

        # Add edges to the graph
        for u, v in self.edge_list:
            g.add_edges([(u, v)])

        # Assign node colors
        color_dict = {
            node: color
            for node, color in zip(range(len(self.node_colors)), self.node_colors)
        }
        g.vs["color"] = [color_dict[i] for i in range(g.vcount())]

        self.graph = g


# Use cases:
# count possible atom graphs given a number of nodes
# count possible group graphs given a number of nodes

# generate all graphs
