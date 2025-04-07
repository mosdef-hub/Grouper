import pytest
from rdkit import Chem

from Grouper import Atom, AtomGraph, Group, GroupGraph, exhaustive_generate
from Grouper.tests.base_test import BaseTest
from Grouper.counting import NodeColoredGraph
import numpy as np
import sympy
from sympy.combinatorics.perm_groups import PermutationGroup
import os
from typing import List, Tuple, Set
import networkx as nx
import igraph
import itertools

class TestCounting(BaseTest):
    @staticmethod
    def parse_vcolg_txt(line: str, node_defs: Set[int]) -> Tuple[int, List[int], List[Tuple[int, int]]]:
        split_pos = line.find("  ")
        if split_pos == -1:
            raise ValueError("Invalid nauty output line...")

        node_description_str = line[:split_pos]
        edge_description_str = line[split_pos + 2:]

        node_description = node_description_str.split()
        edge_description = edge_description_str.split()

        edge_list = []
        for i in range(0, len(edge_description), 2):
            edge_list.append((int(edge_description[i]), int(edge_description[i + 1])))

        n_vertices = int(node_description[0])
        colors = [int(x) for x in node_description[2:]]

        max_color = max(colors, default=0)
        if len(node_defs) < max_color + 1:
            raise ValueError("Number of nodes in node_defs does not match the number of nodes in the nauty_output_file...")

        return n_vertices, colors, edge_list
    
    def test_graph_init(self):
        pass
    def test_graph_visualization(self):
        pass

    def test_benzene_substitution(self):
        edge_list = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
        node_colors = [0] * 6 # All nodes are the same color
        graph = NodeColoredGraph(edge_list, node_colors)
        aut_group = np.array(graph.graph.get_automorphisms_vf2())
        permutation_group = PermutationGroup(aut_group)
        cycle_index_polynomial = graph.get_cycle_index(permutation_group)
        print(cycle_index_polynomial)
        n_types = 2
        n_nodes = len(node_colors)
        types = sympy.symbols(f"t1:{n_types+1}")
        p_g = sympy.factor(
            cycle_index_polynomial.subs(
                [
                    (
                        sympy.symbols(f"s_{i}"),
                        np.sum([types[j] ** i for j in range(n_types)]),
                    )
                    for i in range(1, n_nodes + 1)
                ]
            )
        )
        print(p_g)
        n_unique_graphs = graph.apply_enumeration_theorem(
            cycle_index_polynomial, 2
        )
        assert n_unique_graphs == 13, f"Expected 13 substitutional isomers for benzene, got {n_unique_graphs}"

    def test_atom_graph_count(self):
        n_nodes = 3
        n_atoms = 2

        # os.system(f"geng -c {n_nodes} > geng_output.txt")
        # with open("geng_output.txt", "r") as f:
        #     lines = f.readlines()
        # os.system("rm geng_output.txt")
        # []
        # for line in lines:

        graphs = [b"BW", b"Bw"] # these are the only graphs with 3 nodes: the triangle and the line
        total_unique_graphs = 0
        for g in graphs:
            G = nx.from_graph6_bytes(g)
            iG = igraph.Graph.from_networkx(G)
            aut_group = np.array(iG.get_automorphisms_vf2())
            dummy_graph = NodeColoredGraph(G.edges, [0,0,0])
            dummy_graph.graph = iG
            permutation_group = PermutationGroup(aut_group)
            print("Permutation group:", permutation_group)
            cycle_index_polynomial = dummy_graph.get_cycle_index(permutation_group)
            n_unique_graphs = dummy_graph.apply_enumeration_theorem(
                cycle_index_polynomial, n_atoms
            )
            total_unique_graphs += n_unique_graphs
        
        assert total_unique_graphs == 10, f"Expected 10 unique graphs, got {total_unique_graphs}"
        print("Total unique graphs:", total_unique_graphs)
        G = nx.from_graph6_bytes(graphs[1])
        iG = igraph.Graph.from_networkx(G)
        dummy_graph = NodeColoredGraph(G.edges, [0,0,0])
        atom_types = {'t1': 'C', 't2': 'O'}
        inventory = dummy_graph.generate_inventory(atom_types)
        print("Inventory:", inventory)

    def test_exhaustive_generate_atom_graph(self):

        def is_valid_molecule(edge_list, id_to_node, colors):
            """
            Check if the graph is valid based on the node definitions.
            """
            try:
                aG = AtomGraph()
                for color in colors:
                    aG.add_node(id_to_node[color].pattern)
                for (src,dst) in edge_list:
                    aG.add_edge(src, dst)
                return True
            except Exception as e:
                return False

        node_types = {
            Group("oxygen", "O", [0, 0]),
            Group("carbon", "C", [0, 0, 0,0]),
            Group("nitrogen", "N", [0, 0,0]),
        }
        id_to_node = {i:node for i, node in enumerate(node_types)}

        exhausted_space = exhaustive_generate(
            n_nodes=3,
            node_defs=node_types,
            num_procs=1,
            vcolg_output_file="",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
        )
        total_unique_graphs = 0
        # Count number of lines in vcolg output file
        with open("vcolg_out.txt", "r") as f:
            lines = f.readlines()
            print(lines)
            for line in lines:
                n_vertices, colors, edge_list = self.parse_vcolg_txt(line, set(range(len(node_types))))
                if is_valid_molecule(edge_list, id_to_node, colors):
                    total_unique_graphs += 1

        assert len(exhausted_space) == total_unique_graphs, f"Expected {len(exhausted_space)} unique graphs, got {total_unique_graphs}"


    def test_group_graph_count(self):

        def compute_edge_permutations(edge_orbits):
            """
            Computes edge permutations given edge orbits.

            :param edge_orbits: List of integers where each unique number represents an orbit.
            :return: List of edge permutations.
            """
            orbit_map = {}  # Maps orbit ID to edge indices
            for i, orbit in enumerate(edge_orbits):
                if orbit not in orbit_map:
                    orbit_map[orbit] = []
                orbit_map[orbit].append(i)

            edge_permutations = [list(range(len(edge_orbits)))]  # Start with identity permutation

            # Generate permutations within each orbit group
            for orbit_group in orbit_map.values():
                if len(orbit_group) > 1:  # Only permute if there are multiple edges in an orbit
                    new_perms = []
                    for perm in itertools.permutations(orbit_group):
                        perm_map = {old: new for old, new in zip(orbit_group, perm)}
                        permuted_edges = [perm_map[i] if i in perm_map else i for i in range(len(edge_orbits))]
                        new_perms.append(permuted_edges)
                    edge_permutations.extend(new_perms)

            return np.array(edge_permutations)

        node_types = {
            0: Group("benzene", "C1=CC=CC=C1", [0,1,2,3,4,5]),
            1: Group("ester", "C(O)=O", [0,1])
        }
        # node_colored_graphs = [ # these the node colored graphs of the 3 node line graph with 2 colors 
        #     "3 2 0 0 0  0 2 1 2",
        #     "3 2 0 0 1  0 2 1 2",
        #     "3 2 1 0 0  0 2 1 2",
        #     "3 2 1 0 1  0 2 1 2",
        #     "3 2 1 1 0  0 2 1 2",
        #     "3 2 1 1 1  0 2 1 2",
        # ]
        node_colored_graphs = [
            "3 3 0 0 0  0 1 0 2 1 2"
        ]

        # Calculate possible maximum degree for each node
        # Calculate the possible connections for each node given the maximum degree
        # Apply the enumeration theorem to count the number of unique graphs
        # Apply list constraints to filter

        # solutions for b-b-b
        # solution = [b0_b1**2, b0_b2**2, b0_b3**2, b1_b2**2, b1_b3**2, b2_b3**2, b0_b0*b0_b1, b0_b0*b0_b2, b0_b0*b0_b3, b0_b0*b1_b1, b0_b0*b1_b2, b0_b0*b1_b3, b0_b0*b2_b2, b0_b0*b2_b3, b0_b0*b3_b3, b0_b1*b0_b2, b0_b1*b0_b3, b0_b1*b1_b1, b0_b1*b1_b2, b0_b1*b1_b3, b0_b1*b2_b2, b0_b1*b2_b3, b0_b1*b3_b3, b0_b2*b0_b3, b0_b2*b1_b1, b0_b2*b1_b2, b0_b2*b1_b3, b0_b2*b2_b2, b0_b2*b2_b3, b0_b2*b3_b3, b0_b3*b1_b1, b0_b3*b1_b2, b0_b3*b1_b3, b0_b3*b2_b2, b0_b3*b2_b3, b0_b3*b3_b3, b1_b1*b1_b2, b1_b1*b1_b3, b1_b1*b2_b2, b1_b1*b2_b3, b1_b1*b3_b3, b1_b2*b1_b3, b1_b2*b2_b2, b1_b2*b2_b3, b1_b2*b3_b3, b1_b3*b2_b2, b1_b3*b2_b3, b1_b3*b3_b3, b2_b2*b2_b3, b2_b2*b3_b3, b2_b3*b3_b3]
        for line in node_colored_graphs[:1]:
            n_vertices, colors, edge_list = self.parse_vcolg_txt(line, {0, 1})

            G = NodeColoredGraph(edge_list, colors)
            g = GroupGraph()
            for i in range(n_vertices):
                g.add_node(node_types[colors[i]].type, node_types[colors[i]].pattern, node_types[colors[i]].hubs)

            # print("Edge permutations:", edge_permutations)
            # permutation_group = PermutationGroup(edge_permutations)
            # print("Permutation group:", permutation_group)

            # Compute cycle index polynomial and apply the enumeration theorem
            cycle_index_polynomial = G.get_cycle_index(permutation_group)
            possible_graphs = G.apply_enumeration_theorem(cycle_index_polynomial, len(edge_types)+1)
            inventory = G.generate_inventory({f't{i+1}':e for i,e in enumerate(edge_types.keys())}, cycle_index_polynomial, 2)
            print("Inventory:", inventory)
            filtered_inventory = G.filter_inventory(inventory, edge_types, colors=colors, edge_list=edge_list)
            print("Filtered inventory:", filtered_inventory)
        print("Total unique graphs:", len(filtered_inventory))
        assert True