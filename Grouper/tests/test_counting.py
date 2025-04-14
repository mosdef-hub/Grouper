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
from tqdm import tqdm

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
            
        all_atoms = {
            Group("oxygen", "O", [0, 0]),
            Group("carbon", "C", [0, 0, 0, 0]),
            Group("nitrogen", "N", [0, 0, 0]),
            Group("fluorine", "F", [0])
        }
        def powerset(iterable):
            """
            powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
            """
            xs = list(iterable)
            # note we return an iterator rather than a list
            return itertools.chain.from_iterable(itertools.combinations(xs,n) for n in range(len(xs)+1))



        for combo in powerset(all_atoms):
            if len(combo) == 0:
                continue
            combo = set(combo)

            id_to_node = {i:node for i, node in enumerate(combo)}

            exhausted_space = exhaustive_generate(
                n_nodes=4,
                node_defs=combo,
                num_procs=-1,
                vcolg_output_file="",
                positive_constraints={},
                negative_constraints=set(),
                config_path="",
            )
            total_unique_graphs = 0
            # Count number of lines in vcolg output file
            with open("vcolg_out.txt", "r") as f:
                lines = f.readlines()
                for line in lines:
                    n_vertices, colors, edge_list = self.parse_vcolg_txt(line, set(range(len(combo))))
                    if is_valid_molecule(edge_list, id_to_node, colors):
                        total_unique_graphs += 1

            assert len(exhausted_space) == total_unique_graphs, f"Expected {len(exhausted_space)} unique graphs, got {total_unique_graphs}"


    def test_group_graph_count(self):
        for n_nodes in range(2, 5):
            node_types = {
                # 0: Group("benzene", "C1=CC=CC=C1", [0,1,2,3,4,5]),
                # 0: Group('alkene', 'C=C', [0, 0, 1,1]),
                0: Group("ester", "C(O)=O", [0,1]), 
                # 2: Group("amine", "N", [0,0,0]),
            }
            # generate uncolored graphs
            os.system(f"geng -c {n_nodes} > geng_output.txt")
            # generate colored graphs
            os.system(f"vcolg geng_output.txt -T -m{len(node_types)} > vcolg_out.txt")
            # calculate hub reps
            pattern_inventory = set()
            # determine every possible port pairing for each edge
            with open("vcolg_out.txt", "r") as f:
                lines = f.readlines()
                for line in tqdm(lines, desc="Generating patterns", total=len(lines)):
                    unfiltered_edge_inventory = {} # (src,dst) -> list of all possible port pairings
                    n_vertices, colors, edge_list = self.parse_vcolg_txt(line, set(range(len(node_types))))
                    gG = GroupGraph()
                    # add nodes
                    for i in range(n_vertices):
                        gG.add_node(node_types[colors[i]].type, node_types[colors[i]].pattern, node_types[colors[i]].hubs)
                    # Generate GroupGraph with all possible bonding patterns
                    for (src, dst) in edge_list:
                        unfiltered_edge_inventory[(src, dst)] = []
                        # get all possible port pairings
                        src_ports = node_types[colors[src]].ports
                        dst_ports = node_types[colors[dst]].ports
                        combos = itertools.product(src_ports, dst_ports)
                        # print(list(combos))
                        unfiltered_edge_inventory[(src, dst)] = list(combos)

                    for ports in itertools.product(*unfiltered_edge_inventory.values()):

                        # if set(list(ports)) == set([(0,0), (0,1), (1,2)]):
                        #     import pdb
                        #     pdb.set_trace()

                        try:
                            gG.clear_edges()
                            for i, (src, dst) in enumerate(edge_list):
                                src_port, dst_port = ports[i]
                                # src_hub, dst_hub = node_types[colors[src]].hubs[src_port], node_types[colors[dst]].hubs[dst_port]
                                # print("Adding edge:", (src, src_hub), (dst, dst_hub))
                                gG.add_edge( (src, src_port), (dst, dst_port))

                            if gG.to_smiles() == "O=C1OOC(=O)C(=O)OOC1=O":
                                import pdb
                                pdb.set_trace()
                                
                            pattern_inventory.add(gG.to_smiles())
                        except Exception as e:
                            continue

            exhausted_space = exhaustive_generate(
                n_nodes=n_nodes,
                node_defs=set(node_types.values()),
                num_procs=-1,
                vcolg_output_file="",
                positive_constraints={},
                negative_constraints=set(),
                config_path="",
            )

            # find missing graphs in the exhaustive space
            missing_graphs = set()
            for graph in exhausted_space:
                if graph.to_smiles() not in pattern_inventory:
                    missing_graphs.add(graph.to_smiles())
                    print(graph)
            if missing_graphs:
                print("Missing graphs in pattern inventory:", missing_graphs)

            exhausted_space = set(g.to_smiles() for g in exhausted_space)
            missing_graphs = set()
            for graph in pattern_inventory:
                if graph not in exhausted_space:
                    missing_graphs.add(graph)
            if missing_graphs:
                print("Missing graphs in exhausted space:", missing_graphs)

            assert len(exhausted_space) == len(pattern_inventory), f"Expected {len(exhausted_space)} unique graphs, got {len(pattern_inventory)}"