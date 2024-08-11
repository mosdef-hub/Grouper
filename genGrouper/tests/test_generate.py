from mbuild.tests.base_test import BaseTest
from genGrouper import GroupGraph
from genGrouper.generate import generate_group_graph_space
# from pysmiles import write_smiles
# import pytest

class TestGroupGraph(BaseTest):
    def test_chiral_group_graph_generation(self):
        # Define node types with ports
        node_types = {
            'CC': ['C11', 'C12', 'C21', 'C22',],
            'OH': ['O1'],
        }
        node_types_to_smiles = {
            'CC': 'CC',
            'OH': 'O',
        }
        node_port_to_atom_index = {
            'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
            'OH': {'O1': 0},
        }
        out = generate_group_graph_space(3, node_types, node_port_to_atom_index)

        def check_if_graph_has_ports(graph, ports):
            all_ports = set()
            for e in graph.edges(data=True):
                for p in e[2]['ports']:
                    all_ports.add(p[0].split('.')[1])
                    all_ports.add(p[1].split('.')[1])
            return all_ports == ports

        cc_2_oh = []
        for g in out:
            n_OH, n_CC = 0, 0
            for name, info in g.nodes(data=True):
                if info['type'] == 'OH':
                    n_OH += 1
                if info['type'] == 'CC':
                    n_CC += 1
            if n_OH == 2 and n_CC == 1:
                cc_2_oh.append(g)

        # cis TODO verifiy actually cis 
        assert any([True if check_if_graph_has_ports(g, ports={'C11', 'C22', 'O1'}) else False for g in cc_2_oh])

        # trans
        assert any([True if check_if_graph_has_ports(g, ports={'C11', 'C21', 'O1'}) else False for g in cc_2_oh])

        # terminal alkene
        assert any([True if check_if_graph_has_ports(g, ports={'C11', 'C12', 'O1'}) else False for g in cc_2_oh])

