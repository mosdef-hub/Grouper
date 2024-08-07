from mbuild.tests.base_test import BaseTest
from molGrouper.group_graph import GroupGraph
from molGrouper.generate import generate_group_graph_space
from pysmiles import write_smiles
import pytest

class TestGroupGraph(BaseTest):

    @pytest.fixture(autouse=True)
    def test_creation(self):
        # Define node types with ports
        self.node_types = {
            'NH2': ['N1'], # amine
            'CO': ['C1', 'C2'], # carbonyl
            'CC': ['C11', 'C12', 'C21', 'C22'], # alkene
        }
        self.node_type_to_smiles = {
            'NH2': 'N',
            'CO': 'C=O',
            'CC': 'C=C',
        }
        self.node_port_to_atom_index = {
            'NH2': {'N1': 0},
            'CO': {'C1': 0, 'C2': 0},
            'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
        }

        # Create an instance of GroupGraph for testing
        self.graph = GroupGraph(self.node_types)

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
        out = generate_group_graph_space(3, node_types)

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

