from mbuild.tests.base_test import BaseTest
from molGrouper.group_graph import GroupGraph
from group_selfies import Group
from molGrouper.io import has_mbuild, has_torch
import pytest

class TestGroupGraph(BaseTest):

    @pytest.fixture(autouse=True)
    def test_creation(self):
        # Define node types with ports
        node_types = {
            'type1': ['port1', 'port2'],
            'type2': ['port3', 'port4'],
        }

        # Create an instance of GroupGraph for testing
        self.graph = GroupGraph(node_types)

    def test_add_node(self):
        self.graph.add_node('node1', 'type1')
        assert 'node1' in self.graph.nodes
        assert self.graph.nodes['node1']['type'] == 'type1'
        assert self.graph.nodes['node1']['ports'] == ['port1', 'port2']

        self.graph.add_node('node2', 'type1')
        assert 'node2' in self.graph.nodes
        assert self.graph.nodes['node2']['type'] == 'type1'
        assert self.graph.nodes['node2']['ports'] == ['port1', 'port2']

        self.graph.add_node('node3', 'type2')
        assert 'node3' in self.graph.nodes
        assert self.graph.nodes['node3']['type'] == 'type2'
        assert self.graph.nodes['node3']['ports'] == ['port3', 'port4']

    def test_add_edge(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type1')
        self.graph.add_node('node3', 'type2')

        self.graph.add_edge('node1', 'port1', 'node2', 'port1')
        assert ('node1', 'node2') in self.graph.edges
        assert self.graph.edges[('node1', 'node2')]['ports'][0] == ['node1.port1', 'node2.port1']

        self.graph.add_edge('node3', 'port3', 'node2', 'port2')
        assert ('node3', 'node2') in self.graph.edges
        assert self.graph.edges[('node3', 'node2')]['ports'][0] == ['node3.port3', 'node2.port2']

    def test_make_undirected(self):
        # self.graph.add_node('node1', 'type1')
        # self.graph.add_node('node2', 'type2')
        # self.graph.add_edge('node1', 'port1', 'node2', 'port3')
        # self.graph.make_undirected()

        # # Check if the graph is undirected
        # assert ('node1', 'node2') in self.graph.edges
        # assert ('node2', 'node1') in self.graph.edges
        # assert self.graph.edges['node1', 'node2']['ports'][1] == ['node2.port3', 'node1.port1']
        pass

    def test_n_free_ports(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type2')
        assert self.graph.n_free_ports('node1') ==  2

        # Connect a edge and recheck
        self.graph.add_edge('node1', 'port1', 'node2', 'port3')
        assert self.graph.n_free_ports('node1') == 1
        assert self.graph.n_free_ports('node2') == 1

    @pytest.mark.skipif(not has_mbuild, reason="mBuild package not installed")
    def test_compound_to_group_graph(self):
        import mbuild as mb
        mol = mb.load('CCCCCCCC', smiles=True) # octane molecule
        groups = [Group('c3', 'C([H])([H])([H])(*1)'), Group('c2', 'C([H])([H])(*1)(*1)')]
        graph = GroupGraph()
        group_graph = graph.from_mbuild(mol, groups)
        print(group_graph)


    def test_group_graph_to_vector(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type2')
        self.graph.add_edge('node1', 'port1', 'node2', 'port3')
        vector_form = self.graph.to_vector()
        assert vector_form == [1, 1]

    @pytest.mark.skipif(not has_torch, reason="torch package not installed")
    def test_group_graph_to_pyG(self):
        import torch
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type2')
        self.graph.add_edge('node1', 'port1', 'node2', 'port3')
        group_featurizer = lambda node: torch.tensor([1, 0])

        data = self.graph.to_PyG_Data(group_featurizer)
        assert torch.equal(data.x, torch.tensor([ [1,0], [1,0] ], dtype=torch.float32)) # node features should just be identity
        assert torch.equal(data.edge_index, torch.tensor([ [0], [1] ], dtype=torch.float32)) # graph is directed, node1 -> node2
        assert torch.equal(data.edge_attr, torch.tensor([ [1,0,1,0] ], dtype=torch.float32)) # edge features are one-hot encoded port

    def test_group_graph_to_atomic_graph(self):
        graph = GroupGraph(
            node_types = {
                'NH2': ['N1'], # amine
                'CO': ['C1', 'C2'], # carbonyl
                'CC': ['C11', 'C12', 'C21', 'C22'], # alkene
            }
        )
        node_type_to_smiles = {
            'NH2': 'N([H])[H]',
            'CO': 'C=O',
            'CC': 'C=C',
        }
        node_port_to_atom_index = {
            'NH2': {'N1': 0},
            'CO': {'C1': 0, 'C2': 0},
            'CC': {'C11': 0, 'C12': 0, 'C21': 2, 'C22': 2},
        }
        graph.add_node('node1', 'NH2')
        graph.add_node('node2', 'CO')
        graph.add_edge('node1', 'N1', 'node2', 'C1')
        molecular_graph = graph.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
        print(molecular_graph)

