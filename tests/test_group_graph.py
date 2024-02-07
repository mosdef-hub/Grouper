from mbuild.tests.base_test import BaseTest
from molGrouper.group_graph import GroupGraph
from group_selfies import Group
import mbuild as mb
import torch
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

    def test_str_representation(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type2')
        self.graph.add_edge('node1', 'port1', 'node2', 'port3')

        expected_str = "Nodes: [('node1', ['port1', 'port2']), ('node2', ['port3', 'port4'])] \nEdges: [('node1', 'node2', [['node1.port1', 'node2.port3']])]"
        assert str(self.graph) == expected_str

    def test_compound_to_group_graph(self):
        mol = mb.load('CCCCCCCC', smiles=True) # octane molecule
        groups = [Group('c3', 'C([H])([H])([H])(*1)'), Group('c2', 'C([H])([H])(*1)(*1)')]
        graph = GroupGraph()
        group_graph = graph.from_compound(mol, groups)
        print(group_graph)
        

    def test_group_graph_to_vector(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type2')
        self.graph.add_edge('node1', 'port1', 'node2', 'port3')
        vector_form = self.graph.to_vector()
        assert vector_form == [1, 1]

    def test_group_graph_to_pyG(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type2')
        self.graph.add_edge('node1', 'port1', 'node2', 'port3')
        group_featurizer = lambda node: torch.tensor([1, 0])
        
        data = self.graph.to_data(group_featurizer, max_n_attachments=2)
        assert torch.equal(data.x, torch.tensor([ [1,0], [1,0] ], dtype=torch.float32)) # node features should just be identity
        assert torch.equal(data.edge_index, torch.tensor([ [0], [1] ], dtype=torch.float32)) # graph is directed, node1 -> node2
        assert torch.equal(data.edge_attr, torch.tensor([ [1,0,1,0] ], dtype=torch.float32)) # edge features are one-hot encoded port

if __name__ == "__main__":
    test = TestGroupGraph()

    test.test_add_node()
    test.test_add_edge()
    test.test_make_undirected()
    test.test_n_free_ports()
    test.test_str_representation()
    test.test_Compound_to_group_graph()