from molGrouper.group_graph import GroupGraph
from molGrouper.io import has_mbuild, has_torch
from molGrouper.post_process import substitute_chiral_smiles
from mbuild.tests.base_test import BaseTest
from group_selfies import Group
import networkx as nx
from pysmiles import write_smiles
import pytest

class TestGroupGraph(BaseTest):

    @pytest.fixture(autouse=True)
    def test_creation(self):
        # Define node types with ports
        node_types = {
            'type1': ['port1', 'port2'],
            'type2': ['port3', 'port4'],
            'r1' : ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
            'd6' : ['C1', 'C2'], # alkane
            'c' : ['C1', 'C2', 'C3', 'C4'],
            'r3': ['C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22'],
            's1': ['C1'],
            'sulfonate_ester': ['C1', 'O1'],
        }
        self.node_types = node_types

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

    def test_add_edge_with_invalid_nodes(self):
        with pytest.raises(KeyError):
            self.graph.add_edge('node1', 'port1', 'node2', 'port1')
    
    def test_add_edge_with_invalid_ports(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type1')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node1', 'port1', 'node2', 'port3')
            
    def test_add_edge_with_too_many_ports(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type1')
        self.graph.add_node('node3', 'type2')
        self.graph.add_edge('node1', 'port1', 'node2', 'port1')
        self.graph.add_edge('node1', 'port2', 'node2', 'port2')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node3', 'port3', 'node2', 'port2')

    def test_add_edge_with_same_ports(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type1')
        self.graph.add_node('node3', 'type2')
        self.graph.add_edge('node1', 'port1', 'node2', 'port1')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node1', 'port1', 'node2', 'port1')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node1', 'port1', 'node2', 'port2')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node1', 'port2', 'node2', 'port1')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node1', 'port2', 'node2', 'port1')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node3', 'port1', 'node2', 'port1')
        with pytest.raises(AttributeError):
            self.graph.add_edge('node1', 'port1', 'node3', 'port1')

    def test_equal(self):
        graph1 = GroupGraph(self.node_types)
        graph2 = GroupGraph(self.node_types)
        assert graph1 == graph2

        graph1.add_node('node1', 'type1')
        assert graph1 != graph2

        graph2.add_node('node1', 'type1')
        assert graph1 == graph2

        graph1.add_node('node2', 'type2')
        assert graph1 != graph2

        graph2.add_node('node2', 'type2')
        assert graph1 == graph2

        graph1.add_edge('node1', 'port1', 'node2', 'port3')
        assert graph1 != graph2

        graph2.add_edge('node1', 'port1', 'node2', 'port3')
        assert graph1 == graph2
        
    def test_in(self):
        graph1 = GroupGraph(self.node_types)
        graph2 = GroupGraph(self.node_types)

        assert graph1 in [graph1]

        graph1.add_node('node1', 'type1')
        assert graph1 not in [graph2]

        graph2.add_node('node1', 'type1')
        assert graph1 in [graph2]

        graph1.add_node('node2', 'type2')
        assert graph1 not in [graph2]

        graph2.add_node('node2', 'type2')
        assert graph1 in [graph2]

        graph1.add_edge('node1', 'port1', 'node2', 'port3')
        assert graph1 not in [graph2]

        graph2.add_edge('node1', 'port1', 'node2', 'port3')
        assert graph1 in [graph2]

        assert graph1 in [graph1, graph2]

    def test_rings(self):
        self.graph.add_node('n0', 'r1')
        self.graph.add_node('n1', 'r1')
        self.graph.add_node('n2', 'r1')

        self.graph.add_edge('n0', 'C1', 'n2', 'C1')
        with pytest.raises(AttributeError):
            self.graph.add_edge('n1', 'C1', 'n2', 'C1')

    def test_to_molecular_graph(self):
        node_type_to_smiles = {
            'r1' : 'c1ccccc1',
            'd6': 'C',
            'c': 'C',
            'r3': 'C1CCCCC1'
        }
        node_port_to_atom_index = {
            'r1' : {'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4, 'C6': 5},
            'd6' : {'C1': 0, 'C2': 0},
            'c': {'C1': 0, 'C2': 0, 'C3': 0, 'C4': 0},
            'r3': {'C11': 0, 'C12': 0, 'C13': 1, 'C14': 1, 'C15': 2, 'C16': 2, 'C17': 3, 'C18': 3, 'C19': 4, 'C20': 4, 'C21': 5, 'C22': 5}
        }
        node_port_to_atom_index = {
            'r1' : {'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4, 'C6': 5},
            'd6' : {'C1': 0, 'C2': 0},
            'c': {'C1': 0, 'C2': 0, 'C3': 0, 'C4': 0},
            'r3': {'C11': 0, 'C12': 0, 'C13': 1, 'C14': 1, 'C15': 2, 'C16': 2, 'C17': 3, 'C18': 3, 'C19': 4, 'C20': 4, 'C21': 5, 'C22': 5}
        }
        self.graph.add_node('n0', 'r1')
        self.graph.add_node('n1', 'r1')
        self.graph.add_node('n2', 'r1')
        self.graph.add_edge('n0', 'C2', 'n2', 'C2')
        self.graph.add_edge('n1', 'C1', 'n2', 'C1')
        molecular_graph = self.graph.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
        assert set(molecular_graph.edges) == set([(0, 1), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5), (1, 13), (6, 7), (6, 11), (7, 8), (8, 9), (9, 10), (10, 11), (6, 12), (12, 13), (12, 17), (13, 14), (14, 15), (15, 16), (16, 17)])

        self.graph = GroupGraph(
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
            'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
        }
        self.graph.add_node('node1', 'NH2')
        self.graph.add_node('node2', 'CO')
        self.graph.add_edge('node1', 'N1', 'node2', 'C1')
        molecular_graph = self.graph.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)

        self.graph = GroupGraph( self.node_types )
        node_type_to_smiles = {
            's1': 'C',
            'r1': 'c1ccccc1',
            'c': 'C',
        }
        node_port_to_atom_index = {
            's1': {'C1': 0},
            'r1': {'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4, 'C6': 5},
            'c': {'C1': 0, 'C2': 0, 'C3': 0, 'C4': 0},
        }
        self.graph.add_node('n0', 'c')
        self.graph.add_node('n1', 's1')
        self.graph.add_node('n2', 'r1')
        self.graph.add_edge('n0', 'C4', 'n2', 'C3')
        self.graph.add_edge('n1', 'C1', 'n2', 'C1')
        molecular_graph = self.graph.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
        assert set(molecular_graph.edges) == set([(0, 4), (1, 2), (2, 7), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)])

    def test_n_free_ports(self):
        self.graph.add_node('node1', 'type1')
        self.graph.add_node('node2', 'type2')
        assert self.graph.n_free_ports('node1') ==  2

        # Connect a edge and recheck
        self.graph.add_edge('node1', 'port1', 'node2', 'port3')
        assert self.graph.n_free_ports('node1') == 1
        assert self.graph.n_free_ports('node2') == 1

    def test_chiral_smiles_conversion(self):
        # Define node types with ports
        node_types = {
            'CC': ['C11', 'C12', 'C21', 'C22',],
            'OH': ['O1'],
        }
        node_types_to_smiles = {
            'CC': 'C=C',
            'OH': 'O',
        }
        node_port_to_atom_index = {
            'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
            'OH': {'O1': 0},
        }
        node_type_to_chiral_subs = {
            'CC': {'[C]=[C]' : ['/[C]=[C]/', '/[C]=[C]\\']},
            'OH': {'O': []} # no chiral subs
        }

        # cis
        graph = GroupGraph(node_types)
        graph.add_node('n0', 'OH')
        graph.add_node('n1', 'CC')
        graph.add_node('n2', 'OH')
        graph.add_edge('n0', 'O1', 'n1', 'C11')
        graph.add_edge('n2', 'O1', 'n1', 'C22')
        mG = graph.to_molecular_graph(node_types_to_smiles, node_port_to_atom_index)
        smiles = write_smiles(mG)

        chiral_smiles = substitute_chiral_smiles(smiles, '[C]=[C]', node_type_to_chiral_subs['CC']['[C]=[C]'])
        #cis
        assert '[O]/[C]=[C]\\[O]' in chiral_smiles
        #trans 
        assert '[O]/[C]=[C]/[O]' in chiral_smiles

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
