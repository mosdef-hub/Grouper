from genGrouper import GroupGraph
from genGrouper.io import has_mbuild, has_torch
from genGrouper.post_process import substitute_chiral_smiles
from genGrouper.tests.base_test import BaseTest
from group_selfies import Group
import networkx as nx
from pysmiles import write_smiles
import pytest
import numpy as np

class TestGroupGraph(BaseTest):
    def test_add_node(self):
        # Basic node addition
        graph = GroupGraph()
        graph.add_node('type1', '', [0,1], [0,0])
        assert 0 in [n.id for n in graph.nodes.values()]
        assert set(n.type for n in graph.nodes.values()) == set(['type1'])
        assert set(n.smiles for n in graph.nodes.values()) == set([''])
        assert [n.ports for n in graph.nodes.values()] == [[0,1]]

        # Adding a node with different type and smiles
        graph.add_node('', 'C', [0,1], [0,0])
        assert len(graph.nodes) == 2
        assert set(n.type for n in graph.nodes.values()) == set(['type1', 'C'])
        assert set(n.smiles for n in graph.nodes.values()) == set(['', 'C'])
        assert [n.ports for n in graph.nodes.values()] == [[0,1], [0,1]]

        # Adding a node with only a type
        graph.add_node('type1')
        assert len(graph.nodes) == 3
        assert set(n.type for n in graph.nodes.values()) == set(['type1', 'C', 'type1'])
        assert set(n.smiles for n in graph.nodes.values()) == set(['', 'C', ''])
        assert [n.ports for n in graph.nodes.values()] == [[0,1], [0,1], [0,1]]

    def test_add_edge(self):
        graph = GroupGraph()
        graph.add_node('type1', '', [0,1], [0,0])
        graph.add_node('', 'C', [0,1], [0,0])
        graph.add_node('type1', '', [0,1], [0,0])

        graph.add_edge((0, 0), (1, 0))
        assert (0,0,1,0) in graph.edges   

        graph.add_edge((2, 1), (1, 1))
        assert (2,1,1,1) in graph.edges

    def test_add_edge_with_invalid_nodes(self):
        graph = GroupGraph()
        graph.add_node('node1', '', [0,1], [0,0])
        graph.add_node('node2', '', [0], [0])
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (2, 1))

    def test_add_edge_with_invalid_ports(self):
        graph = GroupGraph()
        graph.add_node('node1', '', [0,1], [0,0])
        graph.add_node('node2', '', [0,1], [0,0])
        with pytest.raises(ValueError):
            graph.add_edge((0, 2), (1, 1))

    def test_add_edge_with_occupied_port(self):
        graph = GroupGraph()
        graph.add_node('node1', '', [0,1], [0,0])
        graph.add_node('node2', '', [0,1], [0,0])
        graph.add_edge((0, 1), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (1, 0))

    def test_add_edge_with_same_ports(self):
        graph = GroupGraph()
        graph.add_node('node1', '', [0,1], [0,0])
        graph.add_node('node2', '', [0,1], [0,0])
        graph.add_node('node3', '', [0,1], [0,0])
        graph.add_edge((0, 1), (1, 1))
        graph.add_edge((1, 0), (2, 1))
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (1, 0))
        with pytest.raises(ValueError):
            graph.add_edge((0, 0), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((1, 1), (2, 0))
        with pytest.raises(ValueError):
            graph.add_edge((2,0), (1,1))
        with pytest.raises(ValueError):
            graph.add_edge((1,1), (2,1))

    def test_equal(self):
        graph1 = GroupGraph()
        graph2 = GroupGraph()
        assert graph1 == graph2
        graph1.add_node('node1', 'type1', [0,1], [0,0])
        assert graph1 != graph2
        graph2.add_node('node1', 'type1', [0,1], [0,0])
        assert graph1 == graph2
        graph1.add_node('node2', 'type2', [0,1], [0,0])
        assert graph1 != graph2
        graph2.add_node('node2', 'type2', [0,1], [0,0])
        assert graph1 == graph2
        graph1.add_edge((0, 0), (1, 0))
        assert graph1 != graph2
        graph2.add_edge((0, 0), (1, 0))
        assert graph1 == graph2

    def test_in(self):
        graph1 = GroupGraph()
        graph2 = GroupGraph()

        assert graph1 in [graph1]

        graph1.add_node('node1', 'type1', [0,1], [0,0])
        assert graph1 not in [graph2]

        graph2.add_node('node1', 'type1', [0,1], [0,0])
        assert graph1 in [graph2]

        graph1.add_node('node2', 'type2', [0,1], [0,0])
        assert graph1 not in [graph2]

        graph2.add_node('node2', 'type2', [0,1], [0,0])
        assert graph1 in [graph2]

        graph1.add_edge((0, 0), (1, 0))
        assert graph1 not in [graph2]

        graph2.add_edge((0, 0), (1, 0))
        assert graph1 in [graph2]

        assert graph1 in [graph1, graph2]

    def test_rings(self):
        graph = GroupGraph()
        graph.add_node('ring', 'c1ccccc1', [0,1,2,3,4,5], [0,1,2,3,4,5])
        graph.add_node('ring')
        graph.add_node('ring')

        graph.add_edge((0, 0), (1, 0))
        with pytest.raises(ValueError):
            graph.add_edge((0, 0), (2, 0))
    @pytest.mark.parametrize("graph_fixture", ["empty_graph", "basic_graph", "single_node_graph"])
    def test_to_atomic_graph(self, request, graph_fixture):
        # Access the graph using request.getfixturevalue
        graph = request.getfixturevalue(graph_fixture)

        # Example assertions for different graphs
        if graph_fixture == "empty_graph":
            assert graph.to_atomic_graph() == None
        elif graph_fixture == "basic_graph":
            atomic_graph = graph.to_atomic_graph()
            assert set([n.type for n in atomic_graph.nodes])
        # elif graph_fixture == "single_node_graph":
        # molecular_graph = graph.to_atomic_graph()
        # molecular_graph = graph.to_atomic_graph(node_type_to_smiles, node_port_to_atom_index)

    def test_n_free_ports(self):
        graph = GroupGraph()
        graph.add_node('node1', 'type1', [0,1], [0,0])
        graph.add_node('node2', 'type2', [0,1,2], [0,0,1])
        assert graph.n_free_ports(0) ==  2

        # Connect a edge and recheck
        graph.add_edge((0, 0), (1, 0))
        assert graph.n_free_ports(0) == 1
        assert graph.n_free_ports(1) == 2

    @pytest.mark.parametrize("graph_fixture", ["empty_graph", "basic_graph", "single_node_graph", "single_edge_graph", "five_member_ring_graph"])
    def test_to_smiles(self, request, graph_fixture):
        graph = request.getfixturevalue(graph_fixture)
        if graph_fixture == "empty_graph":
            assert graph.to_smiles() == ''
        if graph_fixture == "basic_graph":
            assert graph.to_smiles() == 'CN'
        if graph_fixture == "single_node_graph":
            assert graph.to_smiles() == 'C'
        if graph_fixture == "single_edge_graph":
            assert graph.to_smiles() == 'CC'
        if graph_fixture == "five_member_ring_graph":
            assert graph.to_smiles() == 'C1CCCC1'

        graph = GroupGraph()
        print(graph)
        graph.add_node('benzene', 'c1ccccc1', [0,1,2,3,4,5], [0,1,2,3,4,5])
        assert graph.to_smiles() == 'C1=CC=CC=C1' or graph.to_smiles() == 'c1ccccc1'
        graph.add_node('benzene', 'C1=CC=CC=C1', [0,1,2,3,4,5], [0,1,2,3,4,5])
        graph.add_node('benzene', 'C1=CC=CC=C1', [0,1,2,3,4,5], [0,1,2,3,4,5])
        graph.add_node('benzene', 'C1=CC=CC=C1', [0,1,2,3,4,5], [0,1,2,3,4,5])
        graph.add_edge((0, 2), (3, 5))
        graph.add_edge((1, 0), (3, 0))
        graph.add_edge((2, 1), (3, 3))
        assert graph.to_smiles() == 'c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1'
        
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
        graph.add_edge(('n0', 'O1'), ('n1', 'C11'))
        graph.add_edge(('n2', 'O1'), ('n1', 'C22'))
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
        def get_ground_truth(graph):
            hist = {node: 0 for node in graph.node_types}
            for k, v in graph.nodes.items():
                hist[v.type] += 1
            return hist

        graph = GroupGraph()
        graph.add_node('node1', '', [0,1], [0,0])
        graph.add_node('node2', '', [0,1], [0,0])

        group_vector = graph.to_vector()
        true_hist = get_ground_truth(graph)

        assert len(group_vector) == len(true_hist)
        for ntype in group_vector:
            assert group_vector[ntype] == true_hist[ntype]

    @pytest.mark.skipif(not has_torch, reason="torch package not installed")
    def test_group_graph_to_pyG(self, basic_graph):
        import torch
        group_featurizer = lambda node: torch.tensor([1, 0])

        data = basic_graph.to_PyG_Data(group_featurizer)
        assert torch.equal(data.x, torch.tensor([ [1,0], [1,0] ], dtype=torch.float32)) # node features should just be identity
        assert torch.equal(data.edge_index, torch.tensor([ [0], [1] ], dtype=torch.float32)) # graph is directed, node1 -> node2
        assert torch.equal(data.edge_attr, torch.tensor([ [1.,0.,0.,1.] ], dtype=torch.float32)) # edge features are one-hot encoded port


