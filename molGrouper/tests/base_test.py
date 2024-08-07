import pytest
from copy import deepcopy

from molGrouper import GroupGraph
from group_selfies import Group
import networkx as nx

class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture(autouse=True)
    def empty_graph(self):
        return GroupGraph()

    @pytest.fixture(autouse=True)
    def basic_graph(self):
        group_graph = GroupGraph({"C": ["1", "2"]})
        group_graph.add_node("node0", "C")
        group_graph.add_node("node1", "C")
        group_graph.add_edge(("node0", "1"), ("node1", "2"))
        return group_graph

    @pytest.fixture(autouse=True)
    def single_node_graph(self):
        group_graph = GroupGraph({"C":["1", "2"]})
        group_graph.add_node("node0", "C")
        return group_graph

    @pytest.fixture(autouse=True)
    def long_graph(self):
        group_graph = GroupGraph()
        groups = [
            Group('c2', "C(*1)(*1)(*1)(*1)"),
        ]
        smiles = "CCCCCCCCCCC"
        return group_graph.group_graph_from_smiles(smiles, groups)

    @pytest.fixture(autouse=True)
    def branchy_graph(self):
        group_graph = GroupGraph()
        groups = [
            Group('c4', "C(*1)(*1)(*1)(*1)"),
        ]
        smiles = "C(CC(C))(C(C)C)CC"
        return group_graph.group_graph_from_smiles(smiles, groups)

    @pytest.fixture(autouse=True)
    def two_molecule_graph(self, branchy_graph):
        copy_graph = deepcopy(branchy_graph)
        union = nx.operators.union(
            copy_graph, branchy_graph, rename=("1", "2")
        )
        return union


