import pytest

import networkx as nx
import numpy as np
from molGrouper import GroupGraph
from molGrouper.tests.base_test import BaseTest
from group_selfies import Group

class TestGroupGraph(BaseTest):
    def test_nx_visualize_large(self, two_molecule_graph):
        fig = two_molecule_graph.visualize()
        assert fig

    def test_nx_visualize_small(self, single_node_graph):
        fig = single_node_graph.visualize()
        assert fig