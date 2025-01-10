import pytest

import networkx as nx
import numpy as np
from Grouper import GroupGraph
from Grouper.tests.base_test import BaseTest
from Grouper.visualization import visualize
from Grouper.utils import convert_to_nx

class TestGroupGraph(BaseTest):
    @pytest.mark.parametrize("graph_fixture", ["empty_graph", "basic_graph", "single_node_graph"])
    def test_nx_visualize(self, request, graph_fixture):
        graph = request.getfixturevalue(graph_fixture)
        conversion = convert_to_nx(graph)
        if graph_fixture == "empty_graph":
            with pytest.raises(ValueError):
                pos = nx.spring_layout(conversion)
                fig = visualize(graph, pos)
        else:
            pos = nx.spring_layout(conversion)
            fig = visualize(graph, pos)
            assert fig