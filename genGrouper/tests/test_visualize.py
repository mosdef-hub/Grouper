import pytest

import networkx as nx
import numpy as np
from genGrouper import GroupGraph
from genGrouper.tests.base_test import BaseTest
from genGrouper.visualization import nx_visualize
from genGrouper.utils import convert_to_nx

class TestGroupGraph(BaseTest):
    @pytest.mark.parametrize("graph_fixture", ["empty_graph", "basic_graph", "single_node_graph"])
    def test_nx_visualize(self, request, graph_fixture):
        graph = request.getfixturevalue(graph_fixture)
        conversion = convert_to_nx(graph)
        if graph_fixture == "empty_graph":
            with pytest.raises(ValueError):
                fig = nx_visualize(conversion)
        else:
            fig = nx_visualize(conversion)
            assert fig