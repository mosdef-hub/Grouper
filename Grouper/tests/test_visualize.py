import networkx as nx
import pytest

from Grouper import Group, GroupGraph
from Grouper.libraries.Libraries import GroupExtension
from Grouper.tests.base_test import BaseTest
from Grouper.utils import convert_to_nx
from Grouper.visualization import visualize, visualize_group_extension
from Grouper.visualization.visualize_graph import nx_visualize


class TestGroupGraph(BaseTest):
    @pytest.mark.parametrize(
        "graph_fixture", ["empty_graph", "basic_graph", "single_node_graph"]
    )
    def test_visualize_graph(self, request, graph_fixture):
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

    @pytest.mark.parametrize(
        "smarts, smiles, hubs",
        list(
            zip(
                [
                    "[C](-[C])(-[C])(-[C])",
                    "[N+X3](-[C])(=[O])(-[O-])",
                    "[CX3]([C])(-[H])(-[H])",
                    "c1ccccc1",
                    "c1(-[O](-[C]))ccccc1",
                ],
                ["CC", "N(=O)(O)", "[C]([H])([H])", "C", "CO"],
                [[0], [0], [0], [0, 0], [0, 0, 1]],
            )
        ),
    )
    def test_visualize_group_extension(self, smarts, smiles, hubs):
        nt = GroupExtension(Group("desc", smiles, hubs, "SMARTS"), "", smarts, None)
        img = visualize_group_extension(nt)
        assert img

    @pytest.mark.parametrize(
        "smarts, smiles, hubs",
        list(
            zip(
                [
                    "[C](-[C])(-[C])(-[C])",
                    "[NX3](-[C])(-[O-])(-[O-])",
                    "[CX3]([C])(-[H])(-[H])",
                    "c1ccccc1",
                    "c1(-[O](-[C]))ccccc1",
                    "[C]",
                    "bad_smarts",
                ],
                [
                    "CCCC",
                    "N(=O)(O)",
                    "[C]([H])([C])([C])",
                    "N",
                    "C=O",
                    "bad_smiles",
                    "C",
                ],
                [[0], [0], [0], [0, 0], [0, 0, 1], [], []],
            )
        ),
    )
    def test_visualize_group_errors(self, smarts, smiles, hubs):
        with pytest.raises((ValueError, AttributeError)):
            nt = GroupExtension(Group("desc", smiles, hubs, "SMARTS"), "", smarts, None)
            visualize_group_extension(nt)
    
    def test_fail_parse_group_smarts(self):
        with pytest.raises(ValueError, match=r"Could not parse SMARTS\: CN\[IC\]"):
            nt = GroupExtension(Group("desc", "CN[IC]", [0], "SMARTS"), "doi", "CCCC", None)
            visualize_group_extension(nt)

    def test_nx_visualize(self):
        import matplotlib.pyplot as plt
        from Grouper.utils import convert_to_nx
        gG = GroupGraph()
        gG.add_node("node1", "CCC", [0,1,2])
        gG.add_node("node2", "OC", [0])
        gG.add_edge((0,0), (1,0), 2)
        nxGraph = convert_to_nx(gG)
        nxGraph = nxGraph.to_undirected()
        ax, fig = nx_visualize(nxGraph)
        assert isinstance(fig, plt.Figure)
        assert isinstance(ax, plt.Axes)

    def test_cytoscape_visualize(self):
        from Grouper.visualization.visualize_graph import visualize_cytoscape
        from dash import Dash

        gG = GroupGraph()
        gG.add_node("node1", "CCC", [0,1,2])
        gG.add_node("node2", "OC", [0])
        gG.add_edge((0,0), (1,0), 2)
        assert isinstance(visualize_cytoscape(gG), Dash)
