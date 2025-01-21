import networkx as nx
import pytest

from Grouper import Node
from Grouper.libraries.Libraries import NodeTrace
from Grouper.tests.base_test import BaseTest
from Grouper.utils import convert_to_nx
from Grouper.visualization import visualize, visualize_node_trace


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
    def test_visualize_node(self, smarts, smiles, hubs):
        nt = NodeTrace(Node("desc", smiles, hubs), "", smarts, None)
        img = visualize_node_trace(nt)
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
    def test_visualize_node_errors(self, smarts, smiles, hubs):
        # smartsList = ["[C](-[C])(-[C])(-[C])", "[NX3](-[C])(-[O-])(-[O-])", "[CX3]([C])(-[H])(-[H])","c1ccccc1", "c1(-[O](-[C]))ccccc1", "[C]", "bad_smarts"]
        # smilesList = ["CCCC", "N(=O)(O)", "[C]([H])([C])", "N", "C=O", "bad_smiles", "C"]
        # hubsList = [[0], [0], [0], [0,0], [0,0,1], [], []]
        nt = NodeTrace(Node("desc", smiles, hubs), "", smarts, None)
        with pytest.raises((ValueError, AttributeError)) as e:
            visualize_node_trace(nt)
            print(e, smarts, smiles, hubs)
