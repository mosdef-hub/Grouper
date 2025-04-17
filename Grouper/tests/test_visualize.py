import networkx as nx
import pytest

from Grouper import Group
from Grouper.libraries.Libraries import GroupExtension
from Grouper.tests.base_test import BaseTest
from Grouper.utils import convert_to_nx
from Grouper.visualization import visualize, visualize_group_extension


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
        nt = GroupExtension(Group("desc", smiles, hubs, True), "", smarts, None)
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
        # smartsList = ["[C](-[C])(-[C])(-[C])", "[NX3](-[C])(-[O-])(-[O-])", "[CX3]([C])(-[H])(-[H])","c1ccccc1", "c1(-[O](-[C]))ccccc1", "[C]", "bad_smarts"]
        # smilesList = ["CCCC", "N(=O)(O)", "[C]([H])([C])", "N", "C=O", "bad_smiles", "C"]
        # hubsList = [[0], [0], [0], [0,0], [0,0,1], [], []]
        with pytest.raises((ValueError, AttributeError)) as e:
            nt = GroupExtension(Group("desc", smiles, hubs, True), "", smarts, None)
            visualize_group_extension(nt)
            print(e, smarts, smiles, hubs)
