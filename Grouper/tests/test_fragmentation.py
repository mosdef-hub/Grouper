from Grouper import GroupGraph, Node, fragment
from Grouper.tests.base_test import BaseTest

# from pysmiles import write_smiles
# import pytest


class TestGroupGraph(BaseTest):
    def test_fragment(self):
        node_defs = set()
        node_defs.add(Node("hydroxyl", "[OX2H]", [0]))
        node_defs.add(Node("alkene", "C=C", [0, 0, 1, 1]))

        truth = GroupGraph()
        truth.add_node("hydroxyl", "[OH]", [0])
        truth.add_node("alkene", "C=C", [0, 0, 1, 1])
        truth.add_edge((0, 0), (1, 0))

        smiles = "C=CO"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_fragment_2(self):
        node_defs = set()
        node_defs.add(Node("amine", "[NH2]", [0, 0]))
        node_defs.add(Node("alkene", "C=C", [0, 0, 1, 1]))

        truth = GroupGraph()
        truth.add_node("amine", "[NH2]", [0, 0])
        truth.add_node("alkene", "C=C", [0, 0, 1, 1])
        truth.add_edge((0, 0), (1, 0))

        smiles = "C=C[NH2]"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_empty_smiles(self):
        node_defs = set()
        smiles = ""

        out = fragment(smiles, node_defs)

        assert out == GroupGraph()

    def test_invalid_smiles(self):
        node_defs = set()
        smiles = "C==C"  # Invalid SMILES

        try:
            fragment(smiles, node_defs)
        except Exception as e:
            assert isinstance(
                e, ValueError
            )  # or the appropriate exception class for invalid SMILES

    def test_multiple_bonds(self):
        node_defs = set()
        node_defs.add(Node("amine", "[NX3]", [0, 0]))
        node_defs.add(Node("alkene", "C=C", [0, 0, 1, 1]))
        node_defs.add(Node('oxygen', '[O]', [0, 0]))

        truth = GroupGraph()
        truth.add_node("oxygen", "[O]", [0, 0])
        truth.add_node("alkene", "C=C", [0, 0, 1, 1])
        truth.add_node("amine", "[NX3]", [0])
        truth.add_edge((0, 0), (1, 0))
        truth.add_edge((0, 1), (2, 0))

        smiles = "C=CO[NH2]"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_triple_node(self):
        node_defs = set()
        node_defs.add(Node("oxyl", "[O]", [0, 0]))  # oxyl group
        node_defs.add(Node("ester", "[CX3](=O)[OX2H0]", [0, 2]))  # Ester group
        node_defs.add(Node("amine", "[NX3]", [0, 0, 0]))  # Amine group
        node_defs.add(Node("alkene_secondary_amine", "[CX4]([NX3]([CX1]))", [0, 0]))
        node_defs.add(Node("alkene", "[CX4]", [0, 0, 0]))

        truth = GroupGraph()
        truth.add_node("ester", "[CX3](=O)[OX2H0]", [0, 2])
        truth.add_node("oxyl", "[O]", [0, 0])
        truth.add_node("ester", "[CX3](=O)[OX2H0]", [0, 2])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))

        out = fragment("O=COOC(=O)O", node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node("ester", "[CX3](=O)[OX2H0]", [0, 2])
        truth.add_node("alkene_secondary_amine", "[CX4]([NX3]([CX1]))", [0, 0])
        truth.add_node("amine", "N", [0, 0, 0])
        truth.add_edge((0, 1), (2, 0))
        truth.add_edge((1, 1), (2, 1))

        out = fragment("CNCNOC=O", node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node("alkene", "[CX4]", [0, 0, 0])
        truth.add_node("amine", "[NX3]", [0, 0, 0])
        truth.add_node("oxyl", "[O]", [0, 0])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))

        out = fragment("CNO", node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node("oxyl", "[O]", [0, 0])
        truth.add_node("oxyl", "[O]", [0, 0])
        truth.add_node("ester", "[CX3](=O)[OX2H0]", [0, 2])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))

        out = fragment("O=C(O)OO", node_defs)

        assert out == truth

    def test_nodes_made_of_other_nodes(self):
        node_defs = set()
        node_defs.add(Node("oxyl", "[OX2]", [0, 0]))  # oxyl group
        node_defs.add(Node("ester", "[CX3](=O)[OX2H0]", [0, 2]))  # Ester group
        node_defs.add(Node("amine", "[NX3]", [0, 0, 0]))  # Amine group
        node_defs.add(Node(
            "alkene_secondary_amine", "[CX4]([NX3]([CX1]))", [0, 0]
        ))  # can be made of amine and alkene
        node_defs.add(Node("alkene", "[CX4]", [0, 0, 0]))

        truth = GroupGraph()
        truth.add_node("alkene_secondary_amine", "[CX4]([NX3]([CX1]))", [0, 0])
        truth.add_node("alkene", "[CX4]", [0, 0, 0])
        truth.add_node("amine", "[NX3]", [0, 0, 0])
        truth.add_node("alkene", "[CX4]", [0, 0, 0])
        truth.add_node("alkene", "[CX4]", [0, 0, 0])
        truth.add_node("alkene", "[CX4]", [0, 0, 0])
        truth.add_node("alkene", "[CX4]", [0, 0, 0])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))
        truth.add_edge((2, 1), (3, 0))
        truth.add_edge((3, 1), (4, 0))
        truth.add_edge((3, 2), (5, 0))
        truth.add_edge((5, 1), (6, 0))

        out = fragment("CCC(C)NCCNC", node_defs)

        assert out == truth
