from Grouper import GroupGraph, Group, fragment
from Grouper.tests.base_test import BaseTest

# from pysmiles import write_smiles
# import pytest


class TestGroupGraph(BaseTest):
    def test_fragment(self):
        node_defs = set()
        node_defs.add(Group("hydroxyl", "O", [0]))
        node_defs.add(Group("alkene", "C=C", [0, 0, 1, 1]))

        truth = GroupGraph()
        truth.add_node("hydroxyl", "O", [0])
        truth.add_node("alkene", "C=C", [0, 0, 1, 1])
        truth.add_edge((0, 0), (1, 0))

        smiles = "C=CO"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_fragment_2(self):
        node_defs = set()
        node_defs.add(Group("amine", "N", [0, 0]))
        node_defs.add(Group("alkene", "C=C", [0, 0, 1, 1]))

        truth = GroupGraph()
        truth.add_node("amine", "N", [0, 0])
        truth.add_node("alkene", "C=C", [0, 0, 1, 1])
        truth.add_edge((0, 0), (1, 0))

        smiles = "C=CN"

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
        node_defs.add(Group("amine", "N", [0, 0]))
        node_defs.add(Group("alkene", "C=C", [0, 0, 1, 1]))
        node_defs.add(Group('oxygen', 'O', [0, 0]))

        truth = GroupGraph()
        truth.add_node("oxygen", "O", [0, 0])
        truth.add_node("alkene", "C=C", [0, 0, 1, 1])
        truth.add_node("amine", "N", [0])
        truth.add_edge((0, 0), (1, 0))
        truth.add_edge((0, 1), (2, 0))

        smiles = "C=CON"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_triple_node(self):
        node_defs = set()
        node_defs.add(Group("oxyl", "O", [0, 0]))  # oxyl group
        node_defs.add(Group("ester", "C(=O)(O)", [0, 2]))  # Ester group
        node_defs.add(Group("amine", "N", [0, 0, 0]))  # Amine group
        node_defs.add(Group("alkene_secondary_amine", "CNC", [0, 0]))
        node_defs.add(Group("alkene", "C", [0, 0, 0]))

        truth = GroupGraph()
        truth.add_node("ester", "C(=O)O", [0, 2])
        truth.add_node("oxyl", "O", [0, 0])
        truth.add_node("ester", "C(=O)(O)", [0, 2])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))

        out = fragment("O=COOC(=O)O", node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node("ester", "C(=O)O", [0, 2])
        truth.add_node("alkene_secondary_amine", "C(N(C))", [0, 0])
        truth.add_node("amine", "N", [0, 0, 0])
        truth.add_edge((0, 1), (2, 0))
        truth.add_edge((1, 1), (2, 1))

        out = fragment("CNCNOC=O", node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node("alkene", "C", [0, 0, 0])
        truth.add_node("amine", "N", [0, 0, 0])
        truth.add_node("oxyl", "O", [0, 0])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))

        out = fragment("CNO", node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node("oxyl", "O", [0, 0])
        truth.add_node("oxyl", "O", [0, 0])
        truth.add_node("ester", "C(=O)O", [0, 2])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))

        out = fragment("O=C(O)OO", node_defs)

        assert out == truth

    def test_nodes_made_of_other_nodes(self):
        node_defs = set()
        node_defs.add(Group("oxyl", "O", [0, 0]))  # oxyl group
        node_defs.add(Group("ester", "C(=O)O", [0, 2]))  # Ester group
        node_defs.add(Group("amine", "N", [0, 0, 0]))  # Amine group
        node_defs.add(Group(
            "alkene_secondary_amine", "C(N(C))", [0, 0]
        ))  # can be made of amine and alkene
        node_defs.add(Group("alkene", "C", [0, 0, 0]))

        truth = GroupGraph()
        truth.add_node("alkene_secondary_amine", "C(N(C))", [0, 0])
        truth.add_node("alkene", "C", [0, 0, 0])
        truth.add_node("amine", "N", [0, 0, 0])
        truth.add_node("alkene", "C", [0, 0, 0])
        truth.add_node("alkene", "C", [0, 0, 0])
        truth.add_node("alkene", "C", [0, 0, 0])
        truth.add_node("alkene", "C", [0, 0, 0])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))
        truth.add_edge((2, 1), (3, 0))
        truth.add_edge((3, 1), (4, 0))
        truth.add_edge((3, 2), (5, 0))
        truth.add_edge((5, 1), (6, 0))

        out = fragment("CCC(C)NCCNC", node_defs)

        assert out == truth


    def test_other_edge_cases(self):
        node1 = Group("1methyl", "C", [0,0,0])
        node2 = Group("2methyl", "C", [0,0])

        truth = GroupGraph()
        truth.add_node("1methyl", "C", [0])
        truth.add_node("1methyl", "C", [0])
        truth.add_node('2methyl', 'C', [0, 0])
        truth.add_edge((0, 0), (2, 0))
        truth.add_edge((1, 0), (2, 1))

        node_defs = {node2, node1} # {node2}
        out = fragment("CCC", node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node("carbon", "C", [0,0,0,0])
        truth.add_node("ether", "CO", [0, 1])
        truth.add_node("ether", "CO", [0, 1])
        truth.add_edge((0, 0), (2, 0))
        truth.add_edge((2, 1), (1, 0))

        node1 = Group("4methyl", "C", [0,0,0,0])
        node2 = Group("methanol", "CO", [0]) # methanol
        node3 = Group("ether", "CO", [0,1]) # ether

        node_defs = {node3, node2, node1} # {node2}
        graph = fragment("CCOCO", node_defs)

        assert graph == truth

    def test_larger_mols(self):
        node_defs = set()
        node_defs.add(Group('hydroxyl', 'O', [0]))
        node_defs.add(Group('keto', 'O', [0,0]))
        node_defs.add(Group('ester', 'C(=O)O', [0,2]))
        node_defs.add(Group('methyl', 'C', [0,0,0]))
        node_defs.add(Group('t2', 'N', [0,0,0]))
        node_defs.add(Group('secondary_amine', 'CNC', [0,0, 1]))

        truth = GroupGraph()
        truth.add_node("secondary_amine", "CNC", [0, 0])
        truth.add_node("t2", "N", [0, 0, 0])
        truth.add_node("ester", "C(=O)O", [0, 2])
        truth.add_edge((0, 0), (1, 0))
        truth.add_edge((1, 0), (2, 0))

        graph = fragment("CNCNOC=O", node_defs)
