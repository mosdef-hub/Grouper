from Grouper.tests.base_test import BaseTest
from Grouper import GroupGraph, fragment, Node
# from pysmiles import write_smiles
# import pytest

class TestGroupGraph(BaseTest):
    def test_fragment(self):
        node_defs = {}
        node_defs['[OX2H]'] = Node(0, 'hydroxyl', 'O', [0])
        node_defs['C=C'] = Node(1, 'alkene', 'C=C', [0,0,1,1])

        truth = GroupGraph()
        truth.add_node('hydroxyl', '[OH]', [0])
        truth.add_node('alkene', 'C=C', [0,0,1,1])
        truth.add_edge((0,0),(1,0))

        smiles = "C=CO"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_fragment_2(self):
        node_defs = {}
        node_defs['[NH2]'] = Node(0, 'amine', 'N', [0,0])
        node_defs['C=C'] = Node(1, 'alkene', 'C=C', [0,0,1,1])

        truth = GroupGraph()
        truth.add_node('amine', '[NH2]', [0,0])
        truth.add_node('alkene', 'C=C', [0,0,1,1])
        truth.add_edge((0,0),(1,0))

        smiles = "C=C[NH2]"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_empty_smiles(self):
        node_defs = {}
        smiles = ""
        
        out = fragment(smiles, node_defs)
        
        assert out == GroupGraph()

    def test_invalid_smiles(self):
        node_defs = {}
        smiles = "C==C"  # Invalid SMILES
        
        try:
            out = fragment(smiles, node_defs)
        except Exception as e:
            assert isinstance(e, ValueError)  # or the appropriate exception class for invalid SMILES

    def test_multiple_bonds(self):
        node_defs = {}
        node_defs['[OX2]'] = Node(0, 'oxygen', 'O', [0,0])
        node_defs['C=C'] = Node(1, 'alkene', 'C=C', [0,0,1,1])
        node_defs['[NH2]'] = Node(2, 'amine', 'N', [0])

        truth = GroupGraph()
        truth.add_node('oxygen', 'O', [0,0])
        truth.add_node('alkene', 'C=C', [0,0,1,1])
        truth.add_node('amine', 'N', [0])
        truth.add_edge((0,0),(1,0))
        truth.add_edge((0,1),(2,0))

        smiles = "C=CO[NH2]"

        out = fragment(smiles, node_defs)

        assert out == truth

    def test_triple_node(self):
        node_defs = {}
        node_defs['[OX2]'] = Node(0, 'oxyl', 'O', [0,0])  # oxyl group
        node_defs['[CX3](=O)[OX2H0]'] = Node(0, 'ester', 'C(=O)O', [0,2])  # Ester group
        node_defs['[NX3]'] = Node(0, 'amine', 'N', [0,0,0])  # Amine group
        node_defs['[CX4]([NX3]([CX4]))'] = Node(0, 'alkene_secondary_amine', 'CNC', [0,0])
        

        truth = GroupGraph()
        truth.add_node('ester', 'C(=O)O', [0, 2])
        truth.add_node('oxyl', 'O', [0,0])
        truth.add_node('ester', 'C(=O)O', [0, 2])
        truth.add_edge((0, 1), (1, 0))
        truth.add_edge((1, 1), (2, 0))

        out = fragment('O=COOC(=O)O', node_defs)

        assert out == truth

        truth = GroupGraph()
        truth.add_node('ester', 'C(=O)O', [0, 2])
        truth.add_node('alkene_secondary_amine', 'CNC', [0, 0])
        truth.add_node('amine', 'N', [0, 0, 0])
        truth.add_edge((0, 1), (2, 0))
        truth.add_edge((1, 1), (2, 1))

        out = fragment('CNCNOC=O', node_defs)

        assert out == truth


