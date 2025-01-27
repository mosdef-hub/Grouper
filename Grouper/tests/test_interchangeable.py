from Grouper import GroupGraph, Node, fragment
from Grouper.tests.base_test import BaseTest

# from pysmiles import write_smiles
# import pytest


class TestConversion(BaseTest):
    def test_smiles_gg_to_smiles(self):
        test_cases = [
            "O=COOC(=O)O",
            "CNCNOC=O",
            "O=C(O)OC(=O)O",
            "NNO",
            "O=C(O)C(=O)OO",
            "O=COC(=O)OO",
            "CNO",
            "CCO",
            "COC(N)=O",
            "ONO",
            "CNCC(=O)OC",
            "CNC(N)C(=O)O",
            "CC(=O)ON",
            "CNCOC(C)=O",
            "CNC(O)O",
            "CNCOC(=O)CNC",
            "CNC(=O)O",
            "CNOC=O",
            "NCN",
        ]
        node_defs = set()
        node_defs.add(Node('hydroxyl', 'O', [0]))
        node_defs.add(Node('keto', 'O', [0,0]))
        node_defs.add(Node('ester', 'C(=O)O', [0,0,2]))
        node_defs.add(Node('methyl', 'C', [0,0,0]))
        node_defs.add(Node('t2', 'N', [0,0,0]))
        node_defs.add(Node('secondary_amine', 'CNC', [0,0,1, 2,2,2]))

        for smiles in test_cases:
            print(smiles)
            gg = fragment(smiles, node_defs)
            assert smiles == gg.to_smiles()
            