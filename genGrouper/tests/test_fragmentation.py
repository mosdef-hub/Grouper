from genGrouper.tests.base_test import BaseTest
from genGrouper import GroupGraph, fragment, Node
# from pysmiles import write_smiles
# import pytest

class TestGroupGraph(BaseTest):
    def test_fragment(self):
        node_defs = {}
        node_defs['[OX2H]'] = Node(0, 'hydroxyl', 'O', [0], [0])
        node_defs['C=C'] = Node(1, 'alkene', 'C=C', [0,1,2,3], [0,0,1,1])

        truth = GroupGraph()
        truth.add_node('hydroxyl', '[OH]', [0], [0])
        truth.add_node('alkene', 'C=C', [0,1,2,3], [0,0,1,1])
        truth.add_edge((0,0),(1,0))

        smiles = "C=CO"


        out = fragment(smiles, node_defs)
        print(truth)
        print(out)

        assert out == truth


