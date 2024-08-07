from molGrouper.group_graph import GroupGraph
from group_selfies import Group
import networkx as nx
import pytest
import numpy as np

class TestGroupGraph(BaseTest):

    @pytest.fixture(autouse=True)
    def test_creation(self):
        # Define node types with ports
        node_types = {
            'type1': ['port1', 'port2'],
            'type2': ['port3', 'port4'],
            'r1' : ['C1', 'C2', 'C3', 'C4', 'C5', 'C6'],
            'd6' : ['C1', 'C2'], # alkane
            'c' : ['C1', 'C2', 'C3', 'C4'],
            'r3': ['C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22'],
            's1': ['C1'],
            'sulfonate_ester': ['C1', 'O1'],
        }
        self.node_types = node_types

        # Create an instance of GroupGraph for testing
        self.graph = GroupGraph(node_types)

