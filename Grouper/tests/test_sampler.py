import logging

from Grouper import Node, GroupGraph, Sampler
from Grouper.tests.base_test import BaseTest


class TestGeneration(BaseTest):
    def test_simple_exhaustive_generation(self):
        node_defs = [
            {"type": "t2", "smarts": "[N]", "hubs": [0, 0, 0]},
            {"type": "Methyl", "smarts": "[CX3]", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "[CX3](=[OX1])([OX2])", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "[O]", "hubs": [0, 0]},
        ]
        # Create initial graph
        gG = GroupGraph()
        for node_def in node_defs:
            gG.add_node(node_def['type'], node_def['smarts'], node_def['hubs'])
        # Add edges
        gG.add_edge((0,0), (1,0))
        gG.add_edge((1,1), (2,0))
        gG.add_edge((2,1), (3,0))

        # Create sampler
        sampler = Sampler(1.0, 100)

        # Perform monte carlo sampling
        sampler.sample(gG)
