import time
import pytest
from genGrouper import GroupGraph, exhaustive_generate, Node
from genGrouper.tests.base_test import BaseTest
import logging

class TestGeneration(BaseTest):

    @pytest.mark.parametrize("n_nodes", [2,3, 4, 5])
    @pytest.mark.parametrize("num_procs", [1, 2, 4, 8])
    @pytest.mark.parametrize("node_defs", [[
            {"type": "t1", "smiles": "C", "ports": [0,1,2,3], "hubs": [0,0,0,0]},
            {"type": "t2", "smiles": "N", "ports": [0,1,2], "hubs": [0,0,0]},
            {"type": "Methyl", "smiles": "C", "ports": [0,1,2], "hubs": [0,0,0]},
            {"type": "Benzene", "smiles": "C1=CC=CC=C1", "ports": [0,1,2,3,4,5], "hubs": [0,1,2,3,4,5]},
            {"type": "ester", "smiles": "C(=O)O", "ports": [0,1], "hubs": [0,2]},
            {"type": "extra1", "smiles": "O", "ports": [0,1], "hubs": [0,0]}
        ]]
    )
    def test_exhaustive_generate_performance(self, benchmark, n_nodes, num_procs, node_defs):
        # Convert node_defs to the expected format
        node_defs = set(Node(0, n["type"], n["smiles"], n["ports"], n["hubs"]) for n in node_defs)

        write_to_db = False
        verbose = False
        # nauty_path = "/Users/kieran/projects/genGrouper/packages/nauty2_8_8"
        nauty_path = "/raid6/homes/kierannp/projects/molGrouper/packages/nauty2_8_8"
        input_file_path = ""
        positive_constraints = {}
        negative_constraints = set()

        logging.info("Starting benchmark")

        benchmark(
            exhaustive_generate, 
            n_nodes, 
            node_defs, 
            nauty_path, 
            input_file_path, 
            num_procs, 
            positive_constraints, 
            negative_constraints, 
            config_path="",
            verbose=verbose,
        )
        logging.info("Benchmark complete")

