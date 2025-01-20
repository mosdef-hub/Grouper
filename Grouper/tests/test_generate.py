import logging

import pytest

from Grouper import Node, exhaustive_generate
from Grouper.tests.base_test import BaseTest


class TestGeneration(BaseTest):
    @pytest.mark.skip(reason="Too slow for general testing")
    @pytest.mark.parametrize("n_nodes", [2, 3, 4, 5, 6])
    @pytest.mark.parametrize("num_procs", [1, 2, 4, 8, 16])
    @pytest.mark.parametrize(
        "node_defs",
        [
            [
                # {"type": "t1", "smiles": "C", "ports": [0,1,2,3], "hubs": [0,0,0,0]},
                {"type": "t2", "smiles": "N", "hubs": [0, 0, 0]},
                {"type": "Methyl", "smiles": "C", "hubs": [0, 0, 0]},
                # {"type": "Benzene", "smiles": "C1=CC=CC=C1", "ports": [0,1,2,3,4,5], "hubs": [0,1,2,3,4,5]},
                {"type": "ester", "smiles": "C(=O)O", "hubs": [0, 2]},
                {"type": "extra1", "smiles": "O", "hubs": [0, 0]},
            ]
        ],
    )
    def test_exhaustive_generate_performance(
        self, benchmark, n_nodes, num_procs, node_defs
    ):
        # Load nauty path from config file
        import json
        import pathlib
        print(pathlib.Path(__file__).parent.parent.parent / "config.json")
        with open(pathlib.Path(__file__).parent.parent.parent / "config.json", "r") as config_file:

            config = json.load(config_file)
        nauty_path = config.get("nauty_path")
        if not nauty_path:
            raise RuntimeError("Nauty path is not defined in the configuration file.")
        # Convert node_defs to the expected format
        node_defs = set(Node(n["type"], n["smiles"], n["hubs"]) for n in node_defs)

        verbose = False
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
            "",
            verbose,
        )
        logging.info("Benchmark complete")
        
    def test_simple_exhaustive_generation(self):
        node_defs = [
            {"type": "t2", "smarts": "[N]", "hubs": [0, 0, 0]},
            {"type": "Methyl", "smarts": "[CX3]", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "[CX3](=[OX1])([OX2])", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "[O]", "hubs": [0, 0]},
        ]
        # Load nauty path from config file
        import json
        import pathlib
        print(pathlib.Path(__file__).parent.parent.parent / "config.json")
        with open(pathlib.Path(__file__).parent.parent.parent / "config.json", "r") as config_file:
            config = json.load(config_file)
        nauty_path = config.get("nauty_path")
        if not nauty_path:
            raise RuntimeError("Nauty path is not defined in the configuration file.")
        
        # Convert node_defs to the expected format
        node_defs = set(Node(n["type"], n["smarts"], n["hubs"]) for n in node_defs)

        logging.info("Created node_defs")

        verbose = False
        input_file_path = ""
        positive_constraints = {}
        negative_constraints = set()

        logging.info("Starting simple generation")
        print("Starting simple generation")
        exhaustive_generate(
            2,
            node_defs,
            nauty_path,
            input_file_path,
            1, # num_procs
            positive_constraints,
            negative_constraints,
            "",
            verbose,
        )
        print("Simple generation complete")
        logging.info("Simple generation complete")
