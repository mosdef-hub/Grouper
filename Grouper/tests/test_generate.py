import logging
import os
import json
import pathlib

import pytest

from Grouper import Group, exhaustive_generate, random_generate
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
        print(pathlib.Path(__file__).parent.parent.parent / "config.json")
        with open(
            pathlib.Path(__file__).parent.parent.parent / "config.json", "r"
        ) as config_file:
            config = json.load(config_file)
        nauty_path = config.get("nauty_path")
        if not nauty_path:
            raise RuntimeError("Nauty path is not defined in the configuration file.")
        # Convert node_defs to the expected format
        node_defs = set(Group(n["type"], n["smiles"], n["hubs"]) for n in node_defs)

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
            {"type": "Methyl", "smarts": "[C]", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "[C;3](=O)([O;2])", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "[O]", "hubs": [0, 0]},
        ]
        # Load nauty path from config file
        import json
        import pathlib

        packagePath = pathlib.Path(__file__).parent.parent.parent
        print(packagePath / "config.json")
        with open(packagePath / "config.json", "r") as config_file:
            config = json.load(config_file)
        nauty_path = config.get("nauty_path")
        if nauty_path and os.path.exists(nauty_path):
            pass
        elif os.path.exists(packagePath / "packages/nauty"):
            nauty_path = str(packagePath / "packages/nauty")
        else:
            raise RuntimeError(
                f"Nauty path {nauty_path} is not defined in the configuration file or does not exists."
            )

        # Convert node_defs to the expected format
        node_defs = set(Group(n["type"], n["smarts"], n["hubs"], is_smarts=True) for n in node_defs)

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
            1,  # num_procs
            positive_constraints,
            negative_constraints,
            "",
            verbose,
        )
        print("Simple generation complete")
        logging.info("Simple generation complete")

    def test_random_generation(self):
        node_defs = [
            {"type": "t2", "smarts": "N", "hubs": [0, 0, 0]},
            {"type": "Methyl", "smarts": "C", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "C(=O)(O)", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "O", "hubs": [0, 0]},
        ]

        # Load nauty path from config file


        packagePath = pathlib.Path(__file__).parent.parent.parent
        print(packagePath / "config.json")
        with open(packagePath / "config.json", "r") as config_file:
            config = json.load(config_file)
        nauty_path = config.get("nauty_path")
        if nauty_path and os.path.exists(nauty_path):
            pass
        elif os.path.exists(packagePath / "packages/nauty"):
            nauty_path = str(packagePath / "packages/nauty")
        else:
            raise RuntimeError(
                f"Nauty path {nauty_path} is not defined in the configuration file or does not exists."
            )

        # Convert node_defs to the expected format
        node_defs = set(Group(n["type"], n["smarts"], n["hubs"], is_smarts=False) for n in node_defs)

        logging.info("Created node_defs")

        input_file_path = ""
        positive_constraints = {}
        negative_constraints = set()

        logging.info("Starting random generation")
        print("Starting random generation")
        random_generate(
            2, # n_nodes
            node_defs, # node_defs
            5, # n_structures
            -1,  # num_procs
            nauty_path, # nauty_path
            positive_constraints, # positive
            negative_constraints, # negative
        )
        print("Random generation complete")
        logging.info("Random generation complete")