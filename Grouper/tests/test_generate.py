import json
import logging
import pathlib

import pytest

from Grouper import Group, exhaustive_generate, random_generate
from Grouper.tests.base_test import BaseTest


class TestGeneration(BaseTest):
    @pytest.mark.skippable(reason="Too slow for general testing")
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
        # Convert node_defs to the expected format
        node_defs = set(Group(n["type"], n["smiles"], n["hubs"]) for n in node_defs)

        logging.info("Starting benchmark")

        benchmark(
            exhaustive_generate,
            n_nodes,
            node_defs,
            num_procs,
        )
        logging.info("Benchmark complete")

    def test_simple_exhaustive_generation(self):
        node_defs = [
            {"type": "t2", "smarts": "N", "hubs": [0, 0, 0]},
            {"type": "Methyl", "smarts": "C", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "C(=O)O", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "O", "hubs": [0, 0]},
        ]

        # Convert node_defs to the expected format
        node_defs = set(
            Group(n["type"], n["smarts"], n["hubs"], pattern_type="SMILES")
            for n in node_defs
        )
        input_file_path = ""
        positive_constraints = {}
        negative_constraints = set()

        logging.info("Starting simple generation")
        exhaustive_generate(
            2,
            node_defs,
            1,  # num_procs
            input_file_path,
            positive_constraints,
            negative_constraints,
            "",
        )
        logging.info("Simple generation complete")

    def test_random_generation(self):
        node_defs = [
            {"type": "t2", "smarts": "N", "hubs": [0, 0, 0]},
            {"type": "Methyl", "smarts": "C", "hubs": [0, 0, 0]},
            {"type": "ester", "smarts": "C(=O)(O)", "hubs": [0, 2]},
            {"type": "extra1", "smarts": "O", "hubs": [0, 0]},
        ]
        # Convert node_defs to the expected format
        node_defs = set(
            Group(n["type"], n["smarts"], n["hubs"], pattern_type="SMILES")
            for n in node_defs
        )

        logging.info("Created node_defs")

        positive_constraints = {}
        negative_constraints = set()

        logging.info("Starting random generation")
        print("Starting random generation")
        random_generate(
            2,  # n_nodes
            node_defs,  # node_defs
            5,  # n_structures
            -1,  # num_procs
            positive_constraints,  # positive
            negative_constraints,  # negative
        )
        print("Random generation complete")
        logging.info("Random generation complete")
