import logging
import time

import pytest

from Grouper import Group, exhaustive_generate, random_generate
from Grouper.tests.base_test import BaseTest


class TestGeneration(BaseTest):
    def test_positive_constraint_filters_results(self):
        """exhaustive_generate must drop graphs that lack the required node types."""
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0, 0]),
        }
        unconstrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
            vcolg_output_file="",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
        )
        constrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
            vcolg_output_file="",
            positive_constraints={"amine": 1},
            negative_constraints=set(),
            config_path="",
        )
        assert len(constrained) > 0, (
            "Positive constraint of one amine should still leave valid graphs"
        )
        assert len(constrained) < len(unconstrained), (
            "Constrained set should drop the methyl-only graphs"
        )
        for graph in constrained:
            n_amines = sum(
                1 for node in graph.nodes.values() if node.type == "amine"
            )
            assert n_amines >= 1, (
                f"Graph {graph.to_smiles()} violates amine>=1 constraint"
            )

    def test_negative_constraint_filters_results(self):
        """negative_constraints must drop SMILES containing forbidden substrings.

        Compares against an unconstrained run on the same inputs to keep
        the test non-vacuous: if the unconstrained set contained no
        NN-bearing graphs to begin with, the constrained set passing the
        substring check would prove nothing.
        """
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0, 0]),
        }
        unconstrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
            vcolg_output_file="",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
        )
        constrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
            vcolg_output_file="",
            positive_constraints={},
            negative_constraints={"NN"},
            config_path="",
        )
        unc_smiles = {g.to_smiles() for g in unconstrained}
        assert any("NN" in s for s in unc_smiles), (
            "Vacuous test case: unconstrained set has no NN-bearing graphs "
            "to filter — pick a forbidden substring that actually appears."
        )
        assert len(constrained) < len(unconstrained), (
            f"negative_constraints did not filter anything: "
            f"|constrained|={len(constrained)} == |unconstrained|={len(unconstrained)}"
        )
        for graph in constrained:
            smiles = graph.to_smiles()
            assert "NN" not in smiles, (
                f"Graph {smiles} contains the forbidden substring 'NN'"
            )

    def test_exhaustive_generate_runtime_alarm(self):
        """Explosion alarm: a small generate run must finish under a generous wall budget.

        Catches catastrophic regressions (≫10x slowdowns or runaway enumeration)
        without being flaky on shared CI runners. The threshold is intentionally
        loose; use --no-skip on test_exhaustive_generate_performance for fine-
        grained perf tracking against pytest-benchmark baselines in .benchmarks/.
        """
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0]),
            Group("ester", "C(=O)O", [0, 2]),
            Group("hydroxyl", "O", [0, 0]),
        }
        start = time.perf_counter()
        result = exhaustive_generate(
            n_nodes=4,
            node_defs=node_defs,
            num_procs=1,
            vcolg_output_file="",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
        )
        elapsed = time.perf_counter() - start
        assert len(result) > 0, "n=4 over four small groups should produce graphs"
        assert elapsed < 30.0, (
            f"exhaustive_generate(n=4) took {elapsed:.1f}s "
            f"(alarm threshold 30s; suggests an algorithmic regression)"
        )


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
