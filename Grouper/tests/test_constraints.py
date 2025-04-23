import pytest
from Grouper import Atom, AtomGraph, Group, GroupGraph, exhaustive_generate
from Grouper.tests.base_test import BaseTest

class TestConstraints(BaseTest):
    def setup_method(self):
        # Define some basic Group instances for reuse
        self.group_A = Group(name="A", atoms=[Atom("C")])
        self.group_B = Group(name="B", atoms=[Atom("O")])
        self.group_C = Group(name="C", atoms=[Atom("N")])
        self.node_defs = {self.group_A, self.group_B, self.group_C}

    def test_exactly_one_group(self):
        results = exhaustive_generate(
            n_nodes=3,
            node_defs=self.node_defs,
            num_procs=1,
            vcolg_output_file="/tmp/vcolg",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
            rule_expressions=["count_A + count_B == 1"]
        )
        for g in results:
            count = sum(g.count_group(name) for name in ["A", "B"])
            assert count == 1

    def test_at_least_one_group(self):
        results = exhaustive_generate(
            n_nodes=3,
            node_defs=self.node_defs,
            num_procs=1,
            vcolg_output_file="/tmp/vcolg",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
            rule_expressions=["count_C >= 1"]
        )
        for g in results:
            assert g.count_group("C") >= 1

    def test_forbid_connection(self):
        results = exhaustive_generate(
            n_nodes=3,
            node_defs=self.node_defs,
            num_procs=1,
            vcolg_output_file="/tmp/vcolg",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
            rule_expressions=["not connected_A_B"]
        )
        for g in results:
            assert not g.is_connected("A", "B")

    def test_cardinality_limit(self):
        results = exhaustive_generate(
            n_nodes=3,
            node_defs=self.node_defs,
            num_procs=1,
            vcolg_output_file="/tmp/vcolg",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
            rule_expressions=["count_B <= 2"]
        )
        for g in results:
            assert g.count_group("B") <= 2

    def test_combined_constraint(self):
        results = exhaustive_generate(
            n_nodes=4,
            node_defs=self.node_defs,
            num_procs=1,
            vcolg_output_file="/tmp/vcolg",
            positive_constraints={},
            negative_constraints=set(),
            config_path="",
            rule_expressions=[
                "count_A + count_B == 2",
                "count_C >= 1",
                "not connected_B_C"
            ]
        )
        for g in results:
            count_ab = g.count_group("A") + g.count_group("B")
            assert count_ab == 2
            assert g.count_group("C") >= 1
            assert not g.is_connected("B", "C")
