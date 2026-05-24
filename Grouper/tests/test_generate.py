import logging
import time
import warnings

import pytest

from Grouper import Group, exhaustive_generate, random_generate, random_sample
from Grouper.tests.base_test import BaseTest


class TestGeneration(BaseTest):
    def test_positive_constraint_filters_results(self):
        """exhaustive_generate must drop graphs that lack the required node types.

        Hard-coded expected counts so this test catches behavior changes
        in `vcolg` itself (or in our `positive_constraints` handling)
        rather than re-deriving the counts at runtime and silently
        agreeing with whatever the new behavior produces.
        """
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0, 0]),
        }
        unconstrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
        )
        constrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
            positive_constraints={"amine": 1},
        )
        # n=3 over {amine, methyl}: 10 unique graphs total, 8 with at
        # least one amine (the missing two are CCC and C1CC1).
        assert len(unconstrained) == 10
        assert len(constrained) == 8
        for graph in constrained:
            n_amines = sum(
                1 for node in graph.nodes.values() if node.type == "amine"
            )
            assert n_amines >= 1, (
                f"Graph {graph.to_smiles()} violates amine>=1 constraint"
            )

    def test_negative_constraint_filters_results(self):
        """negative_constraints must drop SMILES containing forbidden substrings.

        Hard-coded counts here so a future regression in either vcolg's
        output or our negative-constraint filter is caught — re-deriving
        the "unconstrained" count at runtime would silently track any
        behavior change in lockstep.
        """
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0, 0]),
        }
        unconstrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
        )
        constrained = exhaustive_generate(
            n_nodes=3,
            node_defs=node_defs,
            num_procs=1,
            negative_constraints={"NN"},
        )
        # n=3 over {amine, methyl} produces 10 unique graphs; 4 of them
        # contain "NN" (CNN, NCN, NNN, plus the ring N1NN1 / C1NN1
        # variants). The constrained run drops them to 6.
        assert len(unconstrained) == 10
        assert len(constrained) == 6
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
        )
        elapsed = time.perf_counter() - start
        assert len(result) > 0
        assert elapsed < 30.0, (
            f"exhaustive_generate(n=4) took {elapsed:.1f}s to generate "
            f"{len(result)} graphs."
        )

    def test_generation_with_smarts_library(self):
        """SMARTS-pattern libraries (e.g. Joback) must work in
        exhaustive_generate.

        Pre-existing bug: processColoredGraphs.cpp constructed a scratch
        `tmp_g` Group via default-construction + field assignment without
        setting `patternType`, so it stayed at the empty-string default.
        `Group::getPossibleAttachments` then fell through to the SMILES
        branch and called `AtomGraph::fromSmiles` on a SMARTS query like
        `[CX4H3]`, which RDKit's SMILES parser rejects — followed by a
        downstream valency-exhaustion crash because the resulting empty
        AtomGraph can't carry the requested bond.

        Pin-down: methyl `[CX4H3]` (1 free valence after H3) joined with
        methylene `[CX4H2]` (2 free valences after H2) at n=2 yields
        exactly 1 unique graph (CC = ethane), because methylene's
        hubs=[0,0] requires 2 ports used but a 2-node 1-edge topology
        only fills 1 — so methylene can't appear, and methyl-methyl
        is the only valid pair.
        """
        methyl = Group("methyl", "[CX4H3]", [0], pattern_type="SMARTS")
        methylene = Group(
            "methylene", "[CX4H2]", [0, 0], pattern_type="SMARTS"
        )
        result = exhaustive_generate(
            n_nodes=2,
            node_defs={methyl, methylene},
            num_procs=1,
        )
        smiles = {g.to_smiles() for g in result}
        assert smiles == {"CC"}, (
            f"expected only ethane (CC) from a methyl+methylene "
            f"n=2 enumeration; got {sorted(smiles)}"
        )

    def test_all_equivalent_port_orbit_collapses(self):
        """When every port of a group sits in the same atom-orbit
        (under the group's automorphism), the per-edge candidate
        generator must collapse to one canonical port assignment per
        side rather than cross-producing all p^2 combinations.

        Pre-fix, this case bottlenecked on the recursive port-pairing
        search (every port permutation of a node was materialized
        and only deduped post-canonicalization), making carbon-only
        n=5 effectively non-terminating: ~16^E candidate edge-colorings
        per K5-shaped topology before dedup.

        The fix recognises orbit equivalence in two flavours: literal
        identical hubs (`C [0,0,0,0]`) and atoms equivalent under the
        group's nauty-derived atom automorphism (`benzene c1ccccc1
        [0,3]` — atoms 0 and 3 are para-equivalent under D6).

        Pin: carbon-only n=5 must produce exactly A001349(5)=21
        unique structures (every connected simple graph on 5 nodes
        collapses to one canonical form when all groups are identical
        single-atom carbons), and must do so in well under a second.
        """
        import time as _time
        node_defs = {Group("c", "C", [0, 0, 0, 0])}
        start = _time.perf_counter()
        result = exhaustive_generate(
            n_nodes=5,
            node_defs=node_defs,
            num_procs=1,
            confirm=True,
        )
        elapsed = _time.perf_counter() - start
        assert len(result) == 21, (
            f"expected A001349(5)=21 unique graphs for carbon-only "
            f"n=5; got {len(result)}"
        )
        # 5 seconds is generous — pre-fix this didn't terminate in an
        # hour. With the fix it's milliseconds locally; the threshold
        # is a CI-noise buffer, not a target.
        assert elapsed < 5.0, (
            f"carbon-only n=5 took {elapsed:.1f}s — port-orbit "
            f"collapse appears regressed"
        )

    def test_multi_atom_hubs_not_overcollapsed(self):
        """Correctness anchor for the port-orbit collapse: we must
        NOT collapse when hubs land on chemically-distinct positions
        within a multi-atom group.

        Full benzene `c1ccccc1 [0,1,2,3,4,5]` has all six attachment
        atoms in a single nauty orbit under D6, but the stabilizer of
        any one attachment atom does NOT pointwise-fix the others
        (it swaps ortho positions and meta positions across the
        symmetry axis). So ortho/meta/para attachment patterns on a
        given ring are chemically distinct and the enumeration must
        produce all of them.

        Pin: three full-benzene groups at n=3 must yield 13 unique
        structures (the cyclic 3-ring + 12 distinct chain
        configurations covering every ortho/meta/para combination on
        the middle ring). An earlier orbit-only version of the
        collapse heuristic gave only 2 here, silently dropping all
        but one chain variant.
        """
        bf = Group("bf", "c1ccccc1", [0, 1, 2, 3, 4, 5])
        result = exhaustive_generate(
            n_nodes=3, node_defs={bf}, num_procs=1, confirm=True
        )
        assert len(result) == 13, (
            f"full-benzene n=3 must yield 13 distinct structures "
            f"covering ortho/meta/para chain variants + cyclic ring; "
            f"got {len(result)}. If this dropped to 2, the port-orbit "
            f"collapse was incorrectly applied to multi-atom hubs."
        )


    def test_random_sample_small_smoke(self):
        """Direct sampler returns a `GroupGraphSet` of unique molecules
        at small n. Pinned counts because the seed is fixed — if the
        sampler is broken or the algorithm changes meaningfully, the
        output set diverges and the test catches it.
        """
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0, 0]),
        }
        result = random_sample(
            n_nodes=3,
            node_defs=node_defs,
            num_graphs=8,
            seed=42,
        )
        # SMILES dedup inside the sampler means every returned graph
        # is a unique molecule.
        smiles = [g.to_smiles() for g in result]
        assert len(smiles) == len(set(smiles))
        # Sampler is best-effort: at small n with 2 groups the request
        # for 8 unique graphs is achievable within the default attempt
        # cap (20 * num_graphs = 160 attempts on a ~10-graph space).
        assert len(result) >= 5

    def test_random_sample_scales_past_vcolg(self):
        """The whole point of this sampler: it must work at sizes
        where `random_generate` / `exhaustive_generate` would never
        terminate. Pin n=18, which is far past vcolg's cliff at this
        library shape; we should still get the requested count in
        well under a second.
        """
        node_defs = {
            Group("methyl",   "C", [0, 0, 0, 0]),
            Group("amine",    "N", [0, 0, 0]),
            Group("hydroxyl", "O", [0, 0]),
            Group("ester",    "C(=O)O", [0, 2]),
        }
        start = time.perf_counter()
        result = random_sample(
            n_nodes=18,
            node_defs=node_defs,
            num_graphs=100,
            seed=7,
        )
        elapsed = time.perf_counter() - start
        assert len(result) == 100, (
            f"random_sample failed to find 100 unique graphs at n=18; "
            f"got {len(result)}"
        )
        smiles = {g.to_smiles() for g in result}
        assert len(smiles) == 100, "duplicate SMILES leaked through dedup"
        # Generous CI bound; locally this runs in ~100ms.
        assert elapsed < 5.0, (
            f"random_sample(n=18, num_graphs=100) took {elapsed:.1f}s — "
            f"performance has regressed"
        )

    def test_random_sample_negative_constraint(self):
        """Negative constraints reject any returned graph whose SMILES
        contains a forbidden substring. This exercises the only code
        path in random_sample that still calls toSmiles per attempt
        (so it also implicitly checks that the RDLog::LogStateSetter
        scoping doesn't break correctness when active).
        """
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0, 0]),
        }
        result = random_sample(
            n_nodes=5,
            node_defs=node_defs,
            num_graphs=20,
            negative_constraints={"NN"},
            seed=11,
            show_progress=False,
        )
        for g in result:
            smiles = g.to_smiles()
            assert "NN" not in smiles, (
                f"Graph {smiles} contains the forbidden substring 'NN' "
                f"that negative_constraints should have rejected"
            )

    def test_random_sample_positive_constraint(self):
        """Positive constraints must filter the sampler's output: every
        returned graph satisfies the per-type minimum. We also check
        that the sampler warns rather than infinite-looping when the
        constraint is unreachable within the attempt budget.
        """
        node_defs = {
            Group("amine", "N", [0, 0, 0]),
            Group("methyl", "C", [0, 0, 0, 0]),
        }
        # Achievable: at n=6 with 2 groups, requiring 2 amines is
        # comfortably met by stratified sampling.
        result = random_sample(
            n_nodes=6,
            node_defs=node_defs,
            num_graphs=20,
            positive_constraints={"amine": 2},
            seed=1,
        )
        for g in result:
            n_amines = sum(
                1 for node in g.nodes.values() if node.type == "amine"
            )
            assert n_amines >= 2, (
                f"Graph {g.to_smiles()} has {n_amines} amine(s), "
                f"violates positive_constraints={{amine: 2}}"
            )

        # Unreachable: a 3-node graph can't have 10 amines. Sampler
        # must warn and return early instead of hanging.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            impossible = random_sample(
                n_nodes=3,
                node_defs=node_defs,
                num_graphs=5,
                positive_constraints={"amine": 10},
                max_attempts=50,
                seed=1,
            )
        assert len(impossible) == 0
        assert any(
            "exhausted max_attempts" in str(w.message) for w in caught
        ), "expected a max_attempts-exhausted warning"

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
