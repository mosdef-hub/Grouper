"""Tests for GroupGraphSet — the wrapper around exhaustive_generate
results that adds DataFrame / SMILES-list / SDF helpers."""

import warnings

import pytest

from Grouper import Group, GroupGraph, exhaustive_generate, random_generate
from Grouper.group_graph_set import GroupGraphSet


# ---------------------------------------------------------------------
# Compatibility: existing code that treats the result as a set must
# keep working unchanged. GroupGraphSet IS a set; these tests pin
# that contract.
# ---------------------------------------------------------------------
def test_is_a_set_subclass():
    """GroupGraphSet must be a subclass of set so isinstance checks
    in existing code keep returning True."""
    results = exhaustive_generate(2, {Group("methyl", "C", [0, 0, 0, 0])})
    assert isinstance(results, GroupGraphSet)
    assert isinstance(results, set)


def test_iteration_and_len():
    """Standard Pythonic patterns: iter, len, in, comprehensions."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    n_iter = sum(1 for _ in results)
    assert n_iter == len(results)
    assert n_iter > 0
    smiles_set = {g.to_smiles() for g in results}
    assert len(smiles_set) == n_iter


def test_random_generate_also_wrapped():
    """random_generate is also wrapped — both entry points return the
    same convenience type."""
    groups = {
        Group("methyl", "C", [0, 0, 0]),
        Group("hydroxyl", "O", [0, 0]),
    }
    results = random_generate(3, groups, 5, 1)
    assert isinstance(results, GroupGraphSet)


# ---------------------------------------------------------------------
# to_smiles_list — quickest of the conversions, used to anchor
# downstream pandas / numpy work.
# ---------------------------------------------------------------------
def test_to_smiles_list():
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    smis = results.to_smiles_list()
    assert isinstance(smis, list)
    assert len(smis) == len(results)
    assert all(isinstance(s, str) for s in smis)


# ---------------------------------------------------------------------
# to_dataframe — the headline method. Always-present columns + behaviour
# under a sort_by argument and an empty result set.
# ---------------------------------------------------------------------
def test_to_dataframe_basic():
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    df = results.to_dataframe()
    assert "smiles" in df.columns
    assert "n_nodes" in df.columns
    assert len(df) == len(results)
    assert df["n_nodes"].eq(2).all()


def test_to_dataframe_sort_by_smiles():
    """sort_by="smiles" returns a deterministic order — important for
    downstream comparisons since sets are unordered."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    df = results.to_dataframe(sort_by="smiles")
    assert list(df["smiles"]) == sorted(df["smiles"])


def test_to_dataframe_include_groupgraph():
    """include_groupgraph=True puts the GroupGraph objects in a column."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    df = results.to_dataframe(include_groupgraph=True)
    assert "graph" in df.columns
    assert all(isinstance(g, GroupGraph) for g in df["graph"])


def test_to_dataframe_empty_set():
    """An empty GroupGraphSet produces an empty DataFrame, not an
    error. Belt-and-suspenders since pandas handles this fine, but
    the assertion catches a regression where we'd start synthesizing
    columns from a nonexistent first row."""
    empty = GroupGraphSet()
    df = empty.to_dataframe()
    assert len(df) == 0


# ---------------------------------------------------------------------
# to_dataframe with property estimators — the killer-app feature. We
# pass a synthetic estimator class to avoid coupling these tests to
# the optional Grouper.properties module (which lives on a sibling
# branch and won't always be present).
# ---------------------------------------------------------------------
class _DummyEstimator:
    """Minimal property-estimator stand-in for tests. Mimics
    PropertyEstimate's keys() / __getitem__ contract."""

    method_name = "dummy"

    def __init__(self, gg):
        self._smiles = gg.to_smiles()
        self._n = len(gg.nodes)

    def keys(self):
        return ["heavy_atoms", "label"]

    def __getitem__(self, k):
        if k == "heavy_atoms":
            return self._n  # not strictly heavy atoms, but unique-per-graph
        if k == "label":
            return f"dummy_{self._smiles}"
        raise KeyError(k)

    @classmethod
    def from_group_graph(cls, gg):
        return cls(gg)


def test_to_dataframe_with_estimator_class():
    """Passing a class directly (not a registry name) computes its
    properties for each graph and prefixes the columns with the
    estimator's method_name."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    df = results.to_dataframe(properties=[_DummyEstimator])
    assert "dummy.heavy_atoms" in df.columns
    assert "dummy.label" in df.columns
    assert df["dummy.heavy_atoms"].eq(2).all()


def test_to_dataframe_unknown_string_estimator_raises():
    """A string name not in the registry must raise — not silently
    produce a DataFrame missing the requested columns."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    with pytest.raises(ValueError, match="unknown property estimator"):
        results.to_dataframe(properties=["nonsense_method_name"])


def test_to_dataframe_skip_failures_warns():
    """skip_failures=True turns property-estimator exceptions into
    warnings + None values, instead of aborting the whole DataFrame."""

    class _ExplodingEstimator:
        method_name = "boom"

        @classmethod
        def from_group_graph(cls, gg):
            raise RuntimeError("kaboom")

    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        df = results.to_dataframe(
            properties=[_ExplodingEstimator], skip_failures=True
        )
    # All molecules failed — every warning should have fired.
    assert len(caught) == len(results)
    assert all("kaboom" in str(w.message) for w in caught)
    # No explosion column should exist (every value would be NaN anyway).
    assert "boom.x" not in df.columns


# ---------------------------------------------------------------------
# to_sdf — batch SDF write through GroupGraphSet, with the same
# property-attaching semantics as to_dataframe.
# ---------------------------------------------------------------------
def test_to_sdf_writes_multimol_file(tmp_path):
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    out = tmp_path / "results.sdf"
    n = results.to_sdf(str(out))
    assert n == len(results)
    assert out.stat().st_size > 0
    # Should be a valid multi-mol SDF readable by RDKit
    from rdkit import Chem
    supplier = Chem.SDMolSupplier(str(out))
    mols = [m for m in supplier if m is not None]
    assert len(mols) == n


def test_to_sdf_with_estimator_attaches_properties(tmp_path):
    """Passing properties=[_DummyEstimator] should attach
    `dummy.heavy_atoms` and `dummy.label` SD tags to every molecule
    in the output SDF."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    out = tmp_path / "with_props.sdf"
    n = results.to_sdf(str(out), properties=[_DummyEstimator])
    assert n == len(results)

    from rdkit import Chem
    supplier = Chem.SDMolSupplier(str(out))
    mols = [m for m in supplier if m is not None]
    for m in mols:
        # Attached properties must round-trip through the SDF.
        assert m.GetProp("dummy.heavy_atoms") == "2"
        assert m.GetProp("dummy.label").startswith("dummy_")


# ---------------------------------------------------------------------
# filter — returns another GroupGraphSet (not a plain set), so chained
# calls keep the convenience methods.
# ---------------------------------------------------------------------
def test_filter_returns_groupgraphset():
    """`results.filter(pred)` must preserve the wrapper type so
    chained calls like `.filter(p1).filter(p2).to_dataframe()` work."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    only_with_O = results.filter(lambda g: "O" in g.to_smiles())
    assert isinstance(only_with_O, GroupGraphSet)
    assert all("O" in g.to_smiles() for g in only_with_O)
    # Chained: the result should still have to_dataframe.
    df = only_with_O.to_dataframe()
    assert all("O" in s for s in df["smiles"])


# ---------------------------------------------------------------------
# repr — should be informative without dumping every element of huge
# result sets into the terminal.
# ---------------------------------------------------------------------
def test_repr_truncates_for_large_sets():
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    r = repr(results)
    assert "GroupGraphSet" in r
    assert str(len(results)) in r


def test_repr_empty():
    assert repr(GroupGraphSet()) == "GroupGraphSet(empty)"


# ---------------------------------------------------------------------
# sample — random subset, with and without seed.
# ---------------------------------------------------------------------
def test_sample_returns_groupgraphset_of_right_size():
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    n = min(2, len(results))
    sub = results.sample(n)
    assert isinstance(sub, GroupGraphSet)
    assert len(sub) == n
    # Sampled elements must all come from the parent set.
    assert sub.issubset(results)


def test_sample_seed_is_reproducible():
    """Same seed -> same sample. Sets have no inherent order, but
    `random.sample(seeded_rng, list(set), k)` is deterministic for a
    given Python build because list(set) iteration order is fixed
    once the set is built."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    if len(results) < 2:
        pytest.skip("need at least 2 elements to test reproducibility")
    a = results.sample(2, seed=42)
    b = results.sample(2, seed=42)
    assert a == b


def test_sample_too_many_raises():
    """Asking for more than len(self) must raise — matches the
    `random.sample` semantic."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    with pytest.raises(ValueError):
        results.sample(len(results) + 5)


# ---------------------------------------------------------------------
# to_csv / to_jsonl — pass-throughs to to_dataframe + serializer.
# ---------------------------------------------------------------------
def test_to_csv_round_trip(tmp_path):
    """CSV must be readable back via pandas with the same columns."""
    import pandas as pd
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    out = tmp_path / "results.csv"
    results.to_csv(str(out))
    df = pd.read_csv(out)
    assert "smiles" in df.columns
    assert "n_nodes" in df.columns
    assert len(df) == len(results)


def test_to_csv_default_sort(tmp_path):
    """Default sort_by='smiles' makes the file deterministic across
    runs. Two calls produce the same bytes."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    p1 = tmp_path / "a.csv"
    p2 = tmp_path / "b.csv"
    results.to_csv(str(p1))
    results.to_csv(str(p2))
    # Read bytes through context managers so file handles get released
    # promptly (and CodeQL stops flagging an asserted open()-no-close).
    contents_p1 = p1.read_text()
    contents_p2 = p2.read_text()
    assert contents_p1 == contents_p2


def test_to_jsonl_round_trip(tmp_path):
    """JSONL must be readable back via pandas with one row per line."""
    import pandas as pd
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0]), Group("hydroxyl", "O", [0, 0])}
    )
    out = tmp_path / "results.jsonl"
    results.to_jsonl(str(out))
    df = pd.read_json(out, lines=True)
    assert "smiles" in df.columns
    assert "n_nodes" in df.columns
    assert len(df) == len(results)


def test_to_csv_with_estimator(tmp_path):
    """CSV with property columns: round-trip preserves the data."""
    import pandas as pd
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    out = tmp_path / "with_props.csv"
    results.to_csv(str(out), properties=[_DummyEstimator])
    df = pd.read_csv(out)
    assert "dummy.heavy_atoms" in df.columns
    assert "dummy.label" in df.columns


def test_to_csv_returns_path(tmp_path):
    """Path is returned for chained ops (e.g. logging the destination)."""
    results = exhaustive_generate(
        2, {Group("methyl", "C", [0, 0, 0, 0])}
    )
    out = tmp_path / "results.csv"
    returned = results.to_csv(str(out))
    assert returned == str(out)
