"""Result-set wrapper for `exhaustive_generate` / `random_generate` /
`exhaustive_fragment`.

The C++ bindings return a Python `set[GroupGraph]`. This module wraps
that set in a `GroupGraphSet` subclass so common downstream operations
read as method calls instead of one-line boilerplate that gets
copy-pasted across every notebook:

    results = exhaustive_generate(3, groups)

    results.to_dataframe()                         # one row per graph
    results.to_dataframe(properties=["joback"])    # + Joback predictions
    results.to_smiles_list()                       # ["CCO", ...]
    results.to_sdf("out.sdf")                      # multi-mol SDF
    results.filter(lambda g: len(g.nodes) > 2)     # GroupGraphSet subset

`GroupGraphSet` is a `set` subclass: every existing pattern (iter, len,
membership tests, comprehensions) keeps working unchanged.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from typing import Any, Callable, List, Optional, Type, Union

# ---------------------------------------------------------------------
# Property-estimator registry (lazy import so this module loads even if
# Grouper.properties isn't installed). Keys are short string names a
# user can pass to to_dataframe / to_sdf; values are import paths
# resolved on first use.
# ---------------------------------------------------------------------
_REGISTRY: dict[str, str] = {
    "joback": "Grouper.properties.joback:JobackEstimate",
    # Future: "unifac": "Grouper.properties.unifac:UnifacEstimate", etc.
}


def _resolve_estimator(spec: Union[str, Type]) -> Type:
    """Turn a string name or a class into a property-estimator class.

    Lazy import means a user without `Grouper.properties` installed can
    still use the rest of the wrapper (smiles list, filter, plain SDF)
    without paying an ImportError.
    """
    if isinstance(spec, str):
        if spec not in _REGISTRY:
            known = ", ".join(repr(k) for k in _REGISTRY)
            raise ValueError(
                f"unknown property estimator {spec!r}; known: {known}. "
                f"Pass a class directly to use a custom estimator."
            )
        module_path, class_name = _REGISTRY[spec].split(":")
        try:
            module = __import__(module_path, fromlist=[class_name])
        except ImportError as e:
            raise ImportError(
                f"property estimator {spec!r} requires {module_path} "
                f"which is not importable: {e}"
            ) from e
        return getattr(module, class_name)
    if isinstance(spec, type):
        return spec
    raise TypeError(
        f"property spec must be a string name or a class; got {type(spec).__name__}"
    )


class GroupGraphSet(set):
    """A `set[GroupGraph]` with batch-conversion methods.

    Returned by `exhaustive_generate`, `random_generate`, and
    `exhaustive_fragment` (instead of a plain `set`). Subclasses `set`
    so every existing pattern — iteration, `len`, `in`, set algebra —
    works unchanged. The added methods are pure conveniences.

    Note that set operations like `.union(other)` from the base class
    return a plain `set`, not a `GroupGraphSet`. Wrap explicitly if you
    need the methods on the union result.
    """

    # ----- conversion to standard data shapes ------------------------

    def to_smiles_list(self) -> List[str]:
        """Return a list of canonical SMILES strings, one per graph.

        Order is not stable across runs (sets are unordered). Sort
        explicitly if you need reproducibility:
        `sorted(results.to_smiles_list())`.
        """
        return [g.to_smiles() for g in self]

    def to_dataframe(
        self,
        properties: Optional[Sequence[Union[str, Type]]] = None,
        include_groupgraph: bool = False,
        skip_failures: bool = False,
        sort_by: Optional[str] = None,
    ):
        """Flatten the set into a pandas DataFrame, one row per graph.

        Always-present columns: `smiles`, `n_nodes`.

        `properties=["joback"]` (or a list of estimator classes) adds
        prefixed columns for each estimator's output — e.g.
        `joback.Tb`, `joback.Tc`, ... so multiple estimators don't
        collide on common field names like `MW` or `n_atoms`.

        `include_groupgraph=True` adds a `graph` column with the
        underlying GroupGraph objects (handy for further processing
        from the DataFrame, at the cost of not being JSON/CSV
        serializable).

        `skip_failures=True` catches exceptions from the property
        estimators (e.g. JobackUnsupportedAtomError on a graph
        containing Si) and writes None into the affected columns.
        Default is to raise so silent garbage doesn't leak through.

        `sort_by="smiles"` (or any column name) returns a sorted
        DataFrame for reproducibility — sets have no inherent order.
        """
        import pandas as pd  # heavy import; defer until used

        estimators = [_resolve_estimator(p) for p in (properties or [])]
        rows = []
        for g in self:
            row = {"smiles": g.to_smiles(), "n_nodes": len(g.nodes)}
            if include_groupgraph:
                row["graph"] = g
            for est in estimators:
                prefix = est.method_name + "."
                try:
                    est_obj = est.from_group_graph(g)
                except Exception as e:
                    if skip_failures:
                        warnings.warn(
                            f"property {est.method_name!r} failed on "
                            f"{row['smiles']!r}: {type(e).__name__}: {e}"
                        )
                        # Leave property columns absent; pandas fills
                        # NaN when the row is concatenated with others.
                        continue
                    raise
                for key in est_obj.keys():
                    row[prefix + key] = est_obj[key]
            rows.append(row)

        df = pd.DataFrame(rows)
        if sort_by is not None and sort_by in df.columns:
            df = df.sort_values(sort_by).reset_index(drop=True)
        return df

    def to_sdf(
        self,
        path: str,
        properties: Optional[Sequence[Union[str, Type]]] = None,
        embed_3d: bool = True,
        skip_failures: bool = True,
        **embed_kwargs,
    ) -> int:
        """Write every graph to a multi-molecule SDF.

        Delegates to `Grouper.exports.to_sdf` with the set as the
        iterable target. If `properties` is given, each molecule's SDF
        entry carries the predicted property values as `> <KEY>`
        blocks — feeds straight into KNIME, RDKit's
        `PandasTools.LoadSDF`, or any tool that reads SD properties.

        Returns the number of molecules successfully written.
        """
        from Grouper.exports import to_sdf

        if properties:
            estimators = [_resolve_estimator(p) for p in properties]

            def attach(g):
                out = {}
                for est in estimators:
                    prefix = est.method_name + "."
                    try:
                        est_obj = est.from_group_graph(g)
                    except Exception:
                        continue  # to_sdf still emits the structure
                    for k in est_obj.keys():
                        out[prefix + k] = est_obj[k]
                return out

            props_arg = attach
        else:
            props_arg = None

        return to_sdf(
            self,
            path=path,
            embed_3d=embed_3d,
            properties=props_arg,
            skip_failures=skip_failures,
            **embed_kwargs,
        )

    # ----- convenience filtering / sampling -------------------------

    def filter(self, predicate: Callable[[Any], bool]) -> GroupGraphSet:
        """Return a new GroupGraphSet containing only graphs for which
        `predicate(graph)` is True. Standard list/set comprehensions
        return a plain set; this returns the wrapper so chained calls
        keep the convenience methods.
        """
        return GroupGraphSet(g for g in self if predicate(g))

    def sample(
        self,
        n: int,
        seed: Optional[int] = None,
    ) -> GroupGraphSet:
        """Return a uniformly-random subset of `n` graphs.

        `seed` makes the sample reproducible — useful when you want
        the same subset across re-runs of a screening notebook. Sets
        are unordered in Python, so the sample's "starting position"
        is whatever the hash table dictates; combined with a seeded
        RNG the result is deterministic for a given Python build.

        Raises ValueError if `n > len(self)` (matching `random.sample`
        semantics so the failure mode is recognizable).
        """
        import random

        rng = random.Random(seed) if seed is not None else random
        # `random.sample` requires a sequence, not a set. Materialize
        # via list — fine for the sample sizes (~thousands) people
        # actually use; if a 10⁷-graph set ever wants this method it
        # warrants a reservoir-sampling pass instead.
        return GroupGraphSet(rng.sample(list(self), n))

    # ----- file output ----------------------------------------------

    def to_csv(
        self,
        path: str,
        properties: Optional[Sequence[Union[str, Type]]] = None,
        sort_by: Optional[str] = "smiles",
        **df_kwargs,
    ) -> str:
        """Write to CSV via `to_dataframe(...).to_csv(...)`.

        Defaults to `sort_by="smiles"` so the on-disk file is
        deterministic across runs (sets are unordered). Pass
        `sort_by=None` to skip sorting if order doesn't matter or if
        you've already sorted upstream.

        Returns the path that was written, so calls can chain into
        further file operations.
        """
        df = self.to_dataframe(properties=properties, sort_by=sort_by, **df_kwargs)
        df.to_csv(path, index=False)
        return path

    def to_jsonl(
        self,
        path: str,
        properties: Optional[Sequence[Union[str, Type]]] = None,
        sort_by: Optional[str] = "smiles",
        **df_kwargs,
    ) -> str:
        """Write to JSON Lines (one JSON object per row) via
        `to_dataframe(...).to_json(orient='records', lines=True)`.

        JSONL is the de-facto interchange format for ML training
        pipelines (HuggingFace datasets, JAX/Flax dataloaders) and is
        a more honest serialization of mixed-type data than CSV.
        """
        df = self.to_dataframe(properties=properties, sort_by=sort_by, **df_kwargs)
        df.to_json(path, orient="records", lines=True)
        return path

    # ----- repr ------------------------------------------------------

    def __repr__(self) -> str:
        n = len(self)
        if n == 0:
            return "GroupGraphSet(empty)"
        # Show first few SMILES so the repr is informative without
        # being unbounded for huge result sets.
        smis = []
        for g in self:
            try:
                smis.append(g.to_smiles())
            except Exception:
                smis.append("<unprintable>")
            if len(smis) >= 5:
                break
        more = f", ... ({n - 5} more)" if n > 5 else ""
        return f"GroupGraphSet({n} graphs: {smis}{more})"
