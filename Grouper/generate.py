# ruff: noqa: F401
"""Public generation entry points.

Three functions, all returning ``GroupGraphSet``:

- ``exhaustive_generate`` — full canonical-orbit enumeration via
  ``geng | vcolg``. Combinatorial, gated by ``_safeguards``.

- ``random_generate`` — uniform-over-orbits sampling on top of the
  ``vcolg`` list. Same safeguard gating because it shares the
  combinatorial pre-step.

- ``random_sample`` — direct random construction without
  ``geng``/``vcolg``. Linear in ``num_graphs``, works at sizes where
  the other two are infeasible. NOT uniform over orbits; the trade-off
  is documented in the C++ binding's docstring.

All three are thin Python wrappers around C++ entry points exposed
via the ``Grouper._Grouper`` pybind module. The wrappers add (a) the
``GrouperResourceLimitError`` pre-flight check for the vcolg-backed
functions and (b) the ``GroupGraphSet`` return coercion that exposes
batch convenience methods (``.to_dataframe()``, ``.filter()``, etc.).
"""

import warnings as _warnings

from Grouper._Grouper import exhaustive_generate as _cpp_exhaustive_generate
from Grouper._Grouper import random_generate as _cpp_random_generate
from Grouper._Grouper import random_sample as _cpp_random_sample

from ._safeguards import check_resource_limits
from .group_graph_set import GroupGraphSet


# Pre-flight wrappers around the two C++ entry points whose runtime
# scales combinatorially with ``n_nodes``. The C++ side itself doesn't
# refuse oversized inputs, so a user pasting
# ``exhaustive_generate(12, big_library)`` into a notebook can lock up
# their machine for hours. These wrappers raise
# ``GrouperResourceLimitError`` with a clear message before any heavy
# work starts. Pass ``confirm=True`` (or raise ``max_vcolg_lines``) to
# override.
#
# Each wrapper also coerces the C++ Python-set return value to a
# ``GroupGraphSet`` (subclass of set), exposing convenience methods
# from ``Grouper.group_graph_set`` — ``.to_dataframe()``,
# ``.to_smiles_list()``, ``.to_sdf()``, ``.filter()``. Existing
# patterns (iteration, ``len``, ``in``) work unchanged because
# GroupGraphSet IS a set.
#
# We deliberately do NOT wrap exhaustive_fragment. It returns a
# std::vector<GroupGraph> → Python list, not a set: fragment results
# are ordered by quality and may legitimately contain duplicates
# (different port assignments producing the same SMILES). Tests
# elsewhere rely on list semantics.
def exhaustive_generate(
    n_nodes,
    node_defs,
    num_procs=-1,
    vcolg_output_file="",
    positive_constraints=None,
    negative_constraints=None,
    config_path="",
    *,
    confirm=False,
    max_vcolg_lines=None,
):
    check_resource_limits(
        n_nodes,
        len(node_defs),
        confirm=confirm,
        max_vcolg_lines=max_vcolg_lines,
    )
    return GroupGraphSet(
        _cpp_exhaustive_generate(
            n_nodes,
            node_defs,
            num_procs,
            vcolg_output_file,
            positive_constraints if positive_constraints is not None else {},
            negative_constraints if negative_constraints is not None else set(),
            config_path,
        )
    )


def random_generate(
    n_nodes,
    node_defs,
    num_graphs,
    num_procs=-1,
    positive_constraints=None,
    negative_constraints=None,
    *,
    confirm=False,
    max_vcolg_lines=None,
):
    check_resource_limits(
        n_nodes,
        len(node_defs),
        confirm=confirm,
        max_vcolg_lines=max_vcolg_lines,
    )
    return GroupGraphSet(
        _cpp_random_generate(
            n_nodes,
            node_defs,
            num_graphs,
            num_procs,
            positive_constraints if positive_constraints is not None else {},
            negative_constraints if negative_constraints is not None else set(),
        )
    )


# Direct random sampler — no vcolg pre-step, so no safeguard gate.
# Detailed semantics and trade-offs (not uniform over orbits, scales
# linearly with num_graphs) live in the C++ binding's docstring.
def random_sample(
    n_nodes,
    node_defs,
    num_graphs,
    num_procs=-1,
    positive_constraints=None,
    negative_constraints=None,
    *,
    extra_edge_prob=0.10,
    color_strategy="stratified",
    max_attempts=0,
    seed=-1,
    show_progress=True,
):
    raw = _cpp_random_sample(
        n_nodes,
        node_defs,
        num_graphs,
        num_procs,
        positive_constraints if positive_constraints is not None else {},
        negative_constraints if negative_constraints is not None else set(),
        extra_edge_prob,
        color_strategy,
        max_attempts,
        seed,
        show_progress,
    )
    if len(raw) < num_graphs:
        effective_max = max_attempts if max_attempts > 0 else 20 * num_graphs
        _warnings.warn(
            f"random_sample exhausted max_attempts={effective_max} after "
            f"finding {len(raw)}/{num_graphs} unique graphs. "
            f"Consider relaxing constraints or raising max_attempts.",
            stacklevel=2,
        )
    return GroupGraphSet(raw)
