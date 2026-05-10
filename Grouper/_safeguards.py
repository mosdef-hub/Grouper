"""Pre-flight resource checks for `exhaustive_generate` / `random_generate`.

The combinatorial space `geng | vcolg` enumerates grows roughly as

    A001349(n) * n_colors^n

where A001349(n) is the number of connected simple graphs on n
vertices. That sequence is double-exponential — at n=10 it's ~12
million, at n=12 it's ~155 billion. A user who pastes
`exhaustive_generate(12, my_groups)` into a notebook expecting it to
"just run" can lock up their laptop, fill swap, and lose work.

This module provides cheap pre-flight estimates and configurable
thresholds so that a too-large run raises a clear error *before* any
heavy work starts, rather than silently consuming hours of CPU and
gigabytes of RAM. Pass `confirm=True` (or raise `max_vcolg_lines`) to
override when the configuration is intentional.
"""

from __future__ import annotations

import warnings


class GrouperResourceLimitError(Exception):
    """Raised when an enumeration's estimated workload exceeds the
    configured safety threshold. Pass `confirm=True` to the calling
    function (or raise the threshold via `max_vcolg_lines`) to
    override.
    """


# OEIS A001349: number of connected simple graphs on n nodes.
# Source: https://oeis.org/A001349
# These are the *un-colored* graphs that geng emits before vcolg
# multiplies them out by colorings.
_CONNECTED_GRAPHS_BY_N: dict[int, int] = {
    1:                  1,
    2:                  1,
    3:                  2,
    4:                  6,
    5:                 21,
    6:                112,
    7:                853,
    8:             11_117,
    9:            261_080,
    10:        11_716_571,
    # n >= 11 is treated as intractable below — the table stops here
    # because by then geng's output alone exceeds practical machine
    # capacity even before vcolg adds colorings.
}

_MAX_KNOWN_N: int = max(_CONNECTED_GRAPHS_BY_N)

# Default thresholds. Tuned against the perf measurements in the
# exhaustive-generate-fixes branch:
#   n=5, c=4 (8K vcolg lines):   ~3 seconds par5
#   n=6, c=4 (190K lines):       ~30 seconds par5
# Extrapolating linearly: ~5K lines/sec on a multi-core laptop.

# Soft limit: emit a UserWarning above this. ~3 minutes wall time on
# typical hardware; long enough to notice, short enough that you
# probably meant it.
_DEFAULT_WARN_VCOLG_LINES: int = 1_000_000

# Hard limit: raise GrouperResourceLimitError above this without
# `confirm=True`. ~3 hours wall time — the user almost certainly
# wants to know if the run will go this long.
_DEFAULT_MAX_VCOLG_LINES: int = 50_000_000

# Lines per second, used for human-readable wall-time estimates in
# error messages. Conservative single-thread-ish number; multi-core
# runs are faster.
_LINES_PER_SECOND_ESTIMATE: int = 5_000


def estimate_vcolg_lines(n_nodes: int, n_colors: int) -> int | None:
    """Conservative upper-bound estimate of vcolg's output line count
    for a given (n_nodes, n_colors) configuration.

    Returns `geng_count(n_nodes) * n_colors^n_nodes`, which is the
    upper bound *before* vcolg deduplicates symmetric colorings.
    Real output is typically 30–70% of this — overestimating is the
    correct side to err on for a safety check.

    Returns None when n_nodes is outside the lookup table (i.e.
    ``> _MAX_KNOWN_N``) — the caller treats those cases as
    intractable.
    """
    if n_nodes not in _CONNECTED_GRAPHS_BY_N:
        return None
    return _CONNECTED_GRAPHS_BY_N[n_nodes] * (n_colors ** n_nodes)


def _format_wall_time(estimated_lines: int) -> str:
    """Turn a line count into a rough human-readable wall-time string
    for use in warnings and error messages."""
    seconds = estimated_lines / _LINES_PER_SECOND_ESTIMATE
    if seconds < 60:
        return f"~{seconds:.0f} seconds"
    if seconds < 3600:
        return f"~{seconds / 60:.0f} minutes"
    if seconds < 86400:
        return f"~{seconds / 3600:.1f} hours"
    return f"~{seconds / 86400:.1f} days"


def check_resource_limits(
    n_nodes: int,
    n_colors: int,
    *,
    confirm: bool = False,
    max_vcolg_lines: int | None = None,
    warn_vcolg_lines: int | None = None,
) -> None:
    """Pre-flight check before exhaustive_generate / random_generate.

    Raises GrouperResourceLimitError if the configuration is too big
    AND `confirm` is False. Emits a UserWarning if the configuration
    is large but under the hard limit.

    `n_colors` is the number of distinct group definitions in
    `node_defs` (vcolg's `-m` argument). This is a misleading-
    looking name only because vcolg writes "colors" for what we'd
    call "group types"; we keep its terminology so the math lines up
    with vcolg's own output.

    Parameters not specified here use the module-level defaults.
    """
    if max_vcolg_lines is None:
        max_vcolg_lines = _DEFAULT_MAX_VCOLG_LINES
    if warn_vcolg_lines is None:
        warn_vcolg_lines = _DEFAULT_WARN_VCOLG_LINES

    # Single shared trailer. The strategies live in the tutorial
    # because they're a multi-paragraph topic; the messages here
    # surface the existence of the menu without trying to inline it.
    strategies_pointer = (
        "See docs/tutorials/runtime_safeguards.ipynb for strategies — "
        "multi-processing (num_procs=-1, on by default), random sampling "
        "(random_generate), constraint pruning, library reduction, "
        "file-based streaming (vcolg_output_file), and patterns for "
        "distributed compute."
    )

    if n_nodes > _MAX_KNOWN_N:
        if confirm:
            return
        raise GrouperResourceLimitError(
            f"n_nodes={n_nodes} is intractable for exhaustive enumeration: "
            f"the connected-graph count alone (OEIS A001349) is beyond "
            f"the reach of a single machine for n>10, even with all "
            f"cores in use. The largest safely-supported size is "
            f"n_nodes={_MAX_KNOWN_N}. "
            f"For larger n, use random_generate(..., num_graphs=N) to "
            f"sample, or pass confirm=True if you have a specialized "
            f"setup (precomputed vcolg output, distributed compute) "
            f"and accept the consequences. {strategies_pointer}"
        )

    estimated = estimate_vcolg_lines(n_nodes, n_colors)
    # `estimated` is not None here — n_nodes <= _MAX_KNOWN_N was
    # checked above.

    if estimated > max_vcolg_lines:
        if confirm:
            return
        raise GrouperResourceLimitError(
            f"This configuration (n_nodes={n_nodes}, "
            f"n_colors={n_colors}) would process up to "
            f"~{estimated:,} vcolg lines, exceeding the safety "
            f"threshold of {max_vcolg_lines:,}. "
            f"Estimated wall time: {_format_wall_time(estimated)} "
            f"(multi-processing is already on by default; this is "
            f"the work that overwhelms it). "
            f"To stay under the limit, reduce n_nodes by 1-2 or shrink "
            f"the group library. To proceed anyway, pass confirm=True "
            f"or raise the limit via max_vcolg_lines=<N>. "
            f"{strategies_pointer}"
        )

    if estimated > warn_vcolg_lines:
        warnings.warn(
            f"Large enumeration: ~{estimated:,} vcolg lines "
            f"(estimated wall time: {_format_wall_time(estimated)} on "
            f"multi-core hardware). Pass confirm=True to silence this "
            f"warning, or max_vcolg_lines=<N> to adjust the threshold. "
            f"{strategies_pointer}",
            stacklevel=3,
        )
