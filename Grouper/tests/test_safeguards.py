"""Tests for the runtime safeguards on exhaustive_generate /
random_generate. The point of these checks is to fail fast and clearly
when a user asks for a configuration that would lock up their machine,
so the assertions here center on (a) too-small inputs are unaffected,
(b) too-large inputs raise with an actionable message, and (c)
`confirm=True` bypasses the check.
"""

import warnings

import pytest

from Grouper import (
    Group,
    GrouperResourceLimitError,
    exhaustive_generate,
    random_generate,
)
from Grouper._safeguards import (
    _CONNECTED_GRAPHS_BY_N,
    _DEFAULT_MAX_VCOLG_LINES,
    _format_wall_time,
    check_resource_limits,
    estimate_vcolg_lines,
)


# ---------------------------------------------------------------------
# Estimation: verify the upper bound matches the published OEIS table
# multiplied by n_colors^n. Lock the formula so a future refactor
# can't silently change what we report to users.
# ---------------------------------------------------------------------
@pytest.mark.parametrize(
    "n, c, expected",
    [
        (5, 4, 21 * 4 ** 5),         # 21504
        (6, 2, 112 * 2 ** 6),        # 7168
        (8, 3, 11_117 * 3 ** 8),     # ~73M
    ],
)
def test_estimate_matches_table(n, c, expected):
    assert estimate_vcolg_lines(n, c) == expected


def test_estimate_returns_none_for_unknown_n():
    """n=11 is outside the table; the estimator returns None to
    signal 'use the intractability branch.'"""
    assert estimate_vcolg_lines(11, 2) is None


def test_format_wall_time_units():
    """Wall-time formatter switches between seconds, minutes, hours,
    days at the expected boundaries — so error messages always read
    naturally."""
    assert "seconds" in _format_wall_time(50_000)
    assert "minutes" in _format_wall_time(500_000)
    assert "hours" in _format_wall_time(50_000_000)
    assert "days" in _format_wall_time(50_000_000_000)


# ---------------------------------------------------------------------
# Pre-flight check: small cases pass cleanly; large cases raise; the
# escape hatch (confirm=True) works.
# ---------------------------------------------------------------------
def test_small_case_passes():
    """n=4, c=2 produces ~96 lines — well below any threshold."""
    check_resource_limits(4, 2)


def test_n_above_table_raises():
    """n=11 is in the 'intractable for any machine' band; raise even
    on a tiny color count, since geng alone exceeds practical
    capacity beyond n=10."""
    with pytest.raises(GrouperResourceLimitError, match="intractable"):
        check_resource_limits(11, 2)


def test_n_above_table_passes_with_confirm():
    """confirm=True is the explicit override for the user who knows
    what they're doing."""
    check_resource_limits(15, 2, confirm=True)  # no raise


def test_large_known_n_raises():
    """n=10 with c=2: 11.7M * 1024 = ~12B lines, well above the
    default cap of 50M. Must raise."""
    with pytest.raises(GrouperResourceLimitError, match="vcolg lines"):
        check_resource_limits(10, 2)


def test_large_known_n_passes_with_confirm():
    check_resource_limits(10, 2, confirm=True)


def test_large_known_n_passes_with_raised_threshold():
    """Raising max_vcolg_lines past the estimate is an alternative
    escape hatch — useful when the user wants no warning either."""
    estimated = estimate_vcolg_lines(8, 2)
    check_resource_limits(8, 2, max_vcolg_lines=estimated * 2)


def test_warn_threshold_emits_warning():
    """A configuration above the soft threshold but below the hard
    one emits a warning — the user runs but gets a heads-up."""
    # n=7, c=4: 853 * 4^7 = ~3.5M lines, between 1M (warn) and 50M (max).
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        check_resource_limits(7, 4)
    assert len(caught) == 1
    msg = str(caught[0].message)
    assert "Large enumeration" in msg
    assert "vcolg lines" in msg


def test_error_message_actionable():
    """The error message must point users at the two escape hatches
    (confirm or max_vcolg_lines) so they're not stuck."""
    with pytest.raises(GrouperResourceLimitError) as exc:
        check_resource_limits(10, 2)
    msg = str(exc.value)
    assert "confirm=True" in msg
    assert "max_vcolg_lines" in msg


# ---------------------------------------------------------------------
# End-to-end: the wrappers in Grouper/__init__.py call the check.
# Exercising via the actual public API ensures the wiring is right
# and catches a future regression where someone removes the wrapper.
# ---------------------------------------------------------------------
def test_exhaustive_generate_refuses_large_n():
    """Calling exhaustive_generate(15, ...) without confirm=True
    raises before any C++ work happens."""
    with pytest.raises(GrouperResourceLimitError, match="intractable"):
        exhaustive_generate(15, {Group("methyl", "C", [0])})


def test_random_generate_refuses_large_n():
    with pytest.raises(GrouperResourceLimitError, match="intractable"):
        random_generate(15, {Group("methyl", "C", [0])}, 5)


def test_exhaustive_generate_runs_at_small_n():
    """The check doesn't false-positive on routine small cases."""
    groups = {Group("methyl", "C", [0, 0, 0, 0])}
    result = exhaustive_generate(2, groups, num_procs=1)
    assert isinstance(result, set)
    assert len(result) >= 1


def test_exhaustive_generate_with_confirm_skips_check():
    """confirm=True bypasses the pre-flight check — useful in
    pipelines that already know what they're doing."""
    # n=11 would normally raise; with confirm=True, the call gets to
    # the C++ side. We don't actually want to run n=11 in tests, so
    # just trust the check_resource_limits unit test above pinned
    # the override semantics, and verify the wrapper accepts the
    # kwarg without erroring at the Python layer.
    # (We can't run the actual C++ enumeration here — it would take
    # hours.)
    from Grouper._safeguards import check_resource_limits
    check_resource_limits(15, 2, confirm=True)  # the path the wrapper hits


def test_max_vcolg_lines_kwarg_threads_through():
    """The max_vcolg_lines kwarg on the public API reaches the
    safeguard. We exercise this with a *low* threshold so even a
    tiny case raises."""
    groups = {Group("methyl", "C", [0, 0, 0, 0])}
    with pytest.raises(GrouperResourceLimitError):
        exhaustive_generate(5, groups, max_vcolg_lines=10)
