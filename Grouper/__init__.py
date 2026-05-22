# ruff: noqa: F401
from functools import wraps as _wraps

from Grouper._Grouper import (
    Atom,
    AtomGraph,
    Group,
    GroupGraph,
    exhaustive_fragment,
)
from Grouper._Grouper import exhaustive_generate as _cpp_exhaustive_generate
from Grouper._Grouper import random_generate as _cpp_random_generate

from ._safeguards import (
    GrouperResourceLimitError,
    check_resource_limits,
    estimate_vcolg_lines,
)


# Pre-flight wrappers around the two C++ entry points whose runtime
# scales combinatorially with `n_nodes`. The C++ side itself doesn't
# refuse oversized inputs, so a user pasting `exhaustive_generate(12,
# big_library)` into a notebook can lock up their machine for hours.
# These wrappers raise `GrouperResourceLimitError` with a clear
# message before any heavy work starts. Pass `confirm=True` (or
# raise `max_vcolg_lines`) to override.
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
    return _cpp_exhaustive_generate(
        n_nodes,
        node_defs,
        num_procs,
        vcolg_output_file,
        positive_constraints if positive_constraints is not None else {},
        negative_constraints if negative_constraints is not None else set(),
        config_path,
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
    return _cpp_random_generate(
        n_nodes,
        node_defs,
        num_graphs,
        num_procs,
        positive_constraints if positive_constraints is not None else {},
        negative_constraints if negative_constraints is not None else set(),
    )


from .fragmentation import (  # noqa: E402  (intentional: must come after the wrapper defs above)
    fragment,
)
