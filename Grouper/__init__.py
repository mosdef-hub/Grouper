# ruff: noqa: F401
from functools import wraps as _wraps

from Grouper._Grouper import (
    Atom,
    AtomGraph,
    Group,
    GroupGraph,
)
from Grouper._Grouper import exhaustive_fragment as _cpp_exhaustive_fragment
from Grouper._Grouper import exhaustive_generate as _cpp_exhaustive_generate
from Grouper._Grouper import random_generate as _cpp_random_generate

from ._safeguards import (
    GrouperResourceLimitError,
    check_resource_limits,
    estimate_vcolg_lines,
)
from .results import GroupGraphSet


# Pre-flight wrappers around the two C++ entry points whose runtime
# scales combinatorially with `n_nodes`. The C++ side itself doesn't
# refuse oversized inputs, so a user pasting `exhaustive_generate(12,
# big_library)` into a notebook can lock up their machine for hours.
# These wrappers raise `GrouperResourceLimitError` with a clear
# message before any heavy work starts. Pass `confirm=True` (or
# raise `max_vcolg_lines`) to override.
#
# Each wrapper also coerces the C++ Python-set return value to a
# `GroupGraphSet` (subclass of set), exposing the convenience methods
# from Grouper.results — `.to_dataframe()`, `.to_smiles_list()`,
# `.to_sdf()`, `.filter()`. Existing patterns (iteration, `len`, `in`)
# work unchanged because GroupGraphSet IS a set.
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


@_wraps(_cpp_exhaustive_fragment)
def exhaustive_fragment(*args, **kwargs):
    return GroupGraphSet(_cpp_exhaustive_fragment(*args, **kwargs))


from .fragmentation import (  # noqa: E402  (intentional: must come after the wrapper defs above)
    fragment,
)


# Add Pythonic export shortcuts to the C++ GroupGraph class. We do this
# via attribute assignment rather than C++ binding changes so adding
# new format methods doesn't require a rebuild. Each shortcut is a
# thin wrapper around a function in Grouper.exports.
def _gg_to_3d(self, **kw):
    """Generate 3D coordinates via ETKDG + force-field minimization.
    Returns an RDKit Mol with explicit hydrogens. See
    Grouper.exports.to_3d_mol for parameters."""
    from Grouper.exports import to_3d_mol

    return to_3d_mol(self, **kw)


def _gg_to_sdf(self, path=None, **kw):
    """Render as SDF. Pass `path` to write to disk; otherwise returns
    the SDF as a string. See Grouper.exports.to_sdf for parameters."""
    from Grouper.exports import to_sdf

    return to_sdf(self, path, **kw)


def _gg_to_mol(self, path=None, **kw):
    """Render as a V2000 MOL block (single-molecule)."""
    from Grouper.exports import to_mol

    return to_mol(self, path, **kw)


def _gg_to_xyz(self, path=None, **kw):
    """Render as XYZ (always 3D)."""
    from Grouper.exports import to_xyz

    return to_xyz(self, path, **kw)


def _gg_to_pdb(self, path=None, **kw):
    """Render as PDB (always 3D)."""
    from Grouper.exports import to_pdb

    return to_pdb(self, path, **kw)


def _gg_to_inchi(self):
    """Return the IUPAC InChI string (cross-tool deterministic
    identifier)."""
    from Grouper.exports import to_inchi
    return to_inchi(self)


def _gg_to_inchi_key(self):
    """Return the 27-character InChIKey — primary-key form used by
    most chemistry databases (PubChem, ChEMBL, NIST)."""
    from Grouper.exports import to_inchi_key
    return to_inchi_key(self)


def _gg_to_smarts(self):
    """Return the SMARTS pattern for this molecule. Useful as a query
    pattern in substructure searches downstream."""
    from Grouper.exports import to_smarts
    return to_smarts(self)


def _gg_visualize(self, **kw):
    """Render the port-graph (groups as colored nodes, port
    connections as edges) via matplotlib. Distinct from a 2D
    chemical-structure drawing — for that, use `to_smiles()` and
    feed it to RDKit's Draw module.

    Wraps `Grouper.visualization.visualize_graph.visualize`."""
    from Grouper.visualization.visualize_graph import visualize
    return visualize(self, **kw)


GroupGraph.to_3d = _gg_to_3d
GroupGraph.to_sdf = _gg_to_sdf
GroupGraph.to_mol = _gg_to_mol
GroupGraph.to_xyz = _gg_to_xyz
GroupGraph.to_pdb = _gg_to_pdb
GroupGraph.to_inchi = _gg_to_inchi
GroupGraph.to_inchi_key = _gg_to_inchi_key
GroupGraph.to_smarts = _gg_to_smarts
GroupGraph.visualize = _gg_visualize
