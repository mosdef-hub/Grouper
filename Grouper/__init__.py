# ruff: noqa: F401
from Grouper._Grouper import (
    Atom,
    AtomGraph,
    Group,
    GroupGraph,
    exhaustive_fragment,
)

from ._safeguards import (
    GrouperResourceLimitError,
    check_resource_limits,
    estimate_vcolg_lines,
)
from .fragmentation import fragment
from .generate import exhaustive_generate, random_generate, random_sample
from .group_graph_set import GroupGraphSet

__version__ = "0.0.5"


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


def _gg_visualize(self, figsize=(10, 5), structure_size=(400, 400)):
    """Render a side-by-side view of the molecule for this GroupGraph.

    Two panels:

      Left  — 2D chemical structure (atoms and bonds), via RDKit's
              `Draw.MolToImage` on the SMILES of this graph.
      Right — port-graph (groups as colored nodes, ports as edges),
              via `Grouper.visualization.visualize_graph.visualize`.

    The two views answer different questions: the left tells you
    *what molecule* this is in chemistry terms; the right tells you
    *how Grouper sees it* — which groups, which ports, which edges.
    Together they're the easiest way to sanity-check a generated
    structure in a notebook.

    Returns the matplotlib Figure so callers can `plt.show()`,
    `fig.savefig(...)`, or embed it.
    """
    import matplotlib.pyplot as plt
    from rdkit import Chem
    from rdkit.Chem import Draw

    from Grouper.visualization.visualize_graph import visualize

    fig, (ax_struct, ax_graph) = plt.subplots(1, 2, figsize=figsize)

    # Left: chemical structure. RDKit gives us a PIL image; matplotlib
    # imshow happily accepts it.
    mol = Chem.MolFromSmiles(self.to_smiles())
    if mol is not None:
        img = Draw.MolToImage(mol, size=structure_size)
        ax_struct.imshow(img)
    ax_struct.axis("off")
    ax_struct.set_title("Molecule")

    # Right: port graph. visualize() takes its own ax kwarg.
    visualize(self, ax=ax_graph)
    ax_graph.set_title("GroupGraph (ports)")

    fig.tight_layout()
    return fig


GroupGraph.to_3d = _gg_to_3d
GroupGraph.to_sdf = _gg_to_sdf
GroupGraph.to_mol = _gg_to_mol
GroupGraph.to_xyz = _gg_to_xyz
GroupGraph.to_pdb = _gg_to_pdb
GroupGraph.to_inchi = _gg_to_inchi
GroupGraph.to_inchi_key = _gg_to_inchi_key
GroupGraph.to_smarts = _gg_to_smarts
GroupGraph.visualize = _gg_visualize
