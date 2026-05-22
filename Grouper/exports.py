"""3D conformer generation and chemical-file-format export for GroupGraphs.

The functions here turn a GroupGraph (or a SMILES string) into the file
formats that downstream tools — molecular dynamics packages (GROMACS,
LAMMPS, OpenMM), docking programs (AutoDock Vina, Glide), and
visualizers (PyMOL, VMD) — actually consume:

    SDF    multi-molecule with optional 3D coords + arbitrary properties
    MOL    single-molecule V2000 mol block, like SDF without the $$$$
    XYZ    bare element + 3D coords; used by quantum chemistry codes
    PDB    Protein Data Bank format; used by visualizers and docking

3D coordinates are generated via RDKit's ETKDG distance-geometry
embedding followed by force-field minimization. We default to MMFF94
(more accurate where parameters exist) and fall back to UFF (covers
all atom types in Joback's element scope) if MMFF fails.

Methods on GroupGraph (added via monkey-patch in Grouper/__init__.py):

    gG.to_3d(seed=42, force_field="mmff94")  -> rdkit.Chem.Mol
    gG.to_sdf(path=None, embed_3d=True)      -> str | None
    gG.to_mol(path=None, embed_3d=True)      -> str | None
    gG.to_xyz(path=None)                     -> str | None  (always 3D)
    gG.to_pdb(path=None)                     -> str | None  (always 3D)

The same `to_sdf` function also handles a whole library — pass any
iterable of graphs/SMILES instead of a single target:

    from Grouper.exports import to_sdf
    to_sdf(results, path="library.sdf")            # multi-mol SDF
    to_sdf(["CCO", "CC"], path="out.sdf")          # raw SMILES list
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Optional, Union

from rdkit import Chem
from rdkit.Chem import AllChem

if TYPE_CHECKING:
    # Type-only import: GroupGraph comes from a C++ binding and importing
    # it at runtime would create a circular dependency through
    # Grouper.__init__.py. The string-forward-reference in the function
    # signatures below points here, satisfying static analyzers.
    from Grouper import GroupGraph  # noqa: F401


class GrouperError(Exception):
    """Base class for errors raised by Grouper. Catch this to
    distinguish failures originating in Grouper from those bubbling
    up from dependencies (RDKit, NumPy, etc.).

    Future Grouper-raised exceptions should inherit from this so a
    single ``except GrouperError`` handler covers them all.
    """


class EmbedError(GrouperError):
    """Raised when 3D coordinate generation fails for a molecule.

    Common causes: highly strained ring systems where ETKDG can't find
    a feasible geometry, or molecules with element/valence combinations
    RDKit's force fields don't parameterize. The message includes the
    SMILES so the caller can reproduce.
    """


def _smiles_of(target: Union[str, "GroupGraph"]) -> str:
    """Coerce a SMILES str or a GroupGraph (anything with `to_smiles()`)
    to a SMILES string. Duck-typed to avoid importing GroupGraph here
    (would create a circular import via Grouper/__init__.py)."""
    if isinstance(target, str):
        return target
    if hasattr(target, "to_smiles") and callable(target.to_smiles):
        return target.to_smiles()
    raise TypeError(f"expected SMILES str or GroupGraph; got {type(target).__name__}")


def _mol_with_explicit_h(smiles: str) -> Chem.Mol:
    """Parse a SMILES and add explicit hydrogens. Required before 3D
    embedding — implicit Hs aren't placed by ETKDG, so the resulting
    geometry would be missing them."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise EmbedError(f"could not parse SMILES: {smiles!r}")
    return Chem.AddHs(mol)


def to_3d_mol(
    target: Union[str, "GroupGraph"],
    seed: int = 42,
    force_field: str = "mmff94",
    max_iters: int = 200,
) -> Chem.Mol:
    """Embed and minimize 3D coordinates for `target`. Returns an RDKit
    Mol with explicit hydrogens and one conformer.

    `force_field`:
        - "mmff94"  — Merck molecular force field (default; more
                      accurate where parameters exist).
        - "uff"     — Universal Force Field (covers everything but
                      less accurate for typical organics).
        - "none"    — embed only, skip minimization.

    Falls back to UFF automatically if MMFF94 fails (e.g., the molecule
    has an atom type MMFF doesn't parameterize) — typical organics in
    Joback's element scope are well-covered by either.

    Raises `EmbedError` if ETKDG distance-geometry can't find a
    feasible geometry. Common with strained polycyclics or unusual
    bond patterns; the SMILES is included in the message.
    """
    smiles = _smiles_of(target)
    mol = _mol_with_explicit_h(smiles)

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    rc = AllChem.EmbedMolecule(mol, params)
    if rc < 0:
        # Try once more with the legacy ETKDG (sometimes succeeds where
        # v3 fails on small or strained molecules).
        rc = AllChem.EmbedMolecule(mol, randomSeed=seed)
        if rc < 0:
            raise EmbedError(
                f"3D embedding failed for SMILES: {smiles!r}. "
                f"Common cause: strained ring system or unusual valence."
            )

    if force_field == "mmff94":
        try:
            status = AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
            if status == -1:
                # MMFF couldn't find parameters for some atom; fall back.
                AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
        except Exception:
            AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    elif force_field == "uff":
        AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    elif force_field == "none":
        pass
    else:
        raise ValueError(
            f"unknown force_field {force_field!r}; choose 'mmff94', 'uff', or 'none'"
        )

    return mol


def _write_or_return(text: str, path: Optional[str]) -> Optional[str]:
    """Either write to disk or return as string, depending on whether
    a path was passed. Lets every export function support both
    workflows from one signature."""
    if path is None:
        return text
    with open(path, "w") as f:
        f.write(text)
    return None


def to_sdf(
    target,
    path: Optional[str] = None,
    embed_3d: bool = True,
    name: Optional[str] = None,
    properties=None,
    skip_failures: bool = False,
    **embed_kwargs,
):
    """Render one molecule or a library as SDF — single function for
    both cases, since SDF natively supports one-or-many molecules per
    file.

    `target` may be:
      - a single `GroupGraph` or SMILES string, or
      - any iterable of those (e.g. a `GroupGraphSet`, a list of
        SMILES, a generator).

    `path`:
      - `None` (default) — return the SDF text as a string.
      - a filesystem path — write to disk; return the number of
        molecules successfully written.

    `properties`:
      - `None` — no extra property tags.
      - `dict` — apply the same key→value pairs to every molecule.
      - `callable(target) -> dict` — invoked per molecule; useful for
        attaching per-molecule predictions, e.g.::

            to_sdf(
                results, "library.sdf",
                properties=lambda g: dict(JobackEstimate.from_group_graph(g)),
            )

    `name` is set as the SDF `_Name` field for the single-target case.
    For iterables, names default to ``mol_<index>`` and `name` is
    ignored.

    `skip_failures=True` continues past molecules that fail to embed
    (e.g., strained rings ETKDG can't solve) or parse, emitting a
    warning for each. Default is to stop and re-raise.

    `embed_3d=True` (default) runs ETKDG + MMFF94 for 3D coords. Set
    False for a 2D-only SDF.
    """
    single = isinstance(target, str) or hasattr(target, "to_smiles")
    targets = [target] if single else target

    if properties is None:
        props_for = None
    elif callable(properties):
        props_for = properties
    elif isinstance(properties, dict):
        _frozen = dict(properties)

        def props_for(_t):
            return _frozen
    else:
        raise TypeError(
            "properties must be None, a dict, or a callable; "
            f"got {type(properties).__name__}"
        )

    import io

    buf = io.StringIO() if path is None else None
    writer = Chem.SDWriter(buf if path is None else path)
    n_written = 0
    try:
        for i, t in enumerate(targets):
            smiles = _smiles_of(t)
            try:
                if embed_3d:
                    mol = to_3d_mol(smiles, **embed_kwargs)
                else:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        raise EmbedError(f"could not parse SMILES: {smiles!r}")
            except EmbedError as e:
                if skip_failures:
                    warnings.warn(f"skipping #{i} ({smiles!r}): {e}")
                    continue
                raise

            if single and name is not None:
                mol.SetProp("_Name", name)
            elif not single:
                mol.SetProp("_Name", f"mol_{i}")
                mol.SetProp("smiles", smiles)
            if props_for is not None:
                for k, v in props_for(t).items():
                    mol.SetProp(str(k), str(v))
            writer.write(mol)
            n_written += 1
    finally:
        writer.close()

    if path is None:
        return buf.getvalue()
    return n_written


def to_mol(
    target: Union[str, "GroupGraph"],
    path: Optional[str] = None,
    embed_3d: bool = True,
    **embed_kwargs,
) -> Optional[str]:
    """Render the molecule as a V2000 MOL block (single-molecule). For
    multi-molecule output use `to_sdf` with an iterable target."""
    smiles = _smiles_of(target)
    if embed_3d:
        mol = to_3d_mol(smiles, **embed_kwargs)
    else:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise EmbedError(f"could not parse SMILES: {smiles!r}")
    return _write_or_return(Chem.MolToMolBlock(mol), path)


def to_xyz(
    target: Union[str, "GroupGraph"],
    path: Optional[str] = None,
    **embed_kwargs,
) -> Optional[str]:
    """Render the molecule as XYZ. Always 3D — XYZ doesn't have a 2D
    mode. Used as input to quantum-chemistry codes (Gaussian, ORCA,
    NWChem, Psi4) and many MD setup tools."""
    mol = to_3d_mol(_smiles_of(target), **embed_kwargs)
    return _write_or_return(Chem.MolToXYZBlock(mol), path)


def to_pdb(
    target: Union[str, "GroupGraph"],
    path: Optional[str] = None,
    **embed_kwargs,
) -> Optional[str]:
    """Render the molecule as PDB. Always 3D. Used by structural
    biology / visualizer tools (PyMOL, VMD, ChimeraX) and many
    docking programs."""
    mol = to_3d_mol(_smiles_of(target), **embed_kwargs)
    return _write_or_return(Chem.MolToPDBBlock(mol), path)


# ---------------------------------------------------------------------
# Identifier formats: InChI, InChIKey, SMARTS. These don't need 3D
# coordinates — they're 1D string identifiers used for
# database lookup (PubChem, ChEMBL, NIST), substructure querying, and
# stable cross-tool referencing.
# ---------------------------------------------------------------------
def to_inchi(target: Union[str, "GroupGraph"]) -> str:
    """Return the IUPAC InChI string. InChI is a layered, deterministic
    identifier that's standard across cheminformatics databases — two
    drawings of the same molecule produce the same InChI even when
    their canonical SMILES differ across implementations.

    Used to look molecules up in PubChem, ChEMBL, NIST WebBook, etc.
    For shorter primary-key use, prefer `to_inchi_key`.
    """
    smiles = _smiles_of(target)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise EmbedError(f"could not parse SMILES: {smiles!r}")
    return Chem.MolToInchi(mol)


def to_inchi_key(target: Union[str, "GroupGraph"]) -> str:
    """Return the 27-character InChIKey hash of the InChI. This is the
    primary-key form most chemistry databases use — fixed length,
    URL-safe, and deterministic across implementations.

    Format: `XXXXXXXXXXXXXX-YYYYYYYYYY-Z` where the dashes separate
    the structural-skeleton hash, the stereochemistry/isotope hash,
    and a final disambiguation block.
    """
    smiles = _smiles_of(target)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise EmbedError(f"could not parse SMILES: {smiles!r}")
    return Chem.MolToInchiKey(mol)


def to_smarts(target: Union[str, "GroupGraph"]) -> str:
    """Return the SMARTS pattern for the molecule, as produced by
    RDKit's ``Chem.MolToSmarts``. Useful when a generated structure
    becomes a query pattern for substructure searching elsewhere in
    a pipeline.

    Note: a molecule has infinitely many valid SMARTS renderings;
    this returns RDKit's canonical form, which may differ from the
    output of other toolkits (OpenEye, Daylight, etc.). For exact
    string equality across tools, normalize through a shared
    canonicalizer.
    """
    smiles = _smiles_of(target)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise EmbedError(f"could not parse SMILES: {smiles!r}")
    return Chem.MolToSmarts(mol)
