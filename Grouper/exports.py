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

For batch output (a whole `exhaustive_generate` result to one SDF):

    from Grouper.exports import write_sdf
    write_sdf(list_of_graphs, "out.sdf")
"""

from __future__ import annotations

import warnings
from typing import Iterable, Optional, Union

from rdkit import Chem
from rdkit.Chem import AllChem


class EmbedError(Exception):
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
    raise TypeError(
        f"expected SMILES str or GroupGraph; got {type(target).__name__}"
    )


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
            f"unknown force_field {force_field!r}; "
            f"choose 'mmff94', 'uff', or 'none'"
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
    target: Union[str, "GroupGraph"],
    path: Optional[str] = None,
    embed_3d: bool = True,
    name: Optional[str] = None,
    properties: Optional[dict] = None,
    **embed_kwargs,
) -> Optional[str]:
    """Render the molecule as an SDF block.

    `path=None` returns the SDF string; passing a path writes to disk
    and returns None.

    `embed_3d=True` (default) generates 3D coordinates via ETKDG +
    MMFF94. Set False for a 2D-only SDF (faster, smaller output —
    appropriate when downstream tools will re-embed themselves).

    `properties` is a dict of arbitrary key→value pairs that get
    written into the SDF as `> <KEY>` blocks — useful for piping
    Joback predictions or ML scores alongside the structure.
    """
    smiles = _smiles_of(target)
    if embed_3d:
        mol = to_3d_mol(smiles, **embed_kwargs)
    else:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise EmbedError(f"could not parse SMILES: {smiles!r}")

    if name is not None:
        mol.SetProp("_Name", name)
    if properties:
        for k, v in properties.items():
            mol.SetProp(str(k), str(v))

    # Route through SDWriter so the `> <KEY>` property tags get emitted.
    # MolToMolBlock drops them; writing to a single-mol SDF and stripping
    # gives us a string when no path is requested.
    if path is not None:
        with Chem.SDWriter(path) as writer:
            writer.write(mol)
        return None

    import io
    buf = io.StringIO()
    writer = Chem.SDWriter(buf)
    writer.write(mol)
    writer.close()
    return buf.getvalue()


def to_mol(
    target: Union[str, "GroupGraph"],
    path: Optional[str] = None,
    embed_3d: bool = True,
    **embed_kwargs,
) -> Optional[str]:
    """Render the molecule as a V2000 MOL block (single-molecule). For
    multi-molecule output use SDF (`to_sdf` or `write_sdf`)."""
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
    """Return the SMARTS pattern for the molecule. Useful when a
    generated structure becomes a query pattern for substructure
    searching elsewhere in a pipeline."""
    smiles = _smiles_of(target)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise EmbedError(f"could not parse SMILES: {smiles!r}")
    return Chem.MolToSmarts(mol)


def write_sdf(
    targets: Iterable[Union[str, "GroupGraph"]],
    path: str,
    embed_3d: bool = True,
    properties_fn: Optional[callable] = None,
    skip_failures: bool = False,
    **embed_kwargs,
) -> int:
    """Write a sequence of GroupGraphs (or SMILES) to a multi-molecule
    SDF file. Returns the number of molecules successfully written.

    `properties_fn(target)` is an optional callable that returns a
    dict of properties to attach to each molecule in the SDF — useful
    for embedding predicted properties alongside the structures::

        from Grouper.properties import JobackEstimate
        write_sdf(
            results, "screened.sdf",
            properties_fn=lambda g: dict(JobackEstimate.from_group_graph(g)),
        )

    `skip_failures=True` keeps going past molecules that fail to embed
    (e.g., strained rings ETKDG can't solve), logging a warning for
    each. Default is to stop and re-raise.
    """
    n_written = 0
    with Chem.SDWriter(path) as writer:
        for i, target in enumerate(targets):
            smiles = _smiles_of(target)
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

            mol.SetProp("_Name", f"mol_{i}")
            mol.SetProp("smiles", smiles)
            if properties_fn is not None:
                for k, v in properties_fn(target).items():
                    mol.SetProp(str(k), str(v))
            writer.write(mol)
            n_written += 1
    return n_written
