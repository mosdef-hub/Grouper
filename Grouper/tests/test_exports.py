"""Tests for the chemical-format export module.

The fragmentation/Joback tests check that we *interpret* SMILES
correctly. These tests check that we *emit* SMILES (and friends) in
the formats downstream tools consume — round-tripping through RDKit
to confirm the output parses back to a structurally-equivalent
molecule.
"""

import pytest

from rdkit import Chem

from Grouper import GroupGraph
from Grouper.exports import (
    EmbedError,
    to_3d_mol,
    to_mol,
    to_pdb,
    to_sdf,
    to_xyz,
    write_sdf,
)


# ---------------------------------------------------------------------
# 3D embedding
# ---------------------------------------------------------------------
def test_to_3d_mol_returns_mol_with_one_conformer():
    """ETKDG + force-field minimization must return a mol with 3D
    coordinates. The exact coordinates depend on RDKit version, so we
    only check that a conformer exists and is non-zero."""
    mol = to_3d_mol("CCO")
    assert mol.GetNumConformers() == 1
    conf = mol.GetConformer()
    # n-ethanol has 9 atoms after AddHs (2 heavy C + 1 O + 6 H)
    assert mol.GetNumAtoms() == 9
    # Coords should not all be at origin
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    assert not all(c.x == 0 and c.y == 0 and c.z == 0 for c in coords)


def test_to_3d_mol_seed_is_reproducible():
    """Same seed -> same coordinates (deterministic embedding)."""
    a = to_3d_mol("CCO", seed=42).GetConformer()
    b = to_3d_mol("CCO", seed=42).GetConformer()
    for i in range(9):
        assert abs(a.GetAtomPosition(i).x - b.GetAtomPosition(i).x) < 1e-9


def test_to_3d_mol_force_field_options():
    """All three force_field choices should produce a valid mol."""
    for ff in ["mmff94", "uff", "none"]:
        m = to_3d_mol("CCCC", force_field=ff)
        assert m.GetNumConformers() == 1


def test_to_3d_mol_invalid_smiles_raises():
    with pytest.raises(EmbedError, match="could not parse"):
        to_3d_mol("not a smiles")


# ---------------------------------------------------------------------
# Format round-trips: emit, then parse back, then check the structure
# matches. RDKit's canonical SMILES is the equivalence relation.
# ---------------------------------------------------------------------
@pytest.mark.parametrize("smiles", ["CCO", "CC(=O)C", "c1ccccc1", "CC(C)C"])
def test_sdf_round_trip(smiles):
    """SDF output must parse back to a structurally-equivalent mol."""
    sdf = to_sdf(smiles)
    assert sdf.endswith("$$$$\n"), "SDF must end with $$$$ terminator"
    mol = Chem.MolFromMolBlock(sdf)  # MolBlock parser handles SDF body
    assert mol is not None
    canon_in = Chem.CanonSmiles(smiles)
    canon_out = Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(mol)))
    assert canon_in == canon_out


@pytest.mark.parametrize("smiles", ["CCO", "CC(=O)C", "c1ccccc1"])
def test_mol_round_trip(smiles):
    """V2000 MOL block must parse back to the same molecule."""
    text = to_mol(smiles)
    mol = Chem.MolFromMolBlock(text)
    assert mol is not None
    assert Chem.CanonSmiles(smiles) == Chem.CanonSmiles(
        Chem.MolToSmiles(Chem.RemoveHs(mol))
    )


@pytest.mark.parametrize("smiles", ["CCO", "CC(=O)C", "c1ccccc1"])
def test_xyz_round_trip(smiles):
    """XYZ has no bond information — round-trip checks the atom count
    and elements rather than the SMILES."""
    xyz = to_xyz(smiles)
    lines = xyz.strip().split("\n")
    n_atoms_declared = int(lines[0])
    # Lines after the header + comment line are atom records:
    # `<symbol> <x> <y> <z>`. RDKit's XYZ writer emits exactly
    # n_atoms_declared atom lines.
    atom_lines = lines[2:]
    assert len(atom_lines) == n_atoms_declared
    # Every atom line must have an element symbol + 3 floats.
    for line in atom_lines:
        parts = line.split()
        assert len(parts) == 4
        sym = parts[0]
        assert sym.isalpha()
        for v in parts[1:]:
            float(v)  # raises if not a valid float


def test_pdb_contains_atom_records():
    """PDB output must include ATOM / HETATM records."""
    pdb = to_pdb("CCO")
    has_atom = any(line.startswith(("ATOM", "HETATM")) for line in pdb.split("\n"))
    assert has_atom


# ---------------------------------------------------------------------
# Disk output paths
# ---------------------------------------------------------------------
def test_to_sdf_writes_file(tmp_path):
    """Path argument writes to disk and returns None; file content
    matches the in-memory string output."""
    out = tmp_path / "out.sdf"
    result = to_sdf("CCO", path=str(out))
    assert result is None
    text_on_disk = out.read_text()
    text_in_memory = to_sdf("CCO")
    assert text_on_disk == text_in_memory


def test_to_xyz_writes_file(tmp_path):
    out = tmp_path / "out.xyz"
    to_xyz("CCO", path=str(out))
    text = out.read_text()
    assert text.startswith("9\n")  # 9 atoms after AddHs


# ---------------------------------------------------------------------
# Property attachment in SDF
# ---------------------------------------------------------------------
def test_sdf_property_block():
    """Properties dict must be embedded as `> <KEY>` blocks."""
    sdf = to_sdf("CCO", name="ethanol", properties={"Tb": 337.5, "source": "joback"})
    assert "<Tb>" in sdf
    assert "337.5" in sdf
    assert "<source>" in sdf
    assert "joback" in sdf


# ---------------------------------------------------------------------
# Multi-molecule SDF
# ---------------------------------------------------------------------
def test_write_sdf_multi_molecule(tmp_path):
    """write_sdf must concatenate molecules with $$$$ separators that
    RDKit's SDMolSupplier can re-iterate."""
    out = tmp_path / "multi.sdf"
    n = write_sdf(["CCO", "CC(=O)C", "c1ccccc1"], path=str(out))
    assert n == 3

    supplier = Chem.SDMolSupplier(str(out), removeHs=False)
    mols = [m for m in supplier if m is not None]
    assert len(mols) == 3
    canons_in = {Chem.CanonSmiles(s) for s in ["CCO", "CC(=O)C", "c1ccccc1"]}
    canons_out = {Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(m))) for m in mols}
    assert canons_in == canons_out


def test_write_sdf_with_properties_fn(tmp_path):
    """A properties_fn callback should attach per-molecule property
    blocks. Verify the round-tripped molecule carries them."""
    out = tmp_path / "with_props.sdf"

    def props_for(smi):
        return {"length": len(smi), "tag": "test"}

    write_sdf(["CCO", "CCCCC"], path=str(out), properties_fn=props_for)

    supplier = Chem.SDMolSupplier(str(out))
    mols = [m for m in supplier if m is not None]
    assert mols[0].GetProp("length") == "3"
    assert mols[0].GetProp("tag") == "test"
    assert mols[1].GetProp("length") == "5"


def test_write_sdf_skip_failures(tmp_path):
    """skip_failures=True should warn and continue past unembed-able
    molecules. Use an obviously bad SMILES alongside a good one."""
    out = tmp_path / "skipped.sdf"
    with pytest.warns(UserWarning, match="skipping"):
        n = write_sdf(
            ["CCO", "not a smiles", "CC"],
            path=str(out),
            skip_failures=True,
        )
    assert n == 2


# ---------------------------------------------------------------------
# GroupGraph monkey-patched methods
# ---------------------------------------------------------------------
def _ethanol_groupgraph():
    """Build n-ethanol from atomic groups for the GG-method tests."""
    gg = GroupGraph()
    gg.add_node("methyl",   "C", [0])
    gg.add_node("methylene","C", [0, 0])
    gg.add_node("hydroxyl", "O", [0])
    gg.add_edge((0, 0), (1, 0))
    gg.add_edge((1, 1), (2, 0))
    return gg


def test_groupgraph_to_3d_method():
    """gG.to_3d() returns an RDKit mol with the expected atom count."""
    gg = _ethanol_groupgraph()
    mol = gg.to_3d()
    # 2 C + 1 O = 3 heavy atoms; AddHs gives 6 implicit H -> 9 total.
    assert mol.GetNumAtoms() == 9
    assert mol.GetNumConformers() == 1


def test_groupgraph_to_sdf_method_round_trip():
    """gG.to_sdf() output round-trips to the same canonical SMILES."""
    gg = _ethanol_groupgraph()
    sdf = gg.to_sdf()
    mol = Chem.MolFromMolBlock(sdf)
    assert Chem.CanonSmiles("CCO") == Chem.CanonSmiles(
        Chem.MolToSmiles(Chem.RemoveHs(mol))
    )


def test_groupgraph_format_methods_exist():
    """Every promised export method is bound on GroupGraph and callable."""
    gg = _ethanol_groupgraph()
    for method_name in (
        "to_3d", "to_sdf", "to_mol", "to_xyz", "to_pdb",
        "to_inchi", "to_inchi_key", "to_smarts", "visualize",
    ):
        assert hasattr(gg, method_name), f"{method_name} missing on GroupGraph"
        assert callable(getattr(gg, method_name))


# ---------------------------------------------------------------------
# Identifier formats (InChI, InChIKey, SMARTS) — single-line strings,
# no 3D coords needed.
# ---------------------------------------------------------------------
@pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1", "CC(=O)O"])
def test_to_inchi_round_trip(smiles):
    """InChI must round-trip via Chem.MolFromInchi back to the same
    canonical SMILES."""
    from Grouper.exports import to_inchi
    inchi = to_inchi(smiles)
    assert inchi.startswith("InChI=")
    mol = Chem.MolFromInchi(inchi)
    assert mol is not None
    assert Chem.CanonSmiles(smiles) == Chem.CanonSmiles(Chem.MolToSmiles(mol))


def test_to_inchi_key_format():
    """InChIKey must be the standard 27-char dashed form."""
    from Grouper.exports import to_inchi_key
    key = to_inchi_key("CCO")
    parts = key.split("-")
    assert len(parts) == 3
    assert len(parts[0]) == 14
    assert len(parts[1]) == 10
    assert len(parts[2]) == 1
    assert all(p.isalnum() for p in parts)


def test_to_smarts_matches_self():
    """A SMARTS pattern derived from a molecule must match that
    molecule (substructure of itself)."""
    from Grouper.exports import to_smarts
    smarts = to_smarts("CCO")
    pattern = Chem.MolFromSmarts(smarts)
    mol = Chem.MolFromSmiles("CCO")
    assert mol.HasSubstructMatch(pattern)


def test_groupgraph_identifier_methods():
    """The bound shortcut methods produce the same output as the
    free functions."""
    from Grouper.exports import to_inchi, to_inchi_key, to_smarts
    gg = _ethanol_groupgraph()
    smi = gg.to_smiles()
    assert gg.to_inchi() == to_inchi(smi)
    assert gg.to_inchi_key() == to_inchi_key(smi)
    assert gg.to_smarts() == to_smarts(smi)


def test_invalid_smiles_to_identifier_raises():
    """Identifier methods refuse unparseable SMILES rather than
    silently returning RDKit's None."""
    from Grouper.exports import EmbedError, to_inchi, to_inchi_key, to_smarts
    for fn in (to_inchi, to_inchi_key, to_smarts):
        with pytest.raises(EmbedError, match="could not parse"):
            fn("not a smiles")
