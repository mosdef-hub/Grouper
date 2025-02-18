from io import BytesIO
from typing import Optional, Tuple, Union

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw

from Grouper.libraries.Libraries import NodeTrace
from Grouper import Group


def visualize_node_trace(
    node_trace: NodeTrace,
    text: Optional[Union[str, bool]] = "Matching Subgraph",
    sanitize_smiles: bool = False,
    highlight_color: Tuple[float, float, float, float] = (1.0, 0.5, 0.5, 1.0),
    draw_options: dict = None,
) -> Image.Image:
    """A method to visualize a node trace, which definies a node and what that node is connected to.

    Uses RDKit and PIL to create an image of the Group Node and where is definited
    to be located in a molecule. Extra arguments can help tune the look of this output image.

    Parameters
    ----------
    node_trace : Grouper.libraries.Libraries.NodeTrace
        The Node with chemical information (SMARTS) to parse and visualize
    text : str, optional, default="Matching Subgraph"
        Text to place as a legend under the graph. The color matches the highlighted portion,
        which indicates these atoms would actually be contained in the function group Node.
    sanitize_smiles : bool, optional, default=False
        Argument that can be passed to rdkit.Chem.MolFromSmiles. Sanitization can strip hydrogens.
        See https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.SmilesParserParams
    highlight_color : Tuple[float, float, float, float], optional, default=(1.0, 0.5, 0.5, 1.0)
        RGBA color to highlight the atoms in the subgraph. Default is a light red.
    draw_options : dict, optional, default=None
        If False, use default options, specified in this function. Optionally, require your own,
        including the `bond_width` and `highlight_radius` for the drawing.
        Other options explained here https://greglandrum.github.io/rdkit-blog/posts/2023-05-26-drawing-options-explained.html

    Notes
    -----
    This can fail if SMARTS or SMILES are definied incorrectly, or RDKit has trouble parsing them.
    Parsing info can be found here:

    1. Daylight SMARTS: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
    2. Daylight SMILES:  https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
    3. RDKit SMARTS: https://www.rdkit.org/docs/RDKit_Book.html#smarts-support-and-extensions
    4. RDKit SMILES: https://www.rdkit.org/docs/RDKit_Book.html#smiles-support-and-extensions
    """
    mol = Chem.MolFromSmarts(node_trace.smarts)  # -CH2-
    if not mol:
        raise ValueError(f"Could not parse SMARTS: {node_trace.smarts}")

    # Clean ups structures
    mol.UpdatePropertyCache()
    Chem.rdmolops.AddHs(
        mol
    )  # Adding hydrogens sometimes seems to mess up the SMARTS when bond order is specified.

    # initialize SMARTS and SMILES objects for molecule
    # note, sanitize can sometimes improve the implicit/explicit hydrogens found in RDKit. Only turn on if you verify the SMARTS string is matched properly
    smarts_subgraph = Chem.MolFromSmiles(
        node_trace.node.smarts, sanitize=sanitize_smiles
    )
    if not smarts_subgraph:
        raise ValueError(f"Could not parse SMILES: {node_trace.node.smiles}")
    match_atoms = list(mol.GetSubstructMatch(smarts_subgraph))
    if not match_atoms:
        raise ValueError(
            f"Could not find substructure match for node {node_trace.node.type} in SMARTS: {node_trace.smarts}"
        )

    # identify atoms to highlight
    atom_cols = {}
    for at in match_atoms:
        atom_cols[at] = highlight_color

    # Initialize drawing canvas and options
    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(200, 200)
    if not draw_options:  # can pass draw_options
        draw_options = drawer.drawOptions()
        draw_options.dummiesAreAttachments = True
        draw_options.highlightRadius = 0.4
        draw_options.highlightBondWidthMultiplier = 6
    else:
        drawer.SetDrawOptions(draw_options)
    draw_options.setHighlightColour(highlight_color)

    # Drawing steps
    drawer.DrawMolecule(mol, highlightAtoms=match_atoms)
    drawer.FinishDrawing()
    bio = BytesIO(drawer.GetDrawingText())
    img = Image.open(bio)
    draw = ImageDraw.Draw(img)
    height, _ = img.height, img.width

    # Determine where to place label below object
    if text:
        font = ImageFont.load_default()
        _, descent = font.getmetrics()
        text_width = font.getmask(text).getbbox()[2]
        text_height = font.getmask(text).getbbox()[3] + descent
        pos = (0.5 * height - 0.5 * text_width, height - text_height)
        # Set font and text
        text_color = tuple(int(x * 255) for x in highlight_color)
        draw.text(pos, text, font=font, fill=text_color)
    return img


def visualize_group(
    node: Group,
    text: Optional[Union[str, bool]] = "",
    draw_numbers = False,
    sanitize_smiles: bool = False,
    highlight_color: Tuple[float, float, float, float] = (1.0, 0.5, 0.5, 1.0),
    draw_options: dict = None,
) -> Image.Image:
    """A method to visualize a node trace, which definies a node and what that node is connected to.

    Uses RDKit and PIL to create an image of the Group Node and where is definited
    to be located in a molecule. Extra arguments can help tune the look of this output image.

    Parameters
    ----------
    node_trace : Grouper.libraries.Libraries.NodeTrace
        The Node with chemical information (SMARTS) to parse and visualize
    text : str, optional, default="Matching Subgraph"
        Text to place as a legend under the graph. The color matches the highlighted portion,
        which indicates these atoms would actually be contained in the function group Node.
    draw_numbers : bool, optional, default=False
        If True, draw the atom indices on the image.
    sanitize_smiles : bool, optional, default=False
        Argument that can be passed to rdkit.Chem.MolFromSmiles. Sanitization can strip hydrogens.
        See https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.SmilesParserParams
    highlight_color : Tuple[float, float, float, float], optional, default=(1.0, 0.5, 0.5, 1.0)
        RGBA color to highlight the atoms in the subgraph. Default is a light red.
    draw_options : dict, optional, default=None
        If False, use default options, specified in this function. Optionally, require your own,
        including the `bond_width` and `highlight_radius` for the drawing.
        Other options explained here https://greglandrum.github.io/rdkit-blog/posts/2023-05-26-drawing-options-explained.html

    Notes
    -----
    This can fail if SMARTS or SMILES are definied incorrectly, or RDKit has trouble parsing them.
    Parsing info can be found here:

    1. Daylight SMARTS: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
    2. Daylight SMILES:  https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
    3. RDKit SMARTS: https://www.rdkit.org/docs/RDKit_Book.html#smarts-support-and-extensions
    4. RDKit SMILES: https://www.rdkit.org/docs/RDKit_Book.html#smiles-support-and-extensions
    """
    mol = Chem.MolFromSmarts(node.smarts)  # -CH2-
    if not mol:
        raise ValueError(f"Could not parse SMARTS: {node.smarts}")

    # Clean ups structures
    mol.UpdatePropertyCache()
    Chem.rdmolops.AddHs(
        mol
    )  # Adding hydrogens sometimes seems to mess up the SMARTS when bond order is specified.

    # initialize SMARTS and SMILES objects for molecule
    # note, sanitize can sometimes improve the implicit/explicit hydrogens found in RDKit. Only turn on if you verify the SMARTS string is matched properly
    smarts_subgraph = Chem.MolFromSmiles(
        node.smarts, sanitize=sanitize_smiles
    )
    if not smarts_subgraph:
        raise ValueError(f"Could not parse SMILES: {node.node.smiles}")
    match_atoms = list(mol.GetSubstructMatch(smarts_subgraph))
    if not match_atoms:
        raise ValueError(
            f"Could not find substructure match for node {node.node.type} in SMARTS: {node.smarts}"
        )

    # Initialize drawing canvas and options
    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(200, 200)
    if not draw_options:  # can pass draw_options
        draw_options = drawer.drawOptions()
        draw_options.dummiesAreAttachments = True
    else:
        drawer.SetDrawOptions(draw_options)
    draw_options.setHighlightColour(highlight_color)

    if draw_numbers:
        # Add indices to atoms
        for i, atom in enumerate(mol.GetAtoms()):
            atom.SetProp("atomNote", str(i))

    # Drawing steps
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    bio = BytesIO(drawer.GetDrawingText())
    img = Image.open(bio)
    draw = ImageDraw.Draw(img)
    height, _ = img.height, img.width

    # Determine where to place label below object
    if text:
        font = ImageFont.truetype("Chalkduster.ttf", 16)
        _, descent = font.getmetrics()
        text_width = font.getmask(text).getbbox()[2]
        text_height = font.getmask(text).getbbox()[3] + descent
        pos = (0.5 * height - 0.5 * text_width, height - text_height)
        # Set font and text
        text_color = tuple(int(x * 255) for x in highlight_color)
        draw.text(pos, text, font=font, fill=text_color)
    return img
