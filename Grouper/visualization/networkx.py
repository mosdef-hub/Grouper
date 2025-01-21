import math
from io import BytesIO
from typing import Optional, Tuple

import cairosvg
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import colors
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw


def spring_layout(group_graph, iterations=50, k=1.0):
    """
    Computes a spring layout for a generic graph.

    Parameters:
        group_graph: The graph structure with `nodes` and edges
        iterations: Number of iterations for layout optimization
        k: Spring constant (controls node spacing)

    Returns:
        pos: A dictionary mapping nodes to 2D positions
    """
    nodes = list(group_graph.nodes)

    # Initialize positions randomly in 2D space
    pos = {node: np.random.rand(2) for node in nodes}

    # Compute repulsion and attraction forces iteratively
    for _ in range(iterations):
        forces = {node: np.zeros(2) for node in nodes}

        # Repulsive forces (Coulomb's law)
        for i, node1 in enumerate(nodes):
            for j, node2 in enumerate(nodes):
                if i >= j:
                    continue
                delta = pos[node1] - pos[node2]
                distance = np.linalg.norm(delta) + 1e-4  # Avoid division by zero
                force = (k**2 / distance**2) * delta / distance
                forces[node1] += force
                forces[node2] -= force

        # Attractive forces (Hooke's law)
        for edge in group_graph.edges:
            (node1, port1, node2, port2) = edge
            delta = pos[node1] - pos[node2]
            distance = np.linalg.norm(delta) + 1e-4
            force = -(distance**2 / k) * delta / distance
            forces[node1] += force
            forces[node2] -= force

        # Update positions
        for node in nodes:
            pos[node] += 0.01 * forces[node]  # Small step size for stability

    return pos


# Calculate positions for ports around a node, adjusted by zoom factor and distance
def calculate_port_positions(
    center_x, center_y, num_ports, distance_from_center, zoom_factor=0.4
):
    # Adjust distance from center by the zoom factor
    scaled_distance = distance_from_center
    angle_step = 2 * math.pi / num_ports
    port_positions = []
    for i in range(0, -num_ports - 1, -1):
        angle = i * angle_step
        port_x = center_x + scaled_distance * math.cos(angle)
        port_y = center_y + scaled_distance * math.sin(angle)
        port_positions.append((port_x, port_y))
    return port_positions


# Generate node images with ports and labels
def generate_svg(
    num_ports,
    node_type,
    radius_node=50,
    radius_ports=30,
    distance_from_center=48,
    color="lightblue",
):
    svg_header = '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200" viewBox="0 0 200 200">'
    svg_footer = "</svg>"

    center_x = 100
    center_y = 100

    # Main node circle
    main_circle = f'<circle cx="{center_x}" cy="{center_y}" r="{radius_node}" fill="{color}" stroke="black" stroke-width="2"/>'
    main_label = f'<text x="{center_x}" y="{center_y}" font-size="25" text-anchor="middle" alignment-baseline="middle" fill="black">{node_type}</text>'

    # Port positions around the main node
    ports = []
    angle_step = 2 * math.pi / num_ports
    for i in range(num_ports):
        angle = i * angle_step
        port_x = center_x + distance_from_center * math.cos(angle)
        port_y = center_y + distance_from_center * math.sin(angle)

        # Port circle
        port_circle = f'<circle cx="{port_x}" cy="{port_y}" r="{radius_ports}" fill="orange" stroke="black" stroke-width="2"/>'
        ports.append(port_circle)

        # Port number label (centered in the circle)
        port_label = f'<text x="{port_x}" y="{port_y}" font-size="15" text-anchor="middle" alignment-baseline="middle" fill="black">{i}</text>'
        ports.append(port_label)

    svg_content = "\n".join([main_circle, main_label] + ports)
    return f"{svg_header}\n{svg_content}\n{svg_footer}"


# Convert SVG to image for display
def generate_svg_for_node(
    num_ports,
    node_type,
    radius_node=50,
    radius_ports=10,
    distance_from_center=80,
    color="lightblue",
):
    svg_code = generate_svg(
        num_ports, node_type, radius_node, radius_ports, distance_from_center, color
    )
    png_image = cairosvg.svg2png(
        bytestring=svg_code.encode("utf-8"), background_color=None
    )
    image = Image.open(BytesIO(png_image))
    return image


# Manual visualization without using NetworkX for drawing
def visualize(group_graph, pos=None):
    """
    Visualize a graph with optional custom positioning for nodes.

    Parameters:
    - group_graph: A NetworkX graph object to be visualized.
    - pos: A dictionary specifying positions for nodes (optional).
           If None, a default layout will be used.
    """
    if pos is None:
        pos = spring_layout(group_graph, iterations=50, k=1.0)
    fig, ax = plt.subplots(
        figsize=(5, 5),
    )

    # Set axis limits
    x_values, y_values = zip(*pos.values())
    x_margin = (max(x_values) - min(x_values)) * 0.1 + 0.1
    y_margin = (max(y_values) - min(y_values)) * 0.1 + 0.1
    plt.xlim(min(x_values) - x_margin, max(x_values) + x_margin)
    plt.ylim(min(y_values) - y_margin, max(y_values) + y_margin)

    # Color map for nodes
    node_types = set([n.type for n in group_graph.nodes.values()])
    cmap = sns.color_palette("colorblind", as_cmap=True)
    cmap = colors.ListedColormap(cmap)
    color_list = [colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    type_colors = {x: y for x, y in zip(node_types, color_list)}

    # Dictionary to store port positions
    port_positions = {}

    # Dictionary to keep track of used ports for each node
    used_ports = {node: set() for node in group_graph.nodes}

    # Set zoom factor for the image scaling
    zoom_factor = 0.4

    # Draw nodes with images and calculate scaled port positions
    for node, (x, y) in pos.items():
        num_ports = len(group_graph.nodes[node].ports)
        color = type_colors[group_graph.nodes[node].type]
        image = generate_svg_for_node(
            num_ports, group_graph.nodes[node].type, radius_node=80, color=color
        )
        imagebox = OffsetImage(image, zoom=zoom_factor, resample=True)
        ab = AnnotationBbox(imagebox, (x, y), frameon=False, zorder=1)
        ax.add_artist(ab)
        bbox = ab.get_window_extent(renderer=fig.canvas.get_renderer())
        bbox_data = bbox.transformed(ax.transData.inverted())
        radius = max(bbox_data.width, bbox_data.height) / 2
        # Adjust distance from center by the radius of the node image and make slightly smaller since the image is larger
        distance_from_center = radius * 0.85

        port_positions[node] = calculate_port_positions(
            x,
            y,
            num_ports,
            distance_from_center,
        )

    # Manually draw edges between specified ports
    for i, edge in enumerate(group_graph.edges):
        node1, port1, node2, port2 = edge

        # Ensure ports are not reused
        if port1 in used_ports[node1] or port2 in used_ports[node2]:
            print(
                f"Warning: Port {port1} of node {node1} or port {port2} of node {node2} is already used!"
            )
            continue  # Skip this edge if the port is already used

        # Get port positions for both nodes
        port_pos1 = port_positions[node1][port1]
        port_pos2 = port_positions[node2][port2]

        # Draw edge between the specified ports
        ax.plot(
            [port_pos1[0], port_pos2[0]],
            [port_pos1[1], port_pos2[1]],
            color="black",
            linestyle="dashdot",
            linewidth=2,
            zorder=2,
        )
        # ax.text(
        #     0.5 * (port_pos1[0] + port_pos2[0]),
        #     0.5 * (port_pos1[1] + port_pos2[1]),
        #     f'Edge {i}',
        #     fontsize=12
        # )

        # Mark ports as used
        used_ports[node1].add(port1)
        used_ports[node2].add(port2)

    # Hide axes
    ax.axis("off")

    return fig


def visualize_node_trace1(
    nt,
    highlight_color: Tuple[float, float, float] = (0.8, 0.0, 0.8),
    sanitize_smiles: bool = False,
):
    # TODO: flag for hydrogens in SMARTS string
    mol = Chem.MolFromSmarts(nt.smarts)  # -CH2-
    if not mol:
        raise ValueError(f"Could not parse SMARTS: {nt.smarts}")
    mol.UpdatePropertyCache()
    Chem.rdmolops.AddHs(
        mol
    )  # Adding hydrogens sometimes seems to mess up the SMARTS when bond order is specified.
    d = Draw.rdMolDraw2D.MolDraw2DCairo(200, 200)
    patt = Chem.MolFromSmiles(nt.node.smiles, sanitize=sanitize_smiles)
    if not patt:
        raise ValueError(f"Could not parse SMILES: {nt.node.smiles}")
    hit_ats = list(mol.GetSubstructMatch(patt))
    if not hit_ats:
        raise ValueError(
            f"Could not find substructure match for node {nt.node.type} in SMARTS: {nt.smarts}"
        )
    atom_cols = {}
    for i, at in enumerate(hit_ats):
        atom_cols[at] = highlight_color
    # Necessary for Bonds if we want to highlight these at some point
    # hit_bonds = []
    # for bond in patt.GetBonds():
    #     aid1 = hit_ats[bond.GetBeginAtomIdx()]
    #     aid2 = hit_ats[bond.GetEndAtomIdx()]
    #     hit_bonds.append(mol.GetBondBetweenAtoms(aid1,aid2).GetIdx())
    # bond_cols = {}
    # for i, bd in enumerate(hit_bonds):
    #     bond_cols[bd] = (0.0,0.0,0.8)
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(
        d,
        mol,
        highlightAtoms=hit_ats,
        highlightAtomColors=atom_cols,
        # highlightBonds=hit_bonds, # this is also only needed if we want to draw bonds
        # highlightBondColors=bond_cols
    )
    return d, mol, hit_ats


def visualize_node_trace2(
    node_trace,
    text: Optional[str | bool] = "Matching Subgraph",
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
        node_trace.node.smiles, sanitize=sanitize_smiles
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
        font = ImageFont.truetype("Chalkduster.ttf", 16)
        _, descent = font.getmetrics()
        text_width = font.getmask(text).getbbox()[2]
        text_height = font.getmask(text).getbbox()[3] + descent
        pos = (0.5 * height - 0.5 * text_width, height - text_height)
        # Set font and text
        text_color = tuple(int(x * 255) for x in highlight_color)
        draw.text(pos, text, font=font, fill=text_color)
    return img
