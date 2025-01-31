"""A module to visualize the graph structure of GroupGraph for seeing connectivity."""

import math
from io import BytesIO

import cairosvg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import seaborn as sns
from matplotlib import colors
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from PIL import Image


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
        node1, port1, node2, port2, order = edge

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


def nx_visualize(group_graph, layout=nx.spring_layout, **kwargs):
    """Draw on a networkx canvas the current GroupGraph structure.

    Parameters
    ---------
    group_graph : genGrouper.GroupGraph
        The group graph to plot. Will size the canvas based on the size of the group graph.
        Labels will be placed on each node based on the "type" of the node.
    layout : function, optional, default nx.spring_layout
        Function to apply to the group_graph in order to determine position of nodes on canvas.
    kwargs : dict, optional, default None
        Extra parameters passed to layout function. Default kwargs are passed to nx.spring_layout.

    Returns
    -------
    fig : matplotlib.pyplot.Figure
        Drawn canvas of group_graph

    """
    if not kwargs:
        kwargs = dict(seed=1, k=0.5)

    components = nx.connected_components(group_graph)
    offset = 0
    pos = {}
    for component in components:
        component_subgraph = group_graph.subgraph(component)
        component_pos = layout(
            component_subgraph, center=(offset, 0), **kwargs
        )  # Seed layout for reproducibility
        pos.update(component_pos)
        offset += 0.5 * np.sqrt(
            component_subgraph.number_of_nodes()
        )  # add horizontal space between graphs

    labels = dict(group_graph.nodes.data("type"))
    largest_label = max([len(label) for label in labels.values()])
    node_size = 700 * largest_label + 200

    n_nodes = group_graph.number_of_nodes()
    dimensions = int(np.sqrt(n_nodes)) * 2 + 1 + largest_label // 3
    plt.figure(figsize=(dimensions, dimensions))

    # edge_colors = range(len(group_graph.edges())) # TODO: Add colors based on node edges
    node_types = set(dict(group_graph.nodes.data("type")).values())
    cmap = sns.color_palette("colorblind", as_cmap=True)
    cmap = colors.ListedColormap(cmap)
    color_list = [colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    type_colors = {x: y for x, y in zip(node_types, color_list)}
    color_sequence = [
        type_colors[group_graph.nodes[node]["type"]] for node in group_graph.nodes
    ]
    # colors = np.ones(len(group_graph.edges())) * 100

    options = {
        "node_color": color_sequence,
        # "edge_color": edge_colors,
        "width": 4,
        "edge_cmap": plt.cm.hsv,
        "with_labels": True,
        "node_size": node_size,
        "labels": labels,
        "font_weight": "bold",
        "edgecolors": "black",
    }
    fig = nx.draw(group_graph, pos, **options)
    x_values, y_values = zip(*pos.values())
    x_margin = (max(x_values) - min(x_values)) * 0.1 + 0.1
    y_margin = (max(y_values) - min(y_values)) * 0.1 + 0.1
    plt.xlim(min(x_values) - x_margin, max(x_values) + x_margin)
    plt.ylim(min(y_values) - y_margin, max(y_values) + y_margin)
    return fig


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
            (node1, port1, node2, port2, order) = edge
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
