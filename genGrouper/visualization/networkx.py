"""A module to visualize the graph structure of GroupGraph for seeing connectivity."""

import numpy as np
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import typing as t
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import cairosvg
from PIL import Image
from io import BytesIO
from Grouper import GroupGraph
from Grouper.utils import nxGroupGraph
import math

def generate_svg(num_ports, radius_node=50, radius_ports=10, distance_from_center=48, color='lightblue'):
    svg_header = '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200" viewBox="0 0 200 200">'
    svg_footer = '</svg>'
    
    # Center of the main circle
    center_x = 100
    center_y = 100

    # Main node circle
    main_circle = f'<circle cx="{center_x}" cy="{center_y}" r="{radius_node}" fill="{color}" stroke="black" stroke-width="2"/>'

    # Calculate port positions and draw them
    ports = []
    angle_step = 2 * math.pi / num_ports
    for i in range(num_ports):
        angle = i * angle_step
        port_x = center_x + distance_from_center * math.cos(angle)
        port_y = center_y + distance_from_center * math.sin(angle)
        port_circle = f'<circle cx="{port_x}" cy="{port_y}" r="{radius_ports}" fill="orange" stroke="black" stroke-width="2"/>'
        ports.append(port_circle)

    # Combine all SVG elements
    svg_content = "\n".join([main_circle] + ports)
    return f"{svg_header}\n{svg_content}\n{svg_footer}"

def generate_svg_for_node(num_ports, radius_node=50, radius_ports=10, distance_from_center=80, color = 'lightblue'):
    svg_code = generate_svg(num_ports, radius_node, radius_ports, distance_from_center, color)
    png_image = cairosvg.svg2png(bytestring=svg_code.encode('utf-8'), background_color=None)
    image = Image.open(BytesIO(png_image))
    return image

def nx_visualize(group_graph, layout=nx.spring_layout, **kwargs):
    if not kwargs:
        kwargs = dict(seed=1, k=0.5)

    if not isinstance(group_graph, nxGroupGraph):
        raise ValueError("Input must be an nxGroupGraph object.")
    if len(group_graph.nodes) == 0:
        raise ValueError("Input nxGroupGraph must have nodes to visualize.")
    components = nx.connected_components(group_graph)
    offset = 0
    pos = {}
    for component in components:
        component_subgraph = group_graph.subgraph(component)
        component_pos = layout(component_subgraph, center=(offset, 0), **kwargs)
        pos.update(component_pos)
        offset += 0.5 * np.sqrt(component_subgraph.number_of_nodes())

    labels = dict(group_graph.nodes.data("type"))
    largest_label = max([len(label) for label in labels.values()])
    node_size = 700 * largest_label + 200

    n_nodes = group_graph.number_of_nodes()
    dimensions = int(np.sqrt(n_nodes)) * 2 + 1 + largest_label // 3
    fig, ax = plt.subplots(figsize=(dimensions, dimensions))

    node_types = set(dict(group_graph.nodes.data("type")).values())
    cmap = sns.color_palette("colorblind", as_cmap=True)
    cmap = colors.ListedColormap(cmap)
    color_list = [colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    type_colors = {x: y for x, y in zip(node_types, color_list)}
    color_sequence = [type_colors[group_graph.nodes[node]["type"]] for node in group_graph.nodes]

    for c, (node, (x, y)) in zip(color_sequence, pos.items()):
        num_ports = len(group_graph.nodes[node].get("ports")) if 'ports' in group_graph.nodes[node] else 4    # Default to 4 ports if not specified
        image = generate_svg_for_node(num_ports, radius_node=80, color=c)
        
        imagebox = OffsetImage(image, zoom=0.4, resample=True)
        ab = AnnotationBbox(imagebox, (x, y), frameon=False, bboxprops=dict(edgecolor='black', boxstyle="circle,pad=0.3"))
        ax.add_artist(ab)
    
    nx.draw_networkx_edges(group_graph, pos, ax=ax, width=4, edge_color="black")

    x_values, y_values = zip(*pos.values())
    x_margin = (max(x_values) - min(x_values)) * 0.1 + 0.1
    y_margin = (max(y_values) - min(y_values)) * 0.1 + 0.1
    plt.xlim(min(x_values) - x_margin, max(x_values) + x_margin)
    plt.ylim(min(y_values) - y_margin, max(y_values) + y_margin)

    options = {
        "node_color": color_sequence,
        # "edge_color": edge_colors,
        "width": 4,
        "edge_cmap": plt.cm.hsv,
        "with_labels": True,
        "labels":labels,
        "font_weight":"bold",
        "edgecolors":"black"
    }
    nx.draw(group_graph, pos, **options)
    fig.canvas.draw()

    return fig
