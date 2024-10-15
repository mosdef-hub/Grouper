import math
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import seaborn as sns
from matplotlib import colors
from PIL import Image
from io import BytesIO
import cairosvg

# Calculate positions for ports around a node, adjusted by zoom factor and distance
def calculate_port_positions(center_x, center_y, num_ports, distance_from_center, zoom_factor=0.4):
    # Adjust distance from center by the zoom factor
    scaled_distance = distance_from_center 
    angle_step = 2 * math.pi / num_ports
    port_positions = []
    for i in range(0, -num_ports-1, -1):
        angle = i * angle_step
        port_x = center_x + scaled_distance * math.cos(angle)
        port_y = center_y + scaled_distance * math.sin(angle)
        port_positions.append((port_x, port_y))
    return port_positions


# Generate node images with ports and labels
def generate_svg(num_ports, node_type, radius_node=50, radius_ports=30, distance_from_center=48, color='lightblue'):
    svg_header = '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200" viewBox="0 0 200 200">'
    svg_footer = '</svg>'
    
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
def generate_svg_for_node(num_ports, node_type, radius_node=50, radius_ports=10, distance_from_center=80, color='lightblue'):
    svg_code = generate_svg(num_ports, node_type, radius_node, radius_ports, distance_from_center, color)
    png_image = cairosvg.svg2png(bytestring=svg_code.encode('utf-8'), background_color=None)
    image = Image.open(BytesIO(png_image))
    return image


# Manual visualization without using NetworkX for drawing
def visualize(group_graph, pos):
    fig, ax = plt.subplots(figsize=(5, 5), )

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
        image = generate_svg_for_node(num_ports, group_graph.nodes[node].type, radius_node=80, color=color)
        imagebox = OffsetImage(image, zoom=zoom_factor, resample=True)
        ab = AnnotationBbox(imagebox, (x, y), frameon=False, zorder=1)
        ax.add_artist(ab)
        bbox = ab.get_window_extent(renderer=fig.canvas.get_renderer())
        bbox_data = bbox.transformed(ax.transData.inverted())
        radius = max(bbox_data.width, bbox_data.height) / 2
        # Adjust distance from center by the radius of the node image and make slightly smaller since the image is larger
        distance_from_center = radius * .85

        port_positions[node] = calculate_port_positions(x, y, num_ports, distance_from_center, )

    # Manually draw edges between specified ports
    for i, edge in enumerate(group_graph.edges):
        node1, port1, node2, port2 = edge

        # Ensure ports are not reused
        if port1 in used_ports[node1] or port2 in used_ports[node2]:
            print(f"Warning: Port {port1} of node {node1} or port {port2} of node {node2} is already used!")
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
            zorder = 2
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