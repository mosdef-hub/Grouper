"""A module to visualize the graph structure of GroupGraph for seeing connectivity."""

import numpy as np
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def visualize(group_graph, layout=nx.spring_layout, **kwargs):
	"""Draw on a networkx canvas the current GroupGraph structure.

	Parameters
	---------
	group_graph : molGrouper.GroupGraph
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
		kwargs = dict(
			seed=1, k=0.5
		)

    components = nx.connected_components(group_graph)
    offset = 0
    pos = {}
    for component in components:
        component_subgraph = group_graph.subgraph(component)
        component_pos = layout(component_subgraph, center=(offset, 0), **kwargs)  # Seed layout for reproducibility
        pos.update(component_pos)
        offset += 0.5 * np.sqrt(component_subgraph.number_of_nodes()) # add horizontal space between graphs

    labels = dict(group_graph.nodes.data("type"))
    largest_label = max([len(label) for label in labels.values()])
    node_size = 700*largest_label + 200

    n_nodes = group_graph.number_of_nodes()
    dimensions = int(np.sqrt(n_nodes)) * 2 + 1 + largest_label // 3
    plt.figure(figsize=(dimensions, dimensions))

    # edge_colors = range(len(group_graph.edges())) # TODO: Add colors based on node edges
    node_types = set(dict(group_graph.nodes.data("type")).values())
    cmap = sns.color_palette("colorblind", as_cmap=True)
    cmap = colors.ListedColormap(cmap)
    color_list = [colors.rgb2hex(cmap(i)) for i in range(cmap.N)]
    type_colors = {x:y for x, y in zip(node_types, color_list)}
    color_sequence = [type_colors[group_graph.nodes[node]["type"]] for node in group_graph.nodes]
    # colors = np.ones(len(group_graph.edges())) * 100

    options = {
        "node_color": color_sequence,
        # "edge_color": edge_colors,
        "width": 4,
        "edge_cmap": plt.cm.hsv,
        "with_labels": True,
        "node_size":node_size,
        "labels":labels,
        "font_weight":"bold",
        "edgecolors":"black"
    }
    fig = nx.draw(group_graph, pos, **options)
    x_values, y_values = zip(*pos.values())
    x_margin = (max(x_values) - min(x_values)) * 0.1 + 0.1
    y_margin = (max(y_values) - min(y_values)) * 0.1 + 0.1
    plt.xlim(min(x_values) - x_margin, max(x_values) + x_margin)
    plt.ylim(min(y_values) - y_margin, max(y_values) + y_margin)
    return fig
