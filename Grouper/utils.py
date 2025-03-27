from copy import deepcopy
from typing import Any, Callable, Dict, List, Sequence, Tuple

import networkx as nx
import torch_geometric

from Grouper import GroupGraph


def run_performance_eval(
    nauty_path: str = "/Users/kieran/projects/molGrouper/packages/nauty2_8_8",
    node_defs=None,
    n_runs: int = 3,
    max_nodes: int = 6,
    n_cpus: int = 30,
    verbose: bool = False,
):
    """
    Run performance evaluation of the Grouper generation featrure with growing sized graphs.
    """
    import itertools
    import random
    import time

    from Grouper import exhaustive_generate

    if node_defs is None:
        raise ValueError("node_defs must be provided.")

    # Time the performance of the generation algorithm
    performance = {
        combo: {n: [] for n in range(1, max_nodes + 1)}
        for combo in itertools.combinations(node_defs, len(node_defs) - 1)
    }
    for def_combo in performance:
        for n_nodes in range(1, max_nodes + 1):
            for i in range(n_runs):
                random.seed(0)
                start_time = time.time()
                space = exhaustive_generate(
                    n_nodes,
                    set(def_combo),
                    nauty_path=nauty_path,
                    input_file_path="",
                    num_procs=n_cpus,
                    verbose=verbose,
                )
                end_time = time.time()
                performance[def_combo][n_nodes].append(end_time - start_time)
            print(f"Number of generated graphs: {len(space)}")
            print("")
    # plot_performance(performance, max_nodes)
    return performance


def plot_performance(performance, max_nodes):
    import matplotlib.pyplot as plt
    import numpy as np

    # Prepare data for the plot
    fig, ax = plt.subplots()

    # Use a single colormap
    cmap = plt.get_cmap("inferno")
    norm = plt.Normalize()

    # Iterate through the performance data to plot
    for i, combo in enumerate(performance):
        for n_nodes in performance[combo]:
            times = performance[combo][n_nodes]
            mean_time = np.mean(times)
            color = cmap(mean_time)

            ax.scatter(n_nodes, i, c=color)

    ax.set_xlabel("Number of nodes")
    ax.set_ylabel("Node type combination")

    # Adjust y-axis labels to display the node combinations
    ax.set_yticks(range(len(performance)))
    ax.set_yticklabels([str(combo) for combo in performance])

    # Add the color bar
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array([])  # Only needed for the color bar
    plt.colorbar(sm, ax=ax, label="Mean Time (s)")

    plt.savefig("performance.png")


def convert_edges_to_nodetype(G):
    new_G = deepcopy(G)
    for edge in new_G.edges():
        for i in range(len(new_G.edges[edge]["ports"])):
            node_port_str_tuple = new_G.edges[edge]["ports"][i]
            src, dst = node_port_str_tuple
            srcnode, srcport = src.split(".")
            srcnode_type = new_G.nodes[srcnode]["type"]

            dstnode, dstport = dst.split(".")
            dstnode_type = new_G.nodes[dstnode]["type"]

            new_G.edges[edge]["ports"][i] = (
                f"{srcnode_type}.{srcport}",
                f"{dstnode_type}.{dstport}",
            )

    return new_G


class nxGroupGraph(nx.Graph):
    """A graph with ports as parts of nodes that can be connected to other ports."""

    def __init__(self, node_types: Dict[str, List] = None):
        """
        Initialize a GroupGraph.

        Parameters:
        - node_types (dict): Dictionary of node types and their corresponding ports.

        Raises:
        - TypeError: If node_types is not a dictionary.
        - ValueError: If keys in node_types are not of the same type or values are not lists.
        """
        super(nxGroupGraph, self).__init__()
        if node_types is None:
            node_types = {}

        if not isinstance(node_types, dict):
            raise TypeError(
                "node_types must be a dictionary of node types and their ports"
            )
        if len(node_types) != 0:
            holder = next(iter(node_types.keys()))  # raises error if empty
            holder_type = type(holder)
            if not all(isinstance(k, holder_type) for k in node_types.keys()):
                raise ValueError("All keys in node_types must be of the same type")
            if not all(isinstance(v, list) for v in node_types.values()):
                raise ValueError("All values in node_types must be lists")

        self.node_types = node_types

    def __str__(self) -> str:
        """
        Return a string representation of the GroupGraph.

        Returns:
        - str: String representation of nodes and edges.
        """
        return f"Nodes: {['{} ({}) {}'.format(d[0], d[1]['type'], d[1]['ports']) for d in self.nodes.data()]}\nEdges: {','.join(str(tuple(*e[-1])) for e in self.edges.data('ports'))}\n"

    def __repr__(self) -> str:
        """
        Return a representation of the GroupGraph.

        Returns:
        - str: Representation of the GroupGraph.
        """
        return f"GroupGraph({','.join(str(tuple(*e[-1])) for e in self.edges.data('ports'))})"

    def __eq__(self, other) -> bool:
        """
        Check if two GroupGraph objects are equal.

        Parameters:
        - other: Another GroupGraph object.

        Returns:
        - bool: True if equal, False otherwise.
        """
        # return self.nodes(data=True) == other.nodes(data=True) and self.edges(data=True) == other.edges(data=True)
        return nx.utils.misc.graphs_equal(self, other)

    def __ne__(self, other) -> bool:
        """
        Check if two GroupGraph objects are not equal.

        Parameters:
        - other: Another GroupGraph object.

        Returns:
        - bool: True if not equal, False otherwise.
        """
        return not self.__eq__(other)

    def __bool__(self) -> bool:
        """
        Check if the GroupGraph is non-empty.

        Returns:
        - bool: True if non-empty, False otherwise.
        """
        return len(self.nodes) > 0 and len(self.edges) > 0

    def add_node(
        self, nodeID: Any, node_type: str, smarts: str, hubs: Sequence[int]
    ) -> None:
        """
        Add a node to the GroupGraph.

        Parameters:
        - nodeID: Identifier for the node.
        - node_type: Type of the node.

        Raises:
        - Exception: If the node is already present in the Graph.
        """
        # Sanity check to see if the node is already present in Graph
        if nodeID in self.nodes:
            raise AttributeError(f"Node: {nodeID} is already present in Graph")
        super(nxGroupGraph, self).add_node(nodeID)
        self.nodes[nodeID]["type"] = node_type
        self.nodes[nodeID]["ports"] = self.node_types[self.nodes[nodeID]["type"]]
        self.nodes[nodeID]["smarts"] = smarts
        self.nodes[nodeID]["hubs"] = hubs

    def add_edge(
        self, node_port_1: Tuple[Any, Any], node_port_2: Tuple[Any, Any]
    ) -> None:
        """
        Add an edge between two nodes in the GroupGraph.

        Parameters:
        - node_port_1: Tuple of node and port on the first node.
        - node_port_2: Tuple of node and port on the second node.

        Raises:
        - Exception: If a node has no free ports or if nodes or ports are not present in the Graph.
        """
        if self.n_free_ports(node_port_1[0]) <= 0:
            raise AttributeError(f"Node: {node_port_1[0]} has no free ports!")
        if self.n_free_ports(node_port_2[0]) <= 0:
            raise AttributeError(f"Node: {node_port_2[0]} has no free ports!")
        # check if port is already occupied
        for edge in self.edges(data=True):
            for node_port in edge[-1]["ports"]:
                src_node, src_port = node_port[0].split(".")
                dst_node, dst_port = node_port[1].split(".")
                if src_node == node_port_1[0]:
                    if src_port == node_port_1[1]:
                        raise AttributeError(
                            f"Node: {node_port_1[0]}.{node_port_1[1]} is already occupied!"
                        )
                if src_node == node_port_2[0]:
                    if src_port == node_port_2[1]:
                        raise AttributeError(
                            f"Node: {node_port_2[0]}.{node_port_2[1]} is already occupied!"
                        )
                if dst_node == node_port_1[0]:
                    if dst_port == node_port_1[1]:
                        raise AttributeError(
                            f"Node: {node_port_1[0]}.{node_port_1[1]} is already occupied!"
                        )
                if dst_node == node_port_2[0]:
                    if dst_port == node_port_2[1]:
                        raise AttributeError(
                            f"Node: {node_port_2[0]}.{node_port_2[1]} is already occupied!"
                        )

        edge_ports = []

        for n, p in [
            (node_port_1[0], node_port_1[1]),
            (node_port_2[0], node_port_2[1]),
        ]:
            # Sanity check to see if the nodes and ports are present in Graph
            if n not in self.nodes:
                raise AttributeError(f"Node: {p} is not present in Graph")
            if p not in self.nodes(data=True)[n]["ports"]:
                raise AttributeError(f"Port: {p} is incorrect for Node: {n}!")

            edge_ports.append(str(n) + "." + str(p))

        # Add the port points as edge attributes
        if self.has_edge(node_port_1[0], node_port_2[0]):
            self.edges[node_port_1[0], node_port_2[0]]["ports"].append(edge_ports)
        else:
            super(nxGroupGraph, self).add_edge(
                node_port_1[0], node_port_2[0], ports=[edge_ports]
            )

        # remove ports from portsList

    def make_undirected(self):
        """
        Convert the GroupGraph to an undirected graph.

        Adds bonds and updates edge ports accordingly.
        """
        for edge in self.edges:
            node_port = tuple(self.edges[edge]["ports"][0])
            node_s, node_t = edge
            port_s, port_t = edge["ports"][0], edge["ports"][1]
            self.add_bond(node_t, port_t, node_s, port_s)
            self.edges[edge]["ports"].append(f"{node_port[1]}.{node_port[0]}")

    def n_free_ports(self, nodeID: Any) -> int:
        """
        Get the number of free ports on a node.

        Parameters:
        - nodeID: Identifier for the node.

        Returns:
        - int: Number of free ports.
        """
        # num ports - num edges - num edges with node as target
        occupied_ports = 0
        if len(self.edges()) == 0:
            return len(self.nodes[nodeID]["ports"])
        for e in self.edges(data=True):
            for node_port in e[-1]["ports"]:
                if (
                    node_port[0].split(".")[0] == nodeID
                    or node_port[1].split(".")[0] == nodeID
                ):
                    occupied_ports += 1
        return (
            len(self.nodes[nodeID]["ports"]) - occupied_ports
        )  # total number of ports - number of ports with edges

    def to_PyG_Data(
        self,
        node_descriptor_generator: Callable[[str], Sequence[float]],
        max_ports: int = 0,
    ) -> torch_geometric.data.Data:
        """
        Convert the GroupGraph to a PyG Data representation.

        Parameters:
        - node_descriptor_generater (Callable[[str], Sequence[float]]): Callable to generate node descriptors. Takes in a SMARTS string and returns a sequence of floats.

        Returns:
        - torch_geometric.data.Data: PyG Data representation.
        """
        import torch
        from torch_geometric.data import Data

        one_hot_vector = lambda index, num_classes: torch.eye(num_classes)[index]
        max_ports = (
            max(len(v) for v in self.node_types.values())
            if max_ports == 0
            else max_ports
        )

        # Create the node features
        dummy_feature = node_descriptor_generator(self.nodes(data=True)[0]["smarts"])
        node_features = torch.zeros((len(self.nodes), len(dummy_feature))).to(
            torch.float64
        )
        for i, n in enumerate(self.nodes(data=True)):
            node_features[i] = node_descriptor_generator(n[1]["smarts"])

        # Create the edge index
        edge_index = torch.zeros((2, len(self.edges)), dtype=torch.long)
        for i, e in enumerate(self.edges):
            edge_index[0, i] = list(self.nodes).index(e[0])
            edge_index[1, i] = list(self.nodes).index(e[1])

        # Create the edge features
        edge_features = torch.zeros(len(self.edges), max_ports * 2)

        for i, e in enumerate(self.edges(data=True)):
            # get nodes and ports
            node_ports = [node_port.split(".") for node_port in e[2]["ports"][0]]

            # Convert node_ports to int from string
            node_ports = [(int(node), int(port)) for node, port in node_ports]

            # convert the ports to one-hot vectors
            port_s = self.nodes(data=True)[node_ports[0][0]]["ports"]
            port_t = self.nodes(data=True)[node_ports[1][0]]["ports"]
            port_index_s = port_s.index(node_ports[0][1])
            port_index_t = port_t.index(node_ports[1][1])
            edge_features[i] = torch.cat(
                [
                    one_hot_vector(port_index_s, max_ports),
                    one_hot_vector(port_index_t, max_ports),
                ]
            ).to(torch.float64)

        return Data(x=node_features, edge_index=edge_index, edge_attr=edge_features)


def convert_to_nx(G: GroupGraph) -> nxGroupGraph:
    """
    Convert a GroupGraph to a networkx graph.

    Parameters:
    - G (GroupGraph): The GroupGraph to convert.

    Returns:
    - nxGroupGraph: The networkx representation of the GroupGraph.
    """
    if not isinstance(G, GroupGraph):
        raise TypeError("G must be a GroupGraph")
    if G.nodes == {}:
        return nxGroupGraph()
    nxG = nxGroupGraph(G.node_types)
    for node_id, node in G.nodes.items():
        nxG.add_node(node_id, node.type, node.pattern, node.hubs)
    for edge in G.edges:
        src = (edge[0], edge[1])
        dst = (edge[2], edge[3])
        nxG.add_edge(src, dst)
    return nxG
