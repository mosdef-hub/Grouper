
import typing as t
from copy import deepcopy

def run_performance_eval(
        nauty_path: str = "/Users/kieran/projects/molGrouper/packages/nauty2_8_8", 
        node_defs = None,
        n_runs: int = 3, 
        max_nodes: int = 6,
        verbose: bool = False
):
    """
    Run performance evaluation of the genGrouper generation featrure with growing sized graphs.
    """
    import os
    import time
    from genGrouper import exhaustive_generate
    import random
    import itertools

    if node_defs is None:
        raise ValueError("node_defs must be provided.")
    
    # Time the performance of the generation algorithm
    performance = {combo : { n:[] for n in range(1, max_nodes + 1)} for combo in itertools.combinations(node_defs, len(node_defs) - 1)}
    for def_combo in performance:
        for n_nodes in range(1, max_nodes + 1):
            total_time = 0
            for i in range(n_runs):
                random.seed(0)
                start_time = time.time()
                space = exhaustive_generate(
                    n_nodes, 
                    set(def_combo), 
                    nauty_path=nauty_path,
                    input_file_path="",
                    num_procs=16,
                    verbose=verbose)
                end_time = time.time()
                total_time += end_time - start_time
                performance[def_combo][n_nodes].append(total_time)
            print(f"Number of generated graphs: {len(space)}")
            print("")
    # plot_performance(performance, max_nodes)
    return performance
    
def plot_performance(performance, max_nodes):
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools

    # Prepare data for the plot
    fig, ax = plt.subplots()
    
    # Use a single colormap
    cmap = plt.get_cmap('inferno')
    norm = plt.Normalize()

    # Iterate through the performance data to plot
    for i, combo in enumerate(performance):
        for n_nodes in performance[combo]:
            times = performance[combo][n_nodes]
            mean_time = np.mean(times)
            color = cmap(norm(mean_time))

            ax.scatter(
                n_nodes, 
                i, 
                c=color
            )

    ax.set_xlabel("Number of nodes")
    ax.set_ylabel("Node type combination")
    
    # Adjust y-axis labels to display the node combinations
    ax.set_yticks(range(len(performance)))
    ax.set_yticklabels([str(combo) for combo in performance])

    # Add the color bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Only needed for the color bar
    plt.colorbar(sm, ax=ax, label='Mean Time (s)')

    plt.savefig("performance.png")



def convert_edges_to_nodetype(G):
    new_G = deepcopy(G)
    for edge in new_G.edges():
        for i in range(len(new_G.edges[edge]['ports'])):
            node_port_str_tuple = new_G.edges[edge]['ports'][i]
            src, dst = node_port_str_tuple
            srcnode, srcport = src.split('.')
            srcnode_type = new_G.nodes[srcnode]['type']

            dstnode, dstport = dst.split('.')
            dstnode_type = new_G.nodes[dstnode]['type']

            new_G.edges[edge]['ports'][i] = (f'{srcnode_type}.{srcport}', f'{dstnode_type}.{dstport}')

    return new_G


def multi_to_pair(multi: int, max_multi: int) -> t.Tuple[int, int]:
    """
    Convert a linear index `multi` to 2D coordinates (x, y) within a square grid.

    Parameters:
    - multi (int): Linear index to be converted.
    - max_multi (int): Maximum allowed value for the linear index.

    Returns:
    - Tuple[int, int]: Two-dimensional coordinates (x, y) within the grid.
    """
    if 1 <= multi <= max_multi:
        # Calculate x and y values based on the input multi
        x = (multi - 1) % int(max_multi**.5) + 1
        y = (multi - 1) // int(max_multi**.5) + 1

        x, y = x-1, y-1 # convert to 0-indexed
        return x, y
    else:
        raise ValueError("Input multi must be in the range 1 to max_multi.")
    
