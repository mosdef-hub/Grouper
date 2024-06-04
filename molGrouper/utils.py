
import typing as t
from copy import deepcopy

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
    

def run_performance_eval(
        nauty_path: str = "/Users/kieran/projects/molGrouper/packages/nauty2_8_8", 
        n_runs: int = 3, 
        max_nodes: int = 6,
        max_ports: int = 4,
        n_colors: int = 3,
        verbose: bool = False
):
    """
    Run performance evaluation of the molGrouper generation featrure with growing sized graphs.
    """
    import os
    import time
    from molGrouper.edge_coloring import process_nauty_vcolg_output
    import random

    max_nodes = list(range(2, max_nodes))
    generate_random_ports = lambda max_ports: [str(i) for i in range(random.randint(1, max_ports))]
    node_types = {chr(i+65): generate_random_ports(max_ports) for i in range(n_colors)}

    for n in max_nodes:
        total_time = 0
        for i in range(n_runs):
            os.system(f"{nauty_path}/geng {n} -ctf > /Users/kieran/projects/molGrouper/geng_out.txt")
            os.system(f"{nauty_path}/vcolg /Users/kieran/projects/molGrouper/geng_out.txt -T -m{n_colors} > vcolg_out.txt")
            # time the edge_coloring.py
            start = time.time()
            out = process_nauty_vcolg_output('/Users/kieran/projects/molGrouper/vcolg_out.txt', node_types, verbose=False)
            end = time.time()
            total_time += end - start
        print(f"Average time for {n} nodes: {total_time/n_runs}")

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
