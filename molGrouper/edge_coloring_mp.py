import multiprocessing as mp
from itertools import product
from tqdm import tqdm
from molGrouper.group_graph import GroupGraph

def worker(line, node_types, int_to_node_type, node_int_to_port, verbose):
    return _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose)

def process_nauty_vcolg_output_mp(filename, node_types, verbose=False):
    int_to_node_type = {i: k for i, k in enumerate(node_types.keys())}
    node_int_to_port = {k: {j: k for j, k in enumerate(v)} for i, (k, v) in enumerate(node_types.items())}
    
    with open(filename) as f:
        lines = f.readlines()
        f.close()

    # Create a partial function for the worker to include the other parameters
    from functools import partial
    worker_func = partial(worker, node_types=node_types, int_to_node_type=int_to_node_type, node_int_to_port=node_int_to_port, verbose=verbose)
    
    # Use multiprocessing to process lines in parallel
    num_workers = mp.cpu_count()
    with mp.Pool(processes=num_workers) as pool:
        results = list(tqdm(pool.imap(worker_func, lines), total=len(lines), desc="Line progress"))
    
    # Combine results
    all_group_graphs = [graph for sublist in results for graph in sublist]
    
    return all_group_graphs

def _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose=False):
    group_graphs_list = []

    node_description, edge_description = line.split("  ")
    edge_description = edge_description.split(" ")
    node_description = node_description.split(" ") # graph is in the form of (node_type, edges {color})

    # edges are in the form of (node1, node2)
    edge_list = []
    for i in range(0, len(edge_description), 2):
        edge_list.append((int(edge_description[i]), int(edge_description[i + 1]))) # source, destination

    n_vertices = int(node_description[0])
    n_edges = int(node_description[1])
    colors = node_description[2:]

    gG = GroupGraph(node_types)
    # Add nodes
    for i in range(n_vertices):
        gG.add_node(f"n{i}", int_to_node_type[int(colors[i])])
    # Add generate a mapping of each edge color possible for each edge
    edge_index_to_edge_color = {}
    for i, e in enumerate(edge_list):
        src, dst = e[0], e[1]
        src_ports = node_types[int_to_node_type[int(colors[src])]]
        dst_ports = node_types[int_to_node_type[int(colors[dst])]]
        edge_index_to_edge_color[i] = list(product(src_ports, dst_ports))
    # Enumerate all possible combinations of edge colors
    for port_combo in product(*edge_index_to_edge_color.values()):
        try:
            temp_gG = gG.copy()
            for j, e in enumerate(edge_list):
                port1, port2 = port_combo[j]
                temp_gG.add_edge((f'n{e[0]}', port1),  (f'n{e[1]}', port2))
            group_graphs_list.append(temp_gG)
        except Exception as e:
            if verbose:
                print(f"Couldn't produce graph from multig output: {e}")
    return group_graphs_list
