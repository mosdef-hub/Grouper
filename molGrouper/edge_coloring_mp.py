import multiprocessing as mp
from itertools import product
from tqdm import tqdm
from molGrouper.group_graph import GroupGraph
import networkx as nx
from rdkit import Chem
from pysmiles import write_smiles

# def worker(line, node_types, int_to_node_type, node_int_to_port, verbose):
#     return _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose)
def worker(chunk, node_types, int_to_node_type, node_int_to_port, verbose, just_smiles=False, node_type_to_smiles=None, node_port_to_atom_index=None):
    results = []
    if not just_smiles:
        for i, line in enumerate(chunk):
            results.extend(_process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose))
            if i % 100 == 0:
                print(f"Processed {i} lines")
        
    else:
        assert node_type_to_smiles is not None
        assert node_port_to_atom_index is not None
        unique_mols = set()
        groupGraphs = []
        for i, line in enumerate(chunk):
            edge_colored_graphs = _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose)
            for g in edge_colored_graphs:
                mG = g.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
                if Chem.MolFromSmiles(write_smiles(mG)) is not None:
                    canon = Chem.MolToSmiles(Chem.MolFromSmiles(write_smiles(mG)), canonical=True)
                    if canon in unique_mols:
                        continue
                    unique_mols.add(canon)
                    groupGraphs.append(g)
                else:
                    if verbose:
                        print("Rdkit failed from conversion between smiles and molecular graph")
                        print(g)
                    break
            if i % 100 == 0:
                print(f"Processed {i} lines")
        
        results = groupGraphs
        
    print("Finished processing chunk")
        
    return results

def process_nauty_vcolg_mp(filename, node_types, n_processes = -1, verbose=False, just_smiles=False, node_type_to_smiles=None, node_port_to_atom_index=None):
    int_to_node_type = {i: k for i, k in enumerate(node_types.keys())}
    node_int_to_port = {k: {j: k for j, k in enumerate(v)} for i, (k, v) in enumerate(node_types.items())}

    if just_smiles:
        assert node_type_to_smiles is not None
        assert node_port_to_atom_index is not None
    
    if n_processes == -1:
        num_workers = mp.cpu_count()
    else:
        num_workers = n_processes

    with open(filename) as f:
        lines = f.readlines()

    # lines = lines[:30000]

    # Create a partial function for the worker to include the other parameters
    from functools import partial

    chunk_size = len(lines) // num_workers
        # Create chunks of lines
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]

    worker_func = partial(
        worker, 
        node_types=node_types, 
        int_to_node_type=int_to_node_type, 
        node_int_to_port=node_int_to_port, 
        verbose=verbose,
        just_smiles=just_smiles,
        node_type_to_smiles=node_type_to_smiles,
        node_port_to_atom_index=node_port_to_atom_index
        )
    
    # Use multiprocessing to process lines in parallel

    with mp.Pool(processes=num_workers, maxtasksperchild=1000) as pool:
        try:
            # results = list(tqdm(pool.imap(worker_func, lines), total=len(lines), desc="Line progress"))
            results = list(tqdm(pool.imap(worker_func, chunks), total=len(chunks), desc="Chunk progress"))
        except Exception as e:
            print(f"Error in multiprocessing: {e}")
            pool.terminate()
        finally:
            pool.close()
            pool.join()
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

    non_colored_edge_list = []

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
        non_colored_edge_list.append((src, dst))
        edge_index_to_edge_color[i] = list(product(src_ports, dst_ports))
    #if the degree of the node is larger than the number of ports, then discard this graph
    G = nx.Graph(non_colored_edge_list) # README NEED TO CONFIRM THIS WORKS
    for n, degree in G.degree():
        if degree > len(node_types[int_to_node_type[int(colors[n])]]):
            return []
    #TODO only generate edge colors that account for the symmetries of the node_type,
    #  for example an edge between benzene has the symmetryies of the dihedral group D6 (6 rotations and 6 reflections)
    # therefore if there are only 3 molecules with 2 attached to benzene then there are 3 possible edge colors 
    # (3 apart, 2 apart, and 1 apart)
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
