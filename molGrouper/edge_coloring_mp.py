import multiprocessing as mp
from itertools import product
from tqdm import tqdm
import networkx as nx
import rdkit.Chem
from pysmiles import write_smiles
import datetime
import pickle
import pathlib
from glob import glob
import os
import time

from molGrouper.group_graph import GroupGraph
from molGrouper.utils import convert_edges_to_nodetype


# def worker(line, node_types, int_to_node_type, node_int_to_port, verbose):
#     return _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose)

def node_match(n1, n2):
    return n1['type'] == n2['type']

def edge_match(e1, e2):
    return set(e1['ports']) == set(e2['ports'])

def worker(chunk, node_types, int_to_node_type, node_int_to_port, verbose=False, smiles_hash=True, node_type_to_smiles=None, node_port_to_atom_index=None):
    if smiles_hash:
        assert node_type_to_smiles is not None
        assert node_port_to_atom_index is not None
    groupGraphs = []
    basis = set()
    for i, line in tqdm(enumerate(chunk), desc="Chunk progress", total=len(chunk)):
        edge_colored_graphs = _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose)
        if verbose:
            print(f"Lots of colored edge graphs found: {len(edge_colored_graphs)}")
            print(f"Total unique graphs: {len(basis)}")
        if len(edge_colored_graphs) == 0:
            continue
        for g in edge_colored_graphs:
            # Convert graph to a hash either using WL or SMILES
            if smiles_hash:
                try:
                    mG = g.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
                    smiles = write_smiles(mG)
                    mol = rdkit.Chem.MolFromSmiles(smiles)
                    g_hash = rdkit.Chem.MolToSmiles(mol, canonical=True)
                except Exception as e:
                    continue
            else:
                g_hash = nx.weisfeiler_lehman_graph_hash(g, edge_attr="ports", node_attr="type")  

            if g_hash in basis: # TODO this check if the graphs are isomorphic, but not automorphic. Need to check for automorphism, it is possible that the graph is isomorphic but actually represents a different molecule
                continue
            basis.add(g_hash)
            groupGraphs.append(g) # if non-isomorphic, add to the list of unique molecular graphs
        if verbose:
            if i % 100 == 0 and i != 0:
                print(f"Processed {i} lines")

    now = datetime.datetime.now()
    microtime = now.strftime("%f")
    parent = str(pathlib.Path(__file__).parent.absolute())
    with open(f'{parent}/../group_graph_{microtime}.pkl', 'wb') as f:
        pickle.dump(groupGraphs, f)
    
    print("Finished processing chunk at", now)

def process_nauty_vcolg_mp(filename, node_types, n_processes = -1, verbose=False, smiles_hash = False, node_type_to_smiles=None, node_port_to_atom_index=None):
    int_to_node_type = {i: k for i, k in enumerate(node_types.keys())}
    node_int_to_port = {k: {j: k for j, k in enumerate(v)} for i, (k, v) in enumerate(node_types.items())}

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

    # A rule of thumb is we don't want our chunk_size to be larger than 1000
    chunk_size = min(len(lines) // (num_workers), 1000) 
    # Create chunks of lines
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]

    worker_func = partial(
        worker, 
        node_types=node_types, 
        int_to_node_type=int_to_node_type, 
        node_int_to_port=node_int_to_port, 
        verbose=verbose,
        smiles_hash=smiles_hash,
        node_type_to_smiles=node_type_to_smiles,
        node_port_to_atom_index=node_port_to_atom_index
        )
    
    # Use multiprocessing to process lines in parallel
    start = time.time()
    with mp.Pool(processes=num_workers) as pool:
        try:
            # results = list(tqdm(pool.imap(worker_func, lines), total=len(lines), desc="Line progress"))
            tqdm(pool.imap(worker_func, chunks), total=len(chunks), desc="Worker progress")
        except Exception as e:
            print(f"Error in multiprocessing: {e}")
            pool.terminate()
        finally:
            pool.close()
            pool.join()
    end = time.time()
    print(f"Time taken in parallel: {end - start}")
    # Process text outputs into single file
    glob_files = glob(f"{pathlib.Path(__file__).parent.absolute()}/../group_graph_*.pkl")
    groupGraphs = []
    basis = set()

    for file in glob_files:
        with open(file, 'rb') as f:
            lines = pickle.load(f)
        for g in lines:
            if smiles_hash:
                try:
                    mG = g.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
                    smiles = write_smiles(mG)
                    mol = rdkit.Chem.MolFromSmiles(smiles)
                    g_hash = rdkit.Chem.MolToSmiles(mol, canonical=True)
                except Exception as e:
                    continue
            else:
                g_hash = nx.weisfeiler_lehman_graph_hash(g, edge_attr="ports", node_attr="type")

            if g_hash in basis:
                continue
            basis.add(g_hash)
            groupGraphs.append(g) # if non-isomorphic, add to the list of unique molecular graphs

    with open(f"{pathlib.Path(__file__).parent.absolute()}/../basis.pkl", 'wb') as f:
        pickle.dump(groupGraphs, f)

    print(f"{len(groupGraphs)} unique group graphs found")
    
    # Process smiles
    if smiles_hash:
        start = time.time()
        unique_smiles = set()
        for g in groupGraphs:
            mG = g.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
            try:
                smiles = write_smiles(mG)
                unique_smiles.add(smiles)
            except Exception as e:
                continue
            
        with open(f"{pathlib.Path(__file__).parent.absolute()}/../basis_smiles.txt", "w") as f:
            for g in unique_smiles:
                f.write(g + "\n")
        end = time.time()
        print(f"Time taken to convert to smiles: {end - start}")

    # Remove nauty files
    os.remove(f"{pathlib.Path(__file__).parent.absolute()}/../geng_out.txt")
    os.remove(f"{pathlib.Path(__file__).parent.absolute()}/../vcolg_out.txt")


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
