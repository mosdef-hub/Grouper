from molGrouper.edge_coloring_mp import process_nauty_vcolg_mp
from molGrouper.generate import _call_geng, _call_vcolg
import time
from rdkit import Chem
from pysmiles import write_smiles
import argparse
import pickle
import pathlib

if __name__ == "__main__":
    node_types = {
        'carbon': ['C1', 'C2', 'C3', 'C4'],
        'methine': ['C1', 'C2', 'C3'],
        'methylene': ['C1', 'C2'],
        'methyl': ['C1'],
        'hydroxymethyl': ['C1'],
        'amine': ['N1'],
    }
    node_type_to_smiles = {
        'carbon': 'C',
        'methine': 'C',
        'methylene': 'C',
        'methyl': 'C',
        'hydroxymethyl': 'CO',
        'amine': 'CN',
    }
    node_port_to_atom_index = {
        'carbon': {'C1': 0, 'C2': 0, 'C3': 0, 'C4': 0},
        'methine': {'C1': 0, 'C2': 0, 'C3': 0},
        'methylene': {'C1': 0, 'C2': 0},
        'methyl': {'C1': 0},
        'hydroxymethyl': {'C1': 0},
        'amine': {'N1': 1},
    }
    # parse arguments
    parser = argparse.ArgumentParser(description='Exhaustively generate set of molecular graphs')
    parser.add_argument('--n', type=int, default=3, help='Number of nodes in the graph')
    parser.add_argument('--n_cpus', type=int, default=-1, help='Number of cpus to use for multiprocessing')
    args = parser.parse_args()

    parent = str(pathlib.Path(__file__).parent.absolute())

    print(f"Generating all possible molecular graphs with {args.n} nodes\n")
    print(f"Multiprocessing with {args.n_cpus} cpus\n")

    # call nauty
    start = time.time()
    print(args.n, len(max(node_types.values(), key=len)))
    _call_geng(parent_path=parent, n_nodes = args.n, max_edges = len(max(node_types.values(), key=len)))
    _call_vcolg(parent_path=parent, n_colors=len(node_types))
    end = time.time()
    print(f"Time taken for nauty: {end - start}")

    # process nauty output
    start = time.time()
    process_nauty_vcolg_mp(f'{parent}/vcolg_out.txt', node_types, verbose=False, n_processes=args.n_cpus, smiles_hash=True, node_type_to_smiles=node_type_to_smiles, node_port_to_atom_index=node_port_to_atom_index)
    end = time.time()
    print(f"Time taken for process_nauty_vcolg__mp: {end - start}")
    