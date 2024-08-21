from genGrouper import Node, exhaustive_generate
import time
from rdkit import Chem
import argparse
import pickle
import pathlib

if __name__ == "__main__":

    node_defs = set()
    node_defs.add(Node(0, 'carbon', 'C', [0,1,2,3], [0,0,0,0]))
    node_defs.add(Node(0, 'methine', 'C', [0,1,2], [0,0,0]))
    node_defs.add(Node(0, 'methylene', 'C', [0,1], [0,0]))
    node_defs.add(Node(0, 'methyl', 'C', [0], [0]))
    node_defs.add(Node(0, 'hydroxymethyl', 'CO', [0], [0]))
    node_defs.add(Node(0, 'primary_amine', 'CN', [0,1,2], [0,0,0]))
    node_defs.add(Node(0, 'secondary_amine', 'CNC', [0,1], [0,0]))


    # parse arguments
    parser = argparse.ArgumentParser(description='Exhaustively generate set of molecular graphs')
    parser.add_argument('--n', type=int, default=-1, help='Number of nodes in the graph')
    parser.add_argument('--n_cpus', type=int, default=8, help='Number of cpus to use for multiprocessing')
    args = parser.parse_args()

    parent = str(pathlib.Path(__file__).parent.absolute())

    print(f"Generating all possible molecular graphs with {args.n} nodes\n")
    print(f"Multiprocessing with {args.n_cpus} cpus\n")

    # call nauty
    start = time.time()
    result = exhaustive_generate(
        args.n, 
        node_defs, 
        nauty_path="/raid6/homes/kierannp/projects/molGrouper/packages/nauty2_8_8",
        input_file_path="",
        num_procs=args.n_cpus,
        verbose=False
    )
    end = time.time()
    print(f"Time taken for generation: {end - start}")

    print(f"Total graphs: {len(result)}")

    # save set of unique mols
    with open("unique_mols.txt", "w") as f:
        for mol in result:
            f.write(f"{mol}\n")
