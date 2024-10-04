from genGrouper import Node, exhaustive_generate
import time
from rdkit import Chem
import argparse
import pickle
import pathlib

if __name__ == "__main__":

    node_defs = set()
    node_defs.add(Node(0, 'methyl', 'C', [0,1,2], [0,0,0]))
    node_defs.add(Node(0, 't2', 'N', [0,1,2], [0,0,0]))
    node_defs.add(Node(0, 'ester', 'C(=O)O', [0,1], [0,2]))
    node_defs.add(Node(0, 'extra1', 'O', [0,1], [0,0]))
    # node_defs.add(Node(0, 'hydroxymethyl', 'CO', [0], [0]))
    # node_defs.add(Node(0, 'primary_amine', 'CN', [0,1,2], [0,0,0]))
    # node_defs.add(Node(0, 'secondary_amine', 'CNC', [0,1], [0,0]))

    # positive_constraints = {"hydroxyl" : 1, "tertiary_amine" : 1}
    # negative_constraints = {'NN', 'NO', 'NCN', 'NCO', 'OCO'}
    positive_constraints = {}
    negative_constraints = set()


    # parse arguments
    parser = argparse.ArgumentParser(description='Exhaustively generate set of molecular graphs')
    parser.add_argument('--n', type=int, default=-1, help='Number of nodes in the graph')
    parser.add_argument('--n_cpus', type=int, default=8, help='Number of cpus to use for multiprocessing')
    parser.add_argument('--config_path', type=str, default="", help='Path to config file')
    args = parser.parse_args()

    parent = str(pathlib.Path(__file__).parent.absolute())

    print(f"Generating all possible molecular graphs with {args.n} nodes\n")
    print(f"Multiprocessing with {args.n_cpus} cpus\n")
    if len(args.config_path) > 0:
        print(f"Saving to database with config at {args.config_path}\n")

    # call nauty
    start = time.time()
    result = exhaustive_generate(
        args.n, 
        node_defs, 
        nauty_path="/raid6/homes/kierannp/projects/molGrouper/packages/nauty2_8_8",
        input_file_path="",
        num_procs=args.n_cpus,
        positive_constraints=positive_constraints,
        negative_constraints=negative_constraints,
        config_path=args.config_path,
        verbose=False
    )
    end = time.time()
    print(f"Time taken for generation: {end - start}")

    print(f"Total graphs: {len(result)}")
