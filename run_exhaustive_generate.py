import argparse
import pathlib
import time

from Grouper import Group, exhaustive_generate

if __name__ == "__main__":
    node_defs = set()

    node_defs.add(Group("t2", "N", [0, 0, 0]))
    node_defs.add(Group("Methyl", "C", [0, 0, 0]))
    node_defs.add(Group("ester", "C(=O)O", [0, 2]))
    node_defs.add(Group("extra1", "O", [0, 0]))

    # positive_constraints = {"hydroxyl" : 1, "tertiary_amine" : 1}
    # negative_constraints = {'NN', 'NO', 'NCN', 'NCO', 'OCO'}
    positive_constraints = {}
    negative_constraints = set()

    # parse arguments
    parser = argparse.ArgumentParser(
        description="Exhaustively generate set of molecular graphs"
    )
    parser.add_argument("--n", type=int, default=3, help="Number of nodes in the graph")
    parser.add_argument(
        "--n_cpus",
        type=int,
        default=8,
        help="Number of cpus to use for multiprocessing",
    )
    parser.add_argument(
        "--config_path", type=str, default="", help="Path to config file"
    )
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
        num_procs=args.n_cpus,
        vcolg_output_file="",
        positive_constraints=positive_constraints,
        negative_constraints=negative_constraints,
        config_path=args.config_path,
    )
    end = time.time()
    print(f"Time taken for generation: {end - start}")
    print(f"Total graphs: {len(result)}")
