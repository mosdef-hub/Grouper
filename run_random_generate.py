import argparse
import pathlib
import pickle
import time

from Grouper import Group, random_generate

if __name__ == "__main__":
    node_defs = set()
    node_defs.add(Group("carbon", "C", [0, 0, 0, 0]))
    node_defs.add(Group("hydroxymethyl", "CO", [0, 0, 0]))
    node_defs.add(Group("secondary_amine", "CNC", [0, 0, 0, 1, 2, 2, 2]))
    node_defs.add(Group("amine", "N", [0, 0, 0]))
    node_defs.add(Group("hydroxyl", "O", [0]))

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
        "--n_graphs", type=int, default=100, help="Number of graphs to generate"
    )
    parser.add_argument(
        "--n_cpus",
        type=int,
        default=1,
        help="Number of cpus to use for multiprocessing",
    )
    args = parser.parse_args()

    parent = str(pathlib.Path(__file__).parent.absolute())

    # call nauty
    start = time.time()
    result = random_generate(
        args.n,
        node_defs,
        args.n_graphs,
        num_procs=args.n_cpus,
        positive_constraints=positive_constraints,
        negative_constraints=negative_constraints,
    )
    end = time.time()
    print(f"Time taken for generation: {end - start}")

    print(f"Total graphs: {len(result)}")

    # # Save to pickle
    with open(f"{parent}/random_graphs_{args.n}.pkl", "wb") as f:
        pickle.dump(result, f)
