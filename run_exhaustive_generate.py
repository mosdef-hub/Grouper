from genGrouper import Node, exhaustive_generate
from genGrouper.visualization import visualize
from genGrouper.utils import convert_to_nx
import time
from rdkit import Chem
from rdkit.Chem import Draw
import argparse
import pickle
import pathlib
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
import io

if __name__ == "__main__":

    node_defs = set()
    node_defs.add(Node(0, 'carbon', 'C', [0,0,0,0]))
    # node_defs.add(Node(0, 'methine', 'C', [0,1,2], [0,0,0]))
    # node_defs.add(Node(0, 'methylene', 'C', [0,1], [0,0]))
    # node_defs.add(Node(0, 'methyl', 'C', [0], [0]))
    # node_defs.add(Node(0, 'hydroxymethyl', 'CO', [0], [0]))
    # node_defs.add(Node(0, 'primary_amine', 'CN', [0,1,2], [0,0,0]))
    # node_defs.add(Node(0, 'secondary_amine', 'CNC', [0,1], [0,0]))
    node_defs.add(Node(0, 'tertiary_amine', 'N', [0,0,0]))
    node_defs.add(Node(0, 'hydroxyl', 'O',  [0]))

    # positive_constraints = {"hydroxyl" : 1, "tertiary_amine" : 1}
    # negative_constraints = {'NN', 'NO', 'NCN', 'NCO', 'OCO'}
    positive_constraints = {}
    negative_constraints = set()


    # parse arguments
    parser = argparse.ArgumentParser(description='Exhaustively generate set of molecular graphs')
    parser.add_argument('--n', type=int, default=3, help='Number of nodes in the graph')
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
        nauty_path="/Users/kieran/projects/genGrouper/packages/nauty2_8_8",
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


    # save image of each molecule  
    mols = []
    for i, graph in enumerate(result):
        mol = Chem.MolFromSmiles(graph.to_smiles())
        mols.append(mol)
    Chem.Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200,200)).save(f"{parent}/molecules_{args.n}.png")

    # Function to convert matplotlib figure to PIL Image
    def fig_to_img(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        return Image.open(buf)
    
    # save grid image of each graph
    figures  = []
    for graph in result:
        nx_graph = convert_to_nx(graph)
        fig = visualize(nx_graph)
        figures.append(fig_to_img(fig))
        plt.close(fig)

    # Arrange the figures in a grid 5 per row
    fig_width, fig_height = figures[0].size  # Assume all figures have the same size
    composite = Image.new('RGB', (fig_width * 5, fig_height * ((len(result) // 5) + 1)))  # Create blank canvas

    for i, figure in enumerate(figures):
        x = (i % 5) * fig_width
        y = (i // 5) * fig_height
        composite.paste(figure, (x, y))

    # Show or save the composite image
    composite.save(f'graphs_{args.n}.png')






