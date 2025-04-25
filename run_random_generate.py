from Grouper import Group, random_generate, GroupGraph
from Grouper.visualization import visualize
from Grouper.utils import convert_to_nx
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
    node_defs.add(Group('carbon', 'C', [0,0,0,0]))
    node_defs.add(Group('hydroxymethyl', 'CO', [0,0,0]))
    node_defs.add(Group('secondary_amine', 'CNC', [0,0,0,1,2,2,2]))
    node_defs.add(Group('amine', 'N', [0,0,0]))
    node_defs.add(Group('hydroxyl', 'O',  [0]))

    # positive_constraints = {"hydroxyl" : 1, "tertiary_amine" : 1}
    # negative_constraints = {'NN', 'NO', 'NCN', 'NCO', 'OCO'}
    positive_constraints = {}
    negative_constraints = set()


    # parse arguments
    parser = argparse.ArgumentParser(description='Exhaustively generate set of molecular graphs')
    parser.add_argument('--n', type=int, default=3, help='Number of nodes in the graph')
    parser.add_argument('--n_graphs', type=int, default=100, help='Number of graphs to generate')
    parser.add_argument('--n_cpus', type=int, default=1, help='Number of cpus to use for multiprocessing')
    args = parser.parse_args()

    parent = str(pathlib.Path(__file__).parent.absolute())

    # call nauty
    start = time.time()
    result = random_generate(
        args.n, 
        node_defs, 
        args.n_graphs,
        num_procs=args.n_cpus,
        nauty_path="/raid6/homes/kierannp/projects/nauty2_8_9",
        positive_constraints=positive_constraints,
        negative_constraints=negative_constraints,
    )
    end = time.time()
    print(f"Time taken for generation: {end - start}")

    print(f"Total graphs: {len(result)}")


    # save image of each molecule  
    mols = []
    for i, graph in enumerate(result):
        mol = Chem.MolFromSmiles(graph.to_smiles())
        mols.append(mol)
    Chem.Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200,200)).save(f"{parent}/random_molecules_{args.n}.png")

    # Function to convert matplotlib figure to PIL Image
    def fig_to_img(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight')
        buf.seek(0)
        return Image.open(buf)
    
    # save grid image of each graph
    figures  = []
    for graph in result:
        fig = visualize(graph)
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
    composite.save(f'random_graphs_{args.n}.png')


    # # Save to pickle
    with open(f"{parent}/random_graphs_{args.n}.pkl", 'wb') as f:
        pickle.dump(result, f)
