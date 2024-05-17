from molGrouper.group_graph import GroupGraph
from itertools import product
from rdkit import Chem
from pysmiles import write_smiles
from tqdm import tqdm

"""
read in nauty vcolg input 
map color int to node type
from node type number of port determine the number of edge colored needed to be iterated over
save output graph to to datastructure
"""



def process_nauty_vcolg_output(filename, node_types, verbose=False):
    int_to_node_type = {i: k for i, k in enumerate(node_types.keys())}
    node_int_to_port = {k: {j: k for j, k in enumerate(v)} for i, (k,v) in enumerate(node_types.items())}
    all_group_graphs = []
    with open(filename) as f:
        lines = f.readlines()
        f.close()
    for line in tqdm(lines, desc="Line progress"):
        out = _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose)
        if len(out) > 0:
            all_group_graphs.extend(out)
    return all_group_graphs

def _process_nauty_graph_vcolg_output(line, node_types, int_to_node_type, node_int_to_port, verbose=False):
    group_graphs_list = []

    node_description, edge_description = line.split("  ")
    edge_description = edge_description.split(" ")
    node_description = node_description.split(" ") # graph is in the form of (node_type, edges {color})

    # edges are in the form of (node1, node2)
    edge_list = []
    for i in range(0, len(edge_description), 2):
        edge_list.append((int(edge_description[i]), int(edge_description[i+1]))) # source, destination

    n_vertices = int(node_description[0])
    n_edges = int(node_description[1])
    colors = node_description[2:]
    
    gG = GroupGraph(node_types)
    # Add nodes
    for i in range(n_vertices):
        gG.add_node(f"n{i}", int_to_node_type[int(colors[i])])
    # Add generate a mapping of each edge color posible for each edge
    edge_index_to_edge_color = {}
    for i, e in enumerate(edge_list):
        src, dst = e[0], e[1]
        src_ports = node_types[int_to_node_type[int(colors[0])]]
        dst_ports = node_types[int_to_node_type[int(colors[1])]]
        edge_index_to_edge_color[i] = list(product(src_ports, dst_ports))
    # Enumerate all possible combinations of edge colors
    for port_combo in tqdm(product(*edge_index_to_edge_color.values()), leave=False, desc="Edge color progress"):
        try:
            temp_gG = gG.copy()
            for j, e in enumerate(edge_list):
                port1, port2 = port_combo[j]
                temp_gG.add_edge((f'n{e[0]}', port1),  (f'n{e[1]}', port2))
            group_graphs_list.append(temp_gG)
        except:
            if verbose:
                print("Couldn't produce graph from multig output")
    return group_graphs_list





if __name__ == "__main__":
    node_types = {
        'OH': ['O1'],
        'carbon_1111': ['C1', 'C2', 'C3', 'C4'],
        'methine_111': ['C1', 'C2', 'C3'],
        'methylene_11': ['C1', 'C2'],
        'methyl_1': ['C1'],
        'amine0_11': ['N1', 'N2'],
        'amine1_11': ['N1', 'N2'],
        'amine2_11_1': ['N1', 'C1', 'C2'],
        'amine3_1_1': ['N1', 'C1'],
        'amine3': ['N1'],
    }
    node_type_to_smiles = {
        'OH' : 'O',
        'carbon_1111': 'C',
        'methine_111': 'C',
        'methylene_11': 'C',
        'methyl_1': 'C',
        'amine0_11': 'CN',
        'amine1_11': 'CN',
        'amine2_11_1': 'CN',
        'amine3_1_1': 'CN',
        'amine3': 'CN',
    }
    node_port_to_atom_index = {
        'OH' : {'O1': 0},
        'carbon_1111': {'C1': 0, 'C2': 0, 'C3': 0, 'C4': 0},
        'methine_111': {'C1': 0, 'C2': 0, 'C3': 0},
        'methylene_11': {'C1': 0, 'C2': 0},
        'methyl_1': {'C1': 0},
        'amine0_11': {'N1': 1, 'N2': 1},
        'amine1_11': {'N1': 1, 'N2': 1},
        'amine2_11_1': {'N1': 1, 'C1': 0, 'C2': 0},
        'amine3_1_1': {'N1': 1, 'C1': 0},
        'amine3': {'N1': 1},
    }

    out = process_nauty_vcolg_output('/Users/kieran/projects/molGrouper/molGrouper/vcolg_foo.txt', node_types, verbose=False)

    unique_mols = set()
    groupGraphs = []
    for g in out:
        mG = g.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
        if Chem.MolFromSmiles(write_smiles(mG)) is not None:
            canon = Chem.MolToSmiles(Chem.MolFromSmiles(write_smiles(mG)), canonical=True)
            if canon in unique_mols:
                continue
            unique_mols.add(canon)
            groupGraphs.append(g)
        else:
            print("Rdkit failed from conversion between smiles and molecular graph")
            print(g)
            break
    print(f"Unique: {len(unique_mols)}, Total: {len(out)}")