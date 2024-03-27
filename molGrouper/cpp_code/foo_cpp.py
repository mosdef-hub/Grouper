#/usr/bin python3

import groupGraphBindings as ggb
from rdkit import Chem
from pysmiles import write_smiles

# # Example usage
# graph = ggb.GroupGraph({"type1": [1, 2, 3], "type2": [4, 5, 6]})
# graph.addNode(1, "type1")
# graph.addNode(2, "type2")
# graph.addEdge(1, 2, 2, 5)
# graph.addEdge(1, 3, 2, 4)

# graph.make_undirected()
# graph.printGraph()

# print("Free Ports for Node 1:", graph.n_free_ports(1))
# print("Free Ports for Node 2:", graph.n_free_ports(2))

node_types = {
    'NH2': ['N1'], # amine
    'CO': ['C1', 'C2'], # carbonyl
    'CC': ['C11', 'C12', 'C21', 'C22'], # alkene
}
# node_types = {"NH2": [0], "CO": [0, 1], "CC": [0, 1, 2, 3]}
# int_to_node_type = {0: "NH2", 1: "CO", 2: "CC"}
# node_int_to_port = {'NH2': [0], 'CO': [0, 1], 'CC': [0, 1, 2, 3]}

# # Example usage
# ggb.MultigConverter.multig_line_to_graph(
#     line = "3 3 0 1 0 2 1 2",
#     int_to_node_type = int_to_node_type,
#     node_int_to_port = node_int_to_port,
#     node_types = node_types,
#     verbose=False
# )

# int_to_node_type = {i: k for i, k in enumerate(node_types.keys())}
# node_int_to_port = {k: {j: k for j, k in enumerate(v)} for i, (k,v) in enumerate(node_types.items())}

node_types = {
    'NH2': [0], # amine
    'CO': [0, 1], # carbonyl
    'CC': [0, 1, 2, 3], # alkene
}
node_type_to_smiles = {
    'NH2': 'N',
    'CO': 'C=O',
    'CC': 'C=C',
}
node_port_to_atom_index = {
    'NH2': {'N1': 0},
    'CO': {'C1': 0, 'C2': 0},
    'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
}


result = ggb.MultigConverter.parse_multig_file()



from molGrouper.group_graph import GroupGraph


def cpp_to_py(cpp_groupG):
    py_groupG = GroupGraph(node_types)
    for k, n in cpp_groupG.nodes.items():
        py_groupG.add_node(k, n.type)
    for e in cpp_groupG.edges:
        py_groupG.add_edge(e[0], e[1], e[2], e[3])
    return py_groupG

unique_mols = set()
groupGraphs = []
for g in result:
    if 'placeholder1' in g.node_types:
        continue
    new_g = cpp_to_py(g)
    mG = new_g.to_molecular_graph(node_type_to_smiles, node_port_to_atom_index)
    if Chem.MolFromSmiles(write_smiles(mG)) is not None:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(write_smiles(mG)), canonical=True)
        if canon in unique_mols:
            continue
        unique_mols.add(canon)
        groupGraphs.append(g)
print(len(unique_mols))
print(len(result))

