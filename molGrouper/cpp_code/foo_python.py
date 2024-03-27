#/usr/bin python3

# import groupGraphBindings as ggb
from molGrouper.generate import  _process_multig_output, generate_group_graph_space, _call_geng, _call_vcolg, _call_multig
from rdkit import Chem, RDLogger
from pysmiles import write_smiles


RDLogger.DisableLog('rdApp.*')


node_types = {
    'N': ['N1'], # amine
    'CO': ['C1', 'C2'], # carbonyl
    'CC': ['C11', 'C12', 'C21', 'C22'], # alkene
    # 'C': ['C1'], # alkane
    # 'OH': ['O1'], # alcohol
}
node_type_to_smiles = {
    'N': 'N',
    'CO': 'C=O',
    'CC': 'C=C',
    # 'C': 'C',
    # 'OH': 'O'
}
node_port_to_atom_index = {
    'N': {'N1': 0},
    'CO': {'C1': 0, 'C2': 0},
    'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
    # 'C': {'C1': 0, 'C2': 0, 'C3': 0, 'C4': 0},
    # 'OH': {'O1': 0}
}

int_to_node_type = {i: k for i, k in enumerate(node_types.keys())}
node_int_to_port = {k: {j: k for j, k in enumerate(v)} for i, (k,v) in enumerate(node_types.items())}


out = generate_group_graph_space(5, node_types)

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
print(len(unique_mols))
