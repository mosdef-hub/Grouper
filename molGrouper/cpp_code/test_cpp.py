#/usr/bin python3

import molGrouper


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

int_to_node_type = {i: k for i, k in enumerate(node_types.keys())}
node_int_to_port = {k: {j: k for j, k in enumerate(v)} for i, (k,v) in enumerate(node_types.items())}
