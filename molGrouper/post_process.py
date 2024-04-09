

def substitute_chiral_smiles(smiles, original_sub_smiles, subs):
    new_smiles = []
    for s in subs:
        new_smiles.append(smiles.replace(original_sub_smiles, s))
    return new_smiles