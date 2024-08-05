#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <iostream>

int main() {
    RDKit::RWMol mol;

    // Add atoms (e.g., constructing ethanol C2H6O as an example)
    RDKit::Atom *c1 = new RDKit::Atom(6); // Carbon with atomic number 6
    RDKit::Atom *c2 = new RDKit::Atom(6);
    RDKit::Atom *o = new RDKit::Atom(8);  // Oxygen with atomic number 8

    unsigned int idx1 = mol.addAtom(c1, false, true);
    unsigned int idx2 = mol.addAtom(c2, false, true);
    unsigned int idx3 = mol.addAtom(o, false, true);

    // Add hydrogen atoms for simplicity
    for (int i = 0; i < 3; ++i) {
        RDKit::Atom *h = new RDKit::Atom(1); // Hydrogen with atomic number 1
        mol.addAtom(h, false, true);
        mol.addBond(idx1, mol.getNumAtoms() - 1, RDKit::Bond::SINGLE);
    }
    for (int i = 0; i < 2; ++i) {
        RDKit::Atom *h = new RDKit::Atom(1);
        mol.addAtom(h, false, true);
        mol.addBond(idx2, mol.getNumAtoms() - 1, RDKit::Bond::SINGLE);
    }
    RDKit::Atom *h = new RDKit::Atom(1);
    mol.addAtom(h, false, true);
    mol.addBond(idx3, mol.getNumAtoms() - 1, RDKit::Bond::SINGLE);

    // Add bonds between C-C and C-O
    mol.addBond(idx1, idx2, RDKit::Bond::SINGLE);
    mol.addBond(idx2, idx3, RDKit::Bond::SINGLE);

    // Generate and print the SMILES string
    std::string smiles = RDKit::MolToSmiles(mol);
    std::cout << "SMILES: " << smiles << std::endl;

    return 0;
}

