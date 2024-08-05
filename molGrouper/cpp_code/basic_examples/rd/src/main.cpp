#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

int main() {
    // Example: create a benzene molecule
    RWMol *benzene = SmilesToMol("c1ccccc1");

    if (benzene) {
        // Convert the molecule to a SMILES string
        std::string smiles = MolToSmiles(*benzene);
        std::cout << "SMILES: " << smiles << std::endl;

        // Clean up
        delete benzene;
    } else {
        std::cerr << "Failed to create molecule from SMILES" << std::endl;
    }

    return 0;
}

