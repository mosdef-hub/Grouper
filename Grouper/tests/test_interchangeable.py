from Grouper import GroupGraph, Group, fragment, AtomGraph
from Grouper.tests.base_test import BaseTest

# from pysmiles import write_smiles
# import pytest


class TestConversion(BaseTest):
    def test_smiles_gg_to_smiles(self):
        test_cases = [
            "O=COOC(=O)O",
            "CNCNOC=O",
            "O=C(O)OC(=O)O",
            "NNO",
            "O=C(O)C(=O)OO",
            "O=COC(=O)OO",
            "CNO",
            "CCO",
            "COC(N)=O",
            "ONO",
            "CNCC(=O)OC",
            "CNC(N)C(=O)O",
            "CC(=O)ON",
            "CNCOC(C)=O",
            "CNC(O)O",
            "CNCOC(=O)CNC",
            "CNC(=O)O",
            "CNOC=O",
            "NCN",
            "CNCC(NC)OC=O",
            "CNCCN",
            "O=CONC(=O)O",
            "CCN",
            "CNCOC(=O)OC=O",
            "CNCC(=O)ON",
            "CNN",
            "O=CONOC=O",
            "CNC(C)O",
            "CNCC(N)NC",
            "CNC(OC=O)OC=O",
            "CNC(OC=O)C(=O)O",
            "CNCNN",
            "CNCOC(N)=O",
            "CCC",
            "NOC(N)=O",
            "CNC(N)N",
            "NCC(=O)O",
            "CNCNC",
            "CNC(N)O",
            "CNCNCNC",
            "OCO",
            "CNCC(NC)C(=O)O",
            "CNCCO",
            "CNCC(=O)OOC=O",
            "CNC(O)OC=O",
            "O=COOC(=O)OC=O",
            "CNCOC(=O)C(=O)O",
            "CNC(O)C(=O)O",
            "CC(=O)OC(=O)O",
            "CC(=O)OOC=O",
            "O=C(O)CC(=O)O",
            "NC(=O)OOC=O",
            "COC(=O)C(=O)O",
            "COC(=O)OC=O",
            "O=COC(=O)OC(=O)O",
            "O=COOC(=O)C(=O)O",
            "NOC(=O)OC=O",
            "NNC(=O)O",
            "CNCC(=O)OO",
            "CNC(C(=O)O)C(=O)O",
            "CNCC(O)NC",
            "CNCOC(=O)O",
            "COC(C)=O",
            "CNCCOC=O",
            "NCOC=O",
            "CCCNC",
            "O=CONO",
            "CNC(C)C",
            "O=C(O)NC(=O)O",
            "O=C(O)OO",
            "O=C(O)CO",
            "O=COCO",
            "O=C(O)NO",
            "CNCC(=O)OC(=O)O",
            "NC(=O)OC(=O)O",
            "NNN",
            "NOC(=O)C(=O)O",
            "NCO",
            "CNCNC(=O)O",
            "CNC",
            "O=COCOC=O",
            "O=COCC(=O)O",
            "NNOC=O",
            "CNC(C)C(=O)O",
            "CNC(N)OC=O",
            "CNC(C)OC=O",
            "CC(=O)OO",
            "CNCNO",
            "CNC(C)N",
            "CNCC(C)NC",
            "COC(=O)O",
            "O=C(O)OC(=O)C(=O)O",
            "CCC(=O)O",
            "CNCC(CNC)NC",
            "CNCCCNC",
            "CCOC=O",
            "CNCCC(=O)O",
            "NC(=O)OO",
            "NOC(=O)O",
            # 7 nodes in the graph
            "CNCC(=O)ON(C)COC(N)=O",
            "CC(CN(C)C)C(=O)O",
            "CC(CN(C)C)OC=O",
            "O=C(O)C(O)N(O)C(=O)OO",
            "CNCCC(CC(NC)NO)NC",
            "CNC(CO)NNNN",
            "CNCCN(N)C(=O)OC(=O)ON",
            "CNC1OC(=O)C(NC)N(C(=O)O)C1O",
            "CNCC(=O)ON1C(=O)OC(=O)ONC1NC",
            "CNCC(NC)C(CNC)C(NC)C(NC)C(=O)O",
            "CNC(CC(=O)O)C(NC)N(O)OC=O",
            "CNC(CC(=O)O)C(NC)N(O)C(=O)O",
            "CCN(OC=O)C(N)CNC",
            "CNC(NC(NC)C(=O)OO)C(=O)ON",
            "CNCOC(=O)OC(=O)C(CNN)NC",
            "CNCC(=O)OOC(=O)C(NC)NC(O)NC",
            "CNCC(=O)OCNCNC(=O)O",
            "CNCN(NC)NOC(=O)O",
            "CNCN(O)NOC(=O)NO",
            "CNC(CC(=O)O)OC(=O)C(=O)ONOC=O",
            # "CNC1CO(C=O)NNC(=O)OOC1=O", this is actually an invalid SMILES, somehow it got generated though
            "CNC1COC(=O)NNC(=O)OOC1=O",
            "CNCOC(=O)CNC(O)O",
            "CNC1CCC(NC)C1CO",
            "NC1COC(=O)CC(=O)OC1",
            "CC(=O)ONNOC(=O)CC(=O)O",
            "CNCC(=O)OOC(=O)NC(=O)ONO",
            "CN1CNC(N)OC1=O",
            "CNCCC(CCOC=O)CNC",
            "CNCC(NC)C(NC)C(N)C(=O)OOC=O",
            "CNCC(NC)N(OC=O)C(NC)C(=O)OO",
            "COC(=O)CCCNO",
            "CNCOC(=O)C(C)NC(=O)OC",
            "CNCC(NC)N(C(=O)O)C(NC)C(=O)OO",
            "CNCC(CNNOC(=O)OC=O)NC",
            "CNCN(C(=O)OC)C(=O)OOC(C)=O",
            "CN1COC(=O)COC(=O)OC1=O",
        ]
        # test_cases = [
        # #     'CNC(=O)O' # C - N - C(=O)-O
        # 'CNC1CO(C=O)NNC(=O)OOC1=O'
        # ] 
        node_defs = set()
        node_defs.add(Group('hydroxyl', 'O', [0]))
        node_defs.add(Group('keto', 'O', [0,0]))
        node_defs.add(Group('ester', 'C(=O)O', [0,2]))   # - C=OO, -C=OO-,-C=OO-, C=OO-
        node_defs.add(Group('methyl', 'C', [0,0,0]))
        node_defs.add(Group('t2', 'N', [0,0,0]))
        node_defs.add(Group('secondary_amine', 'CNC', [0,0,1,2,2,2]))

        for smiles in test_cases:
            all_gGs = fragment(smiles, node_defs)
            if len(all_gGs) == 0:
                raise ValueError(f"Failed to fragment {smiles}")
            for gg in all_gGs:
                s = gg.to_smiles()
                assert smiles == s
                
    def test_aromatic_smiles(self):
        smiles = "c1ccccc1"
        node_defs = set()
        node_defs.add(Group('ringC', '[cX3]', [0,0,0], "SMARTS"))
        node_defs.add(Group('COOH', '[CX3](=[OX1])[O;H1]', [0], "SMARTS"))
        gg = fragment(smiles, node_defs)[0]
        assert gg.to_smiles() == smiles
