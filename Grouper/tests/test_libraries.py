import pytest
from networkx import is_connected

from Grouper import Group, fragment
from Grouper.libraries.Libraries import (
    BasisSet,
    GroupExtension,
    Joback,
    SaftGammaMie,
    Unifac,
)
from Grouper.utils import convert_to_nx
from Grouper.tests.base_test import BaseTest

Libraries = {
    "base": BasisSet,
    "saftgm": SaftGammaMie,
    "joback": Joback,
    "UNIFAC": Unifac,
}


class TestLibraries(BaseTest):
    @pytest.mark.parametrize(
        "library,n_graphs",
        [("saftgm", 1), ("joback", 41), ("UNIFAC", 72), ("base", 0)],
    )
    def test_build_library(self, library, n_graphs):
        library = Libraries[library]()
        assert library.n_nodes == n_graphs

    def test_node_trace(self):
        library = Libraries["saftgm"]()
        node = Group("-CH3", "[CH3]", [0], "SMARTS")
        nt = library.query_nodes({"node": node})[0]

        assert nt.node == node
        assert nt.extended_smarts == "[CX4H3]"
        assert nt.doi == "https://doi.org/10.1080/00268978800101601"
        assert nt.priority is None

        nt = GroupExtension(node, "", "[CX4H3]", 1)
        assert nt.priority == 1

    def test_add_node(self):
        library = BasisSet()
        library.add_node(Group("-CH3", "[CH3]", [0], "SMARTS"), "", "[CX4H3]", None)
        assert library.n_nodes == 1

    def test_query_node(self):
        library = Libraries["joback"]()
        nt = library.query_nodes(
            {"extended_smarts": "[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]"}
        )[0]
        assert nt.node == Group("=C<", "C", [0, 0, 0])

    def test_list_nodes(self):
        library = Libraries["joback"]()
        nodes = list(library.get_nodes())
        assert len(nodes) == 41 == library.n_nodes

class TestLibrariesFragmentations(BaseTest):  
    def assert_equal_edgeDicts(self, edge1, edge2):
        for e1, edge_count in edge1.items():
            if edge2.get(e1) == edge_count:
                del edge2[e1]
                continue
            elif edge2.get((e1[1], e1[0])) == edge_count:
                del edge2[(e1[1], e1[0])]
                continue
            print(f"Extra edge {e1} in solution was missed.")
            return False 
        if edge2: # if any missed in edge2
            print(f"Generated extra edge {edge2=} in fragmented graph.")
            return False
        return True
        
    test_molecules = [
        "CC", "CC(C)C", "CCC(C)(C)C", "CC(C)(C)C",# alkanes
        "CO", "OCC(O)CO", # alcohols
        "C(=O)O", "CC(=O)OC", "OCCCC(=O)O", "C(=O)N", "C(=O)N(C)C", # carboxyls/amides
        "CN", "CNC", "CN(C)CN"# amines
        "C[O-]", "C(=O)[O-]", "CCS(=O)(=O)[O-]", "CCC#N", "C(=O)[N-]",# acids
        "C(F)(F)(F)F", "CCCCF",# fluorines
        "C1CCCCC1C","C1CCCC(C)C1C", "C1CCC(C)CC1C", "C1CC(C)CCC1C",# rings
        "CC1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]", "C1=CC=C(C=C1)C(=O)O",# phenyls
    ]
    nodesList = [{"CH3":2}, {"CH":1, "CH3":3}, {"CH3":4, "CH2":1, "C":1}, {"CH3":4, "C":1}]
    edgesList = [
        {("CH3","CH3"):1}, {("CH3","CH"):3}, {("CH3","CH2"):1, ("CH3","C"):3, ("CH2","C"):1},
        {("C","CH3"):4}
    ]

    @pytest.fixture(scope="session")
    def unifac_groups(self):
        unifac = Unifac()
        return list(unifac.get_nodes())
    
    @pytest.mark.parametrize("MOLSMILES,SITESDICT,EDGESDICT", zip(test_molecules, nodesList, edgesList))
    def test_unifac_smiles(self, MOLSMILES, SITESDICT, EDGESDICT, unifac_groups):
        from collections import Counter
        # molecule = test_molecules[0]
        from Grouper import fragment
        gGList = fragment(MOLSMILES, unifac_groups[:10], returnHandler="ideal", incompleteGraphHandler="remove")
        firstgG = gGList[0]
        assert len(gGList) == 1 # may have to consider this
        assert is_connected(convert_to_nx(firstgG).to_undirected())
        assert firstgG.to_smiles() == MOLSMILES
        nodesSmarts = Counter([node.type for node in firstgG.nodes.values()])
        assert SITESDICT == nodesSmarts

        nodeEdges = Counter([(firstgG.nodes[edge1].type,firstgG.nodes[edge2].type) for edge1, _, edge2, _, bond_number in firstgG.edges])
        assert self.assert_equal_edgeDicts(EDGESDICT,nodeEdges), (EDGESDICT,nodeEdges)

    # def test_unifac_nodes_and_edges(self, unifac_groups):
    #     molecule = "c1c2ccccc2ccc1"
    #     gGList = fragment(molecule, unifac_groups, returnHandler="ideal", incompleteGraphHandler="keep") # is this what we want?
    #     firstgG = gGList[0]
    #     assert len(gGList) == 1 # may have to consider this
    #     assert is_connected(convert_to_nx(firstgG).to_undirected())
    #     assert firstgG.to_smiles() == molecule
    #     nodesSmarts = [node.pattern for node in firstgG.nodes]
    #     solution = ["C"]
    #     assert all(solution == nodesSmarts)
    #     nodes = ["ACH:8;AC:2"]
    #     edges = {"AC-AC":1, "ACH-AC":4, "ACH-ACH":8}
