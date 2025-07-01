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
from Grouper.tests.base_test import BaseTest
from Grouper.utils import convert_to_nx

Libraries = {
    "base": BasisSet,
    "saftgm": SaftGammaMie,
    "joback": Joback,
    "UNIFAC": Unifac,
}


class TestLibraries(BaseTest):
    @pytest.mark.parametrize(
        "library,n_graphs",
        [("saftgm", 1), ("joback", 41), ("UNIFAC", 92), ("base", 0)],
    )
    def test_build_library(self, library, n_graphs):
        library = Libraries[library]()
        assert library.n_groups == n_graphs

    def test_group_extension(self):
        library = Libraries["saftgm"]()
        group = Group("-CH3", "[CH3]", [0], "SMARTS")
        nt = library.query_groups({"group": group})[0]

        assert nt.group == group
        assert nt.extended_smarts == "[CX4H3]"
        assert nt.doi == "https://doi.org/10.1080/00268978800101601"
        assert nt.priority is None

        nt = GroupExtension(group, "", "[CX4H3]", 1)
        assert nt.priority == 1

    def test_add_group(self):
        library = BasisSet()
        library.add_group(Group("-CH3", "[CH3]", [0], "SMARTS"), "", "[CX4H3]", None)
        assert library.n_groups == 1

    def test_query_group(self):
        library = Libraries["joback"]()
        nt = library.query_groups(
            {"extended_smarts": "[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]"}
        )[0]
        assert nt.group == Group("=C<", "C", [0, 0, 0])

    def test_list_groups(self):
        library = Libraries["joback"]()
        nodes = list(library.get_groups())
        assert len(nodes) == 41 == library.n_groups

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
        if edge2:  # if any missed in edge2
            print(f"Generated extra edge {edge2=} in fragmented graph.")
            return False
        return True

    # TODO: How to handle "CN", which has a UNIFAC Group GCPair(raw"[CX4;H3][NX3;H2]","CH3NH2"), but no ports
    # TODO: More failing molecules that I tried to test -NCC
    unifac_fail_molecules = [
        "C(=O)N",
        "OCC(=O)N(C)C",
        "C(=O)N(C)C",
        "CN",
        "C[O-]",
        "CCC(=O)[O-]",
        "CCS(=O)(=O)[O-]",
        "CCC#N",
        "C(=O)[N-]",
        "CO",
    ]
    test_molecules = [  # NOTE that "C(=O)N", or formaldehyde doesn't work
        "CC",
        "CC(C)C",
        "CCC(C)(C)C",
        "CC(C)(C)C",  # alkanes
        "OCC(O)CO",  # alcohols
        "O=CO",
        "COC(C)=O",
        "OCCCC(=O)O",
        "CC(N)=O",
        "CNC(C)=O",  # carboxyls/amides
        "CNC",
        "CNCNC",  # amines
        "FC(F)F",
        "NCCCCF",  # fluorines
        "CC1CCCCC1",
        "CC1CCCCC1C",
        "CC1CCCC(C)C1",
        "CC1CCC(C)CC1",  # rings
        "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]",
        "O=C(O)c1ccccc1",  # phenyls
    ]
    groupsList = [
        {"CH3": 2},
        {"CH": 1, "CH3": 3},
        {"CH3": 4, "CH2": 1, "C": 1},
        {"CH3": 4, "C": 1},  # alkanes
        {"OH (P)": 2, "CH2": 2, "CH": 1, "OH (S)": 1},  # alcohols
        {"HCOO": 1},
        {"CH3COO": 1, "CH3": 1},
        {"CH2COO": 1, "CH2": 2, "OH (P)": 1},
        {"CONH2": 1, "CH3": 1},
        {"CONHCH3": 1, "CH3": 1},  # carboxyl/amides
        {"CH3": 1, "CH3NH": 1},
        {"CH2": 1, "CH3NH": 2},  # amines
        {"CF3": 1},
        {"CH2NH2": 1, "CH2": 2, "CF": 1},  # fluorines
        {"CY-CH2": 5, "CH3": 1, "CY-CH": 1},
        {"CY-CH2": 4, "CH3": 2, "CY-CH": 2},  # rings
        {"CY-CH2": 4, "CH3": 2, "CY-CH": 2},
        {"CY-CH2": 4, "CH3": 2, "CY-CH": 2},  # rings
        {"ACCH3": 1, "ACH": 2, "ACNO2": 3},
        {"ACH": 5, "AC": 1, "COOH": 1},  # phenyl groups
    ]
    edgesList = [
        {("CH3", "CH3"): 1},
        {("CH3", "CH"): 3},
        {("CH3", "CH2"): 1, ("CH3", "C"): 3, ("CH2", "C"): 1},
        {("C", "CH3"): 4},  # alkanes
        {("OH (P)", "CH2"): 2, ("OH (S)", "CH"): 1, ("CH2", "CH"): 2},  # alcohols
        {},
        {("CH3COO", "CH3"): 1},
        {("CH2COO", "CH2"): 1, ("CH2", "CH2"): 1, ("CH2", "OH (P)"): 1},
        {("CONH2", "CH3"): 1},
        {("CONHCH3", "CH3"): 1},  # carbxoyls/amides
        {("CH3", "CH3NH"): 1},
        {("CH2", "CH3NH"): 2},  # amines
        {},
        {("CH2NH2", "CH2"): 1, ("CH2", "CF"): 1, ("CH2", "CH2"): 1},  # fluorines
        {("CY-CH2", "CY-CH2"): 4, ("CH3", "CY-CH"): 1, ("CY-CH2", "CY-CH"): 2},
        {
            ("CY-CH2", "CY-CH2"): 3,
            ("CH3", "CY-CH"): 2,
            ("CY-CH2", "CY-CH"): 2,
            ("CY-CH", "CY-CH"): 1,
        },  # rings
        {("CY-CH2", "CY-CH2"): 2, ("CH3", "CY-CH"): 2, ("CY-CH2", "CY-CH"): 4},
        {("CY-CH2", "CY-CH2"): 2, ("CH3", "CY-CH"): 2, ("CY-CH2", "CY-CH"): 4},  # rings
        {("ACCH3", "ACNO2"): 2, ("ACNO2", "ACH"): 4},
        {("ACH", "ACH"): 4, ("AC", "ACH"): 2, ("AC", "COOH"): 1},  # phenyl groups
    ]

    @pytest.fixture(scope="session")
    def unifac_groups(self):
        unifac = Unifac()
        return list(unifac.get_groups())

    @pytest.mark.parametrize(
        "MOLSMILES,SITESDICT,EDGESDICT", zip(test_molecules, groupsList, edgesList)
    )
    def test_unifac_smiles(self, MOLSMILES, SITESDICT, EDGESDICT, unifac_groups):
        from collections import Counter

        # molecule = test_molecules[0]
        gGList = fragment(
            MOLSMILES,
            unifac_groups,
            returnHandler="ideal",
            incompleteGraphHandler="remove",
            nodeDefsSorter="size",
        )
        # may have to consider this just to be 1
        assert len(gGList) >= 1, print(
            fragment(
                MOLSMILES,
                unifac_groups,
                returnHandler="ideal",
                incompleteGraphHandler="keep",
                nodeDefsSorter="size",
            )
        )
        firstgG = gGList[0]
        assert is_connected(convert_to_nx(firstgG).to_undirected())
        assert firstgG.to_smiles() == MOLSMILES
        groupsSmarts = Counter([group.type for group in firstgG.nodes.values()])
        assert SITESDICT == groupsSmarts

        groupEdges = Counter(
            [
                (firstgG.nodes[edge1].type, firstgG.nodes[edge2].type)
                for edge1, _, edge2, _, bond_number in firstgG.edges
            ]
        )
        assert self.assert_equal_edgeDicts(EDGESDICT, groupEdges), (
            EDGESDICT,
            groupEdges,
        )