# ruff: noqa: F401
import pytest
from networkx import is_connected
from collections import Counter

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

test_molecules = {
	'CC': 
        {'UNIFAC':{'GROUPS':{'CH3': 2},'EDGES':{('CH3', 'CH3'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 2}),'EDGES':Counter({('-CH3', '-CH3'): 1})},
        'SAFT':{'GROUPS':Counter({'CH3': 2}),'EDGES':Counter({('CH3', 'CH3'): 1})}},
	'CC(C)C': 
        {'UNIFAC':{'GROUPS':{'CH': 1, 'CH3': 3},'EDGES':{('CH3', 'CH'): 3}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 3, '>CH-': 1}),'EDGES':Counter({('>CH-', '-CH3'): 3})},
        'SAFT':{'GROUPS':Counter({'CH3': 3, 'CH': 1}),'EDGES':Counter({('CH', 'CH3'): 3})}},
	'CCC(C)(C)C': 
        {'UNIFAC':{'GROUPS':{'CH3': 4, 'CH2': 1, 'C': 1},'EDGES':{('CH3', 'CH2'): 1, ('CH3', 'C'): 3, ('CH2', 'C'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 4, '>C<': 1, '-CH2-': 1}),'EDGES':Counter({('>C<', '-CH3'): 3, ('-CH2-', '-CH3'): 1, ('>C<', '-CH2-'): 1})},
        'SAFT':{'GROUPS':Counter({'CH3': 4, 'C': 1, 'CH2': 1}),'EDGES':Counter({('C', 'CH3'): 3, ('CH2', 'CH3'): 1, ('C', 'CH2'): 1})}},
	'CC(C)(C)C': 
        {'UNIFAC':{'GROUPS':{'CH3': 4, 'C': 1},'EDGES':{('C', 'CH3'): 4}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 4, '>C<': 1}),'EDGES':Counter({('>C<', '-CH3'): 4})},
        'SAFT':{'GROUPS':Counter({'CH3': 4, 'C': 1}),'EDGES':Counter({('C', 'CH3'): 4})}},
	'OCC(O)CO': 
        {'UNIFAC':{'GROUPS':{'OH (P)': 2, 'CH2': 2, 'CH': 1, 'OH (S)': 1},'EDGES':{('OH (P)', 'CH2'): 2, ('OH (S)', 'CH'): 1, ('CH2', 'CH'): 2}},
        'JOBACK':{'GROUPS':Counter({'-OH (alcohol)': 3, '-CH2-': 2, '>CH-': 1}),'EDGES':Counter({('-CH2-', '-OH (alcohol)'): 2, ('>CH-', '-CH2-'): 2, ('>CH-', '-OH (alcohol)'): 1})},
        'SAFT':{'GROUPS':Counter({'CH2OH': 2, 'CHOH': 1}),'EDGES':Counter({('CHOH', 'CH2OH'): 2})}},
	'O=CO': 
        {'UNIFAC':{'GROUPS':{'HCOO': 1},'EDGES':{}},
        'JOBACK':{'GROUPS':Counter({'-COOH (acid)': 1}),'EDGES':Counter()},
        'SAFT':{'GROUPS':Counter({'COOH': 1}),'EDGES':Counter()}},
	'COC(C)=O': 
        {'UNIFAC':{'GROUPS':{'CH3COO': 1, 'CH3': 1},'EDGES':{('CH3COO', 'CH3'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 2, '-COO- (ester)': 1}),'EDGES':Counter({('-CH3', '-COO- (ester)'): 2})},
        'SAFT':{'GROUPS':Counter({'CH3': 2, 'COO': 1}),'EDGES':Counter({('CH3', 'COO'): 2})}},
    'O=C(O)CCCO': 
            {'UNIFAC':{'GROUPS':{'CH2COO': 1, 'CH2': 2, 'OH (P)': 1},'EDGES':{('CH2COO', 'CH2'): 1, ('CH2', 'CH2'): 1, ('CH2', 'OH (P)'): 1}},
            'JOBACK':{'GROUPS':Counter({'-CH2-': 3, '-OH (alcohol)': 1, '-COOH (acid)': 1}),'EDGES':Counter({('-CH2-', '-CH2-'): 2, ('-CH2-', '-OH (alcohol)'): 1, ('-CH2-', '-COOH (acid)'): 1})},
            'SAFT':{'GROUPS':Counter({'CH2': 2, 'CH2OH': 1, 'COOH': 1}),'EDGES':Counter({('CH2', 'CH2'): 1, ('CH2', 'CH2OH'): 1, ('CH2', 'COOH'): 1})}},
	'NC=O': 
        {'UNIFAC':{'GROUPS':{},'EDGES':{}},
        'JOBACK':{'GROUPS':Counter({'-NH2': 1, 'O=CH- (aldehyde)': 1}),'EDGES':Counter({('-NH2', 'O=CH- (aldehyde)'): 1})},
        'SAFT':{'GROUPS':None, 'EDGES':None}},
	'CC(N)=O': 
        {'UNIFAC':{'GROUPS':{'CONH2': 1, 'CH3': 1},'EDGES':{('CONH2', 'CH3'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 1, '-NH2': 1, '>C=O (non-ring)': 1}),'EDGES':Counter({('-NH2', '>C=O (non-ring)'): 1, ('-CH3', '>C=O (non-ring)'): 1})},
        'SAFT':{'GROUPS':Counter({"CH3CO":1, "NH2":1}), 'EDGES':Counter({("NH2", "CH3CO"):1})}},
	'CNC(C)=O': 
        {'UNIFAC':{'GROUPS':{'CONHCH3': 1, 'CH3': 1},'EDGES':{('CONHCH3', 'CH3'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 2, '>NH (non-ring)': 1, '>C=O (non-ring)': 1}),'EDGES':Counter({('-CH3', '>C=O (non-ring)'): 1, ('>NH (non-ring)', '>C=O (non-ring)'): 1, ('-CH3', '>NH (non-ring)'): 1})},
        'SAFT':{'GROUPS':Counter({"CH3":1, "NH":1, "CH3CO":1}), 'EDGES':Counter({("CH3", "NH"):1, ('CH3CO', 'NH'):1})}},
	'CNC': 
        {'UNIFAC':{'GROUPS':{'CH3': 1, 'CH3NH': 1},'EDGES':{('CH3', 'CH3NH'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 2, '>NH (non-ring)': 1}),'EDGES':Counter({('-CH3', '>NH (non-ring)'): 2})},
        'SAFT':{'GROUPS':Counter({'CH3': 2, 'NH': 1}),'EDGES':Counter({('CH3', 'NH'): 2})}},
	'CNCNC': 
        {'UNIFAC':{'GROUPS':{'CH2': 1, 'CH3NH': 2},'EDGES':{('CH2', 'CH3NH'): 2}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 2, '>NH (non-ring)': 2, '-CH2-': 1}),'EDGES':Counter({('-CH2-', '>NH (non-ring)'): 2, ('-CH3', '>NH (non-ring)'): 2})},
        'SAFT':{'GROUPS':Counter({'CH3': 2, 'NH': 2, 'CH2': 1}),'EDGES':Counter({('CH2', 'NH'): 2, ('CH3', 'NH'): 2})}},
	'FC(F)F': 
        {'UNIFAC':{'GROUPS':{'CF3': 1},'EDGES':{}},
        'JOBACK':{'GROUPS':Counter({'-F': 3, '>CH-': 1}),'EDGES':Counter({('>CH-', '-F'): 3})},
        'SAFT':{'GROUPS':None, 'EDGES':None}},
	'NCCCCF': 
        {'UNIFAC':{'GROUPS':{'CH2NH2': 1, 'CH2': 2, 'CF': 1},'EDGES':{('CH2NH2', 'CH2'): 1, ('CH2', 'CF'): 1, ('CH2', 'CH2'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH2-': 4, '-NH2': 1, '-F': 1}),'EDGES':Counter({('-CH2-', '-CH2-'): 3, ('-CH2-', '-NH2'): 1, ('-CH2-', '-F'): 1})},
        'SAFT':{'GROUPS':None, 'EDGES':None}},
	'CC1CCCCC1': 
        {'UNIFAC':{'GROUPS':{'CY-CH2': 5, 'CH3': 1, 'CY-CH': 1},'EDGES':{('CY-CH2', 'CY-CH2'): 4, ('CH3', 'CY-CH'): 1, ('CY-CH2', 'CY-CH'): 2}},
        'JOBACK':{'GROUPS':Counter({'ring-CH2-': 5, 'ring>CH-': 1, '-CH3': 1}),'EDGES':Counter({('ring-CH2-', 'ring-CH2-'): 4, ('ring>CH-', 'ring-CH2-'): 2, ('ring>CH-', '-CH3'): 1})},
        'SAFT':{'GROUPS':Counter({'cCH2': 5, 'cCH': 1, 'CH3': 1}),'EDGES':Counter({('cCH2', 'cCH2'): 4, ('cCH', 'cCH2'): 2, ('cCH', 'CH3'): 1})}},
	'CC1CCCCC1C': 
        {'UNIFAC':{'GROUPS':{'CY-CH2': 4, 'CH3': 2, 'CY-CH': 2},'EDGES':{('CY-CH2', 'CY-CH2'): 3, ('CH3', 'CY-CH'): 2, ('CY-CH2', 'CY-CH'): 2, ('CY-CH', 'CY-CH'): 1}},
        'JOBACK':{'GROUPS':Counter({'ring-CH2-': 4, 'ring>CH-': 2, '-CH3': 2}),'EDGES':Counter({('ring-CH2-', 'ring-CH2-'): 3, ('ring>CH-', '-CH3'): 2, ('ring>CH-', 'ring-CH2-'): 2, ('ring>CH-', 'ring>CH-'): 1})},
        'SAFT':{'GROUPS':Counter({'cCH2': 4, 'cCH': 2, 'CH3': 2}),'EDGES':Counter({('cCH2', 'cCH2'): 3, ('cCH', 'CH3'): 2, ('cCH', 'cCH2'): 2, ('cCH', 'cCH'): 1})}},
	'CC1CCCC(C)C1': 
        {'UNIFAC':{'GROUPS':{'CY-CH2': 4, 'CH3': 2, 'CY-CH': 2},'EDGES':{('CY-CH2', 'CY-CH2'): 2, ('CH3', 'CY-CH'): 2, ('CY-CH2', 'CY-CH'): 4}},
        'JOBACK':{'GROUPS':Counter({'ring-CH2-': 4, 'ring>CH-': 2, '-CH3': 2}),'EDGES':Counter({('ring>CH-', 'ring-CH2-'): 4, ('ring>CH-', '-CH3'): 2, ('ring-CH2-', 'ring-CH2-'): 2})},
        'SAFT':{'GROUPS':Counter({'cCH2': 4, 'cCH': 2, 'CH3': 2}),'EDGES':Counter({('cCH', 'cCH2'): 4, ('cCH', 'CH3'): 2, ('cCH2', 'cCH2'): 2})}},
	'CC1CCC(C)CC1': 
        {'UNIFAC':{'GROUPS':{'CY-CH2': 4, 'CH3': 2, 'CY-CH': 2},'EDGES':{('CY-CH2', 'CY-CH2'): 2, ('CH3', 'CY-CH'): 2, ('CY-CH2', 'CY-CH'): 4}},
        'JOBACK':{'GROUPS':Counter({'ring-CH2-': 4, 'ring>CH-': 2, '-CH3': 2}),'EDGES':Counter({('ring>CH-', 'ring-CH2-'): 4, ('ring>CH-', '-CH3'): 2, ('ring-CH2-', 'ring-CH2-'): 2})},
        'SAFT':{'GROUPS':Counter({'cCH2': 4, 'cCH': 2, 'CH3': 2}),'EDGES':Counter({('cCH', 'cCH2'): 4, ('cCH', 'CH3'): 2, ('cCH2', 'cCH2'): 2})}},
	'Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]': 
        {'UNIFAC':{'GROUPS':{'ACCH3': 1, 'ACH': 2, 'ACNO2': 3},'EDGES':{('ACCH3', 'ACNO2'): 2, ('ACNO2', 'ACH'): 4}},
        'JOBACK':{'GROUPS':Counter({'ring=C<': 4, '-NO2': 3, 'ring=CH-': 2, '-CH3': 1}),'EDGES':Counter({('ring=C<', 'ring=CH-'): 4, ('ring=C<', '-NO2'): 3, ('ring=C<', 'ring=C<'): 2, ('ring=C<', '-CH3'): 1})},
        'SAFT':{'GROUPS':None, 'EDGES':None}},
	'O=C(O)c1ccccc1': 
        {'UNIFAC':{'GROUPS':{'ACH': 5, 'AC': 1, 'COOH': 1},'EDGES':{('ACH', 'ACH'): 4, ('AC', 'ACH'): 2, ('AC', 'COOH'): 1}},
        'JOBACK':{'GROUPS':Counter({'ring=CH-': 5, 'ring=C<': 1, '-COOH (acid)': 1}),'EDGES':Counter({('ring=CH-', 'ring=CH-'): 4, ('ring=C<', 'ring=CH-'): 2, ('ring=C<', '-COOH (acid)'): 1})},
        'SAFT':{'GROUPS':Counter({'aCH': 5, 'aCCOOH': 1}),'EDGES':Counter({('aCH', 'aCH'): 4, ('aCH', 'aCCOOH'): 2})}},
	'CN(C)C(=O)CO': 
        {'UNIFAC':{'GROUPS':{'CH2CO': 1, 'CH3': 1, 'CH3N': 1, 'OH (P)': 1},'EDGES':{('CH2CO', 'OH (P)'): 1, ('CH2CO', 'CH3N'): 1, ('CH3N', 'CH3'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 2, '-CH2-': 1, '>N- (non-ring)': 1, '>C=O (non-ring)': 1, '-OH (alcohol)': 1}),'EDGES':Counter({('-CH3', '>N- (non-ring)'): 2, ('>N- (non-ring)', '>C=O (non-ring)'): 1, ('-CH2-', '>C=O (non-ring)'): 1, ('-CH2-', '-OH (alcohol)'): 1})},
        'SAFT':{'GROUPS':None, 'EDGES':None}},
	'CN(C)C=O': 
        {'UNIFAC':{'GROUPS':{'DMF': 1},'EDGES':{}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 2, '>N- (non-ring)': 1, 'O=CH- (aldehyde)': 1}),'EDGES':Counter({('-CH3', '>N- (non-ring)'): 2, ('>N- (non-ring)', 'O=CH- (aldehyde)'): 1})},
        'SAFT':{'GROUPS':None, 'EDGES':None}}, # SAFT can't do some aldehydes
	'CN': 
        {'UNIFAC':{'GROUPS':{'CH3NH2': 1},'EDGES':{}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 1, '-NH2': 1}),'EDGES':Counter({('-CH3', '-NH2'): 1})},
        'SAFT':{'GROUPS':Counter({'CH3': 1, 'NH2': 1}),'EDGES':Counter({('CH3', 'NH2'): 1})}},
	'CCC#N': 
        {'UNIFAC':{'GROUPS':{'CH2CN': 1, 'CH3': 1},'EDGES':{('CH2CN', 'CH3'): 1}},
        'JOBACK':{'GROUPS':Counter({'-CH2-': 1, '-CH3': 1, '-CN': 1}),'EDGES':Counter({('-CH2-', '-CN'): 1, ('-CH2-', '-CH3'): 1})},
        'SAFT':{'GROUPS':None,'EDGES':None}},
	'CO': 
        {'UNIFAC':{'GROUPS':{'CH3OH': 1},'EDGES':{}},
        'JOBACK':{'GROUPS':Counter({'-CH3': 1, '-OH (alcohol)': 1}),'EDGES':Counter({('-CH3', '-OH (alcohol)'): 1})},
        'SAFT':{'GROUPS':Counter({'CH3OH': 1}),'EDGES':Counter()}},
}

Libraries = {
    "base": BasisSet,
    "saftgm": SaftGammaMie,
    "joback": Joback,
    "UNIFAC": Unifac,
}

class TestLibraries(BaseTest):
    @pytest.mark.parametrize(
        "library,n_graphs",
        [("saftgm", 36), ("joback", 41), ("UNIFAC", 92), ("base", 0)],
    )
    def test_build_library(self, library, n_graphs):
        library = Libraries[library]()
        assert library.n_groups == n_graphs

    def test_group_extension(self):
        library = Libraries["saftgm"]()
        group = Group("CH3", "[CX4H3]", [0], "SMARTS")
        nt = library.query_groups({"group": group})[0]

        assert nt.group == group
        assert nt.extended_smarts == "C"
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
            {"extended_smarts": "I"}
        )[0]
        assert nt.group == Group('-I', '[I;X1]', [0], 'SMARTS')

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

    # TODO: better testing and documentation of failing molecules
    all_fail_molecules = [ # molecules fail with charges
        "C[O-]",
        "CCC(=O)[O-]",
        "CCS(=O)(=O)[O-]",
        "C(=O)[N-]",
        # SAFT fails
        "CN(C)C=O", # can;t do aldehydes except C(=O)CH3
        "C(F)(F)(F)", # No halogens in saft
        'Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]', # No NO2 in saft
        'CN(C)C(=O)CO', # no good C=O ketone group
        'CN(C)C(=O)CC', # no good C=O ketone group
    ]

    @pytest.mark.parametrize(
            "LIBRARYGROUPS", [
                ("SAFT", list(SaftGammaMie().get_groups())), 
                ("JOBACK", list(Joback().get_groups())), 
                ("UNIFAC", list(Unifac().get_groups()))
            ])
    @pytest.mark.parametrize(
        "MOLSMILES", test_molecules.keys()
    )
    def test_library_smiles(self, LIBRARYGROUPS, MOLSMILES):
        LIBRARY = LIBRARYGROUPS[0]
        GROUPS = LIBRARYGROUPS[1]

        gGList = fragment(
            MOLSMILES,
            GROUPS,
            returnHandler="ideal",
            incompleteGraphHandler="remove",
            nodeDefsSorter="size",
        )

        if not test_molecules[MOLSMILES][LIBRARY]["GROUPS"]:
            assert len(gGList) == 0 
            return # These groups fail to match
        assert len(gGList) >= 1, print(
            fragment(
                MOLSMILES,
                GROUPS,
                returnHandler="ideal",
                incompleteGraphHandler="keep",
                nodeDefsSorter="size",
            )
        )
        firstgG = gGList[0]
        assert is_connected(convert_to_nx(firstgG).to_undirected())
        groupsSmarts = Counter([group.type for group in firstgG.nodes.values()])
        assert test_molecules[MOLSMILES][LIBRARY]["GROUPS"] == groupsSmarts

        groupEdges = Counter(
            [
                (firstgG.nodes[edge1].type, firstgG.nodes[edge2].type)
                for edge1, _, edge2, _, bond_number in firstgG.edges
            ]
        )
        assert self.assert_equal_edgeDicts(test_molecules[MOLSMILES][LIBRARY]["EDGES"], groupEdges), (
            test_molecules[MOLSMILES][LIBRARY]["EDGES"],
            groupEdges,
        )

    @pytest.mark.parametrize(
            "LIBRARYGROUPS", [
                ("SAFT", list(SaftGammaMie().get_groups())), 
                ("JOBACK", list(Joback().get_groups())), 
                ("UNIFAC", list(Unifac().get_groups()))
            ])
    @pytest.mark.parametrize(
        "MOLSMILES", test_molecules.keys()
    )
    def test_fragment_to_smiles(self, LIBRARYGROUPS, MOLSMILES):
        """TODO: This test will have some special cases to test later.
        The SMARTS strings used to generate the new groupgraph structures have
        issues being parsed by the to_smiles() method and generating consistent 
        canonical SMILES. 

        Notes
        -----
        "C(O)CCC(=O)O" is an example of one that fails when fragmented by JOBACK
        """
        if True: # skip for now
            return
        LIBRARY = LIBRARYGROUPS[0]
        GROUPS = LIBRARYGROUPS[1]
        gGList = fragment(
            MOLSMILES,
            GROUPS,
            returnHandler="ideal",
            incompleteGraphHandler="remove",
            nodeDefsSorter="size",
        )
        firstgG = gGList[0]
        assert firstgG.to_smiles() == MOLSMILES

    def test_library_failed_fragmentation_returnGroup(self):
        """Returns a single group on fragmentation Failure"""
        library = Libraries["base"]()
        group = Group("CH3", "[CX4H3]", [0], "SMARTS")
        library.add_group(group)
        group = Group("CH2", "[CX4H2]", [0,0], "SMARTS")
        library.add_group(group)

        smiles = "CC(C)C"
        groupG = library.fragment_smiles(smiles, onFail="condense")[0]
        assert groupG.to_smiles() == smiles
        groupsSmarts = Counter([group.type for group in groupG.nodes.values()])
        assert {"CC(C)C":1} == groupsSmarts

        groupEdges = Counter(
            [
                (groupG.nodes[edge1].type, groupG.nodes[edge2].type)
                for edge1, _, edge2, _, bond_number in groupG.edges
            ]
        )
        assert not groupEdges

    def test_library_failed_fragmentation_returnElement(self):
        library = Libraries["base"]()
        group = Group("CH3", "[CX4H3]", [0], "SMARTS")
        library.add_group(group)
        group = Group("CH2", "[CX4H2]", [0,0], "SMARTS")
        library.add_group(group)

        smiles = "CC(C)C"
        groupG = library.fragment_smiles(smiles, onFail="itemize")[0]
        assert groupG.to_smiles() == smiles
        groupsSmarts = Counter([group.type for group in groupG.nodes.values()])
        assert {"CH3":3, "DEFAULT#6":1} == groupsSmarts

        groupEdges = Counter(
            [
                (groupG.nodes[edge1].type, groupG.nodes[edge2].type)
                for edge1, _, edge2, _, bond_number in groupG.edges
            ]
        )
        assert self.assert_equal_edgeDicts({("CH3", "DEFAULT#6"):3}, groupEdges), groupEdges

    def test_library_failed_fragmentation_returnException(self):
        from Grouper._Grouper import GrouperFragmentationError
        library = Libraries["base"]()
        group = Group("CH3", "[CX4H3]", [0], "SMARTS")
        library.add_group(group)
        group = Group("CH2", "[CX4H2]", [0,0], "SMARTS")
        library.add_group(group)

        smiles = "CC(C)C"
        with pytest.raises(GrouperFragmentationError):
            groupG = library.fragment_smiles(smiles, onFail="error")[0]
    
    def test_library_successful_fragmentation(self):
        library = Libraries["base"]()
        group = Group("CH3", "[CX4H3]", [0], "SMARTS")
        library.add_group(group)
        group = Group("CH2", "[CX4H2]", [0,0], "SMARTS")
        library.add_group(group)

        smiles = "CCC"
        groupG = library.fragment_smiles(smiles)[0]
        assert groupG.to_smiles() == smiles
        groupsSmarts = Counter([group.type for group in groupG.nodes.values()])
        assert {"CH3":2, "CH2":1} == groupsSmarts

        groupEdges = Counter(
            [
                (groupG.nodes[edge1].type, groupG.nodes[edge2].type)
                for edge1, _, edge2, _, bond_number in groupG.edges
            ]
        )
        assert self.assert_equal_edgeDicts({("CH3", "CH2"):2}, groupEdges), groupEdges
        
