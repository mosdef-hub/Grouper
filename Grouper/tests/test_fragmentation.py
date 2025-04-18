"""Tests on fragmentation of a molecule into GroupGraph."""

from Grouper import Group, GroupGraph
from Grouper.fragmentation import (
    _get_first_compatible_tuples,
    _get_hubs_from_string,
    _get_maximal_compatible_tuples,
    _get_next_available_port_from_hub,
    _smarts_with_ports,
    fragment,
)
from Grouper.tests.base_test import BaseTest


class TestGeneralFragmentations(BaseTest):
    def test_overlap_fragment(self):
        smiles = "C(C)CO"

        query1 = "[C]C"
        query2 = "[O]"
        query3 = "[C]"
        queries = [query1, query2, query3]
        solution = 2

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == 2, question

    def test_two_fragmentations(self):
        smiles = "CCOC"
        query1 = "[O]C"
        query2 = "[C]"
        queries = [query1, query2]
        solution = 2

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == 2, question

    def test_larger_fragmentation(self):
        smiles = "C(C)CO"
        query1 = "[C]C"
        query2 = "[O]C"
        query3 = "[C]"
        query4 = "[O]"
        queries = [query1, query2, query3, query4]
        solution = 2

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == solution, question

    def test_branching_fragmentations(self):
        smiles = "C(O)CO"  # test a branch
        query1 = "[O]C"
        query2 = "[C]"
        query3 = "[O]"
        queries = [query1, query2, query3]
        solution = 1

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == solution, question

    def test_optimal_fragmentations_only(self):
        smiles = "COCO"  # middle OC should not be matched
        query1 = "[O]C"
        query2 = "[C]"
        query3 = "[O]"
        queries = [query1, query2, query3]
        solution = 1

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == solution, question

    def test_trim_potential_fragmentations(self):
        smiles = "COCO"  # test when could find wrong initial structure
        query1 = "[O]C"
        queries = [query1]
        solution = 1

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == solution, question

    def test_complex_possible_OC_fragmentations_from_strings(self):
        smiles = "C(O)(O)(O)OC"
        query1 = "[O]C"
        query2 = "[C]C"
        query3 = "C"
        query4 = "O"
        queries = [query1, query2, query3, query4]
        solution = 3

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == solution, question

    def test_fragment_ring(self):
        smiles = "CCC1CCCCC1"
        query1 = "[C;R]"
        query2 = "[C]C"
        queries = [query1, query2]
        solution = 1

        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == solution, question
        assert (
            len(question[0].nodes) == 7
        ), question  # 7 groups, 1 ethyl and 6 ring carbons

    def test_on_larger_molecules(self):
        queries = ["[C;R]", "[C](=O)O", "[C]O", "[C]C", "[O;X1]", "[C]", "[N]"]
        smiles = "CNC(C(C)O)C(N)CCN"
        # smiles = 'C-N-C(C-(C)-O)C-(N)-CC-N'
        solution = 1
        question = fragment(smiles, queries, returnHandler="ideal")
        assert len(question) == solution, question
        assert len(question[0].nodes) == 8, question  # 1 CO, 2 CC, 2 C, 3 N
        assert len(question[0].edges) == 7, (
            len(question[0].edges),
            question,
        )  # 8 nodes - 1 - 0 rings

    def test_incomplete_fragmentation(self):
        pass

        # Test 11: Duplicate groups

    def test_duplicate_groups_in_nodedefs(self):
        pass

    def test_bond_specification_in_nodeDefs(self):
        pass

    def test_full_fragmentation(self):
        node_defs = set()
        node_defs.add(Group("hydroxyl", "O", [0]))
        node_defs.add(Group("alkene", "C=C", [0, 0, 1, 1]))

        solution = GroupGraph()
        solution.add_node("hydroxyl", "O", [0])
        solution.add_node("alkene", "C=C", [0, 0, 1, 1])
        solution.add_edge((0, 0), (1, 0))

        smiles = "C=CO"

        question = fragment(smiles, node_defs)
        assert solution in question

    def test_full_fragmentation_double_bond(self):
        node_defs = set()
        node_defs.add(Group("amine", "N", [0, 0]))
        node_defs.add(Group("alkene", "C=C", [0, 0, 1, 1]))

        solution = GroupGraph()
        solution.add_node("amine", "N", [0, 0])
        solution.add_node("alkene", "C=C", [0, 0, 1, 1])
        solution.add_edge((0, 0), (1, 0))

        smiles = "C=CN"

        question = fragment(smiles, node_defs)
        assert solution in question

    def test_invalid_smiles(self):
        node_defs = set()
        smiles = ""

        question = fragment(smiles, node_defs)
        assert question == []

        node_defs = set()
        smiles = "C==C"  # Invalid SMILES

        try:
            fragment(smiles, node_defs)
        except Exception as e:
            assert isinstance(
                e, ValueError
            )  # or the appropriate exception class for invalid SMILES

    def test_specific_graph_alkene(self):
        node_defs = set()
        node_defs.add(Group("amine", "N", [0, 0]))
        node_defs.add(Group("alkene", "C=C", [0, 0, 1, 1]))
        node_defs.add(Group("oxygen", "O", [0, 0]))

        solution = GroupGraph()
        solution.add_node("oxygen", "O", [0, 0])
        solution.add_node("alkene", "C=C", [0, 0, 1, 1])
        solution.add_node("amine", "N", [0])
        solution.add_edge((0, 0), (1, 0))
        solution.add_edge((0, 1), (2, 0))

        smiles = "C=CON"

        question = fragment(smiles, node_defs)
        assert solution in question

    def test_multi_fragment_COO(self):
        node_defs = set()
        node_defs.add(Group("oxyl", "O", [0, 0]))  # oxyl group
        node_defs.add(
            Group("ester", "[C](=O)[O]", [0, 2], pattern_type="SMARTS")
        )  # Ester group
        node_defs.add(Group("amine", "N", [0, 0, 0]))  # Amine group
        node_defs.add(Group("alkene_secondary_amine", "CNC", [0, 0, 1, 2, 2]))
        node_defs.add(Group("alkene", "C", [0, 0, 0]))

        solution = GroupGraph()
        solution.add_node("ester", "C(=O)O", [0, 2])
        solution.add_node("oxyl", "O", [0, 0])
        solution.add_node("ester", "C(=O)(O)", [0, 2])
        solution.add_edge((0, 1), (1, 0))
        solution.add_edge((1, 1), (2, 0))

        question = fragment("O=COOC(=O)O", node_defs)
        assert solution in question

        solution = GroupGraph()
        solution.add_node("ester", "C(=O)O", [0, 2])
        solution.add_node("alkene_secondary_amine", "C(N(C))", [0, 0])
        solution.add_node("amine", "N", [0, 0, 0])
        solution.add_edge((0, 1), (2, 0))
        solution.add_edge((1, 1), (2, 1))

        question = fragment("CNCNOC=O", node_defs)
        assert solution in question

        solution = GroupGraph()
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("amine", "N", [0, 0, 0])
        solution.add_node("oxyl", "O", [0, 0])
        solution.add_edge((0, 1), (1, 0))
        solution.add_edge((1, 1), (2, 0))

        question = fragment("CNO", node_defs)
        assert solution in question

        solution = GroupGraph()
        solution.add_node("oxyl", "O", [0, 0])
        solution.add_node("oxyl", "O", [0, 0])
        solution.add_node("ester", "C(=O)O", [0, 2])
        solution.add_edge((0, 1), (1, 0))
        solution.add_edge((1, 1), (2, 0))

        question = fragment("O=C(O)OO", node_defs)
        assert solution in question

    def test_fragment_match_hubs(self):
        node_defs = set()
        node_defs.add(Group("oxyl", "O", [0, 0], "SMARTS"))  # oxyl group
        node_defs.add(Group("ester", "C(=O)O", [0, 2], "SMARTS"))  # Ester group
        node_defs.add(Group("amine", "N", [0, 0, 0], "SMARTS"))  # Amine group
        node_defs.add(
            Group("alkene_secondary_amine", "C(N(C))", [0, 0], "SMARTS")
        )  # can be made of amine and alkene
        node_defs.add(Group("alkene", "C", [0, 0, 0], "SMARTS"))

        solution = GroupGraph()
        solution.add_node("alkene_secondary_amine", "C(N(C))", [0, 0], "SMARTS")
        solution.add_node("alkene", "C", [0, 0, 0], "SMARTS")
        solution.add_node("amine", "N", [0, 0, 0], "SMARTS")
        solution.add_node("alkene", "C", [0, 0, 0], "SMARTS")
        solution.add_node("alkene", "C", [0, 0, 0], "SMARTS")
        solution.add_node("alkene", "C", [0, 0, 0], "SMARTS")
        solution.add_node("alkene", "C", [0, 0, 0], "SMARTS")
        solution.add_edge((0, 1), (1, 0))
        solution.add_edge((1, 1), (2, 0))
        solution.add_edge((2, 1), (3, 0))
        solution.add_edge((3, 1), (4, 0))
        solution.add_edge((3, 2), (5, 0))
        solution.add_edge((5, 1), (6, 0))

        question = fragment("CCC(C)NCCNC", node_defs, matchHubs=True)
        assert question[0] == solution
        # without matchHubs, no matches are found
        question = fragment("CCC(C)NCCNC", node_defs, matchHubs=False)
        # assert question == [] # TODO: This test fails in linux only

        node1 = Group("1methyl", "C", [0, 0, 0])
        node2 = Group("2methyl", "C", [0, 0])

        solution = GroupGraph()
        solution.add_node("1methyl", "C", [0])
        solution.add_node("1methyl", "C", [0])
        solution.add_node("2methyl", "C", [0, 0])
        solution.add_edge((0, 0), (2, 0))
        solution.add_edge((1, 0), (2, 1))

        node_defs = {node2, node1}  # {node2}
        question = fragment("CCC", node_defs)
        assert solution in question, question

        solution = GroupGraph()
        solution.add_node("carbon", "C", [0, 0, 0, 0])
        solution.add_node("ether", "CO", [0, 1])
        solution.add_node("ether", "CO", [0, 1])
        solution.add_edge((0, 0), (2, 0))
        solution.add_edge((2, 1), (1, 0))

        node1 = Group("4methyl", "C", [0, 0, 0, 0])
        node2 = Group("methanol", "CO", [0, 0, 0, 1])  # methanol
        node3 = Group("ether", "CO", [0, 0, 0, 1])  # ether

        node_defs = {node3, node2, node1}  # {node2}
        question = fragment("CCOCO", node_defs)
        assert solution in question, question


class TestFragmentationUtilities(BaseTest):
    def test_compatible_query_substructures(self):
        substructsList = [(0, 1), (1, 2), (2, 3)]
        solution = [((0, 1), (2, 3))]
        test = _get_maximal_compatible_tuples(substructsList)
        assert test == solution, test

        substructsList = [(0, 1), (0, 2), (0, 3), (4, 5)]
        solution = [((0, 3), (4, 5)), ((0, 2), (4, 5)), ((0, 1), (4, 5))]
        test = _get_maximal_compatible_tuples(substructsList)
        assert test == solution, test

        substructsList = ((0, 1), (0, 2), (0, 3), (4, 5), (5, 6))
        solution = [
            ((0, 3), (5, 6)),
            ((0, 2), (5, 6)),
            ((0, 1), (5, 6)),
            ((0, 3), (4, 5)),
            ((0, 2), (4, 5)),
            ((0, 1), (4, 5)),
        ]
        test = _get_maximal_compatible_tuples(substructsList)
        assert test == solution, test

        substructsList = ((0, 1), (2, 3), (4, 5), (1, 2), (3, 4), (5, 6))
        solution = [
            ((0, 1), (3, 4), (5, 6)),
            ((0, 1), (2, 3), (5, 6)),
            ((1, 2), (3, 4), (5, 6)),
            ((0, 1), (2, 3), (4, 5)),
        ]
        test = _get_maximal_compatible_tuples(substructsList)
        assert test == solution, test

    def test_first_query_substructures(self):
        substructsList = ((0, 1), (0, 2), (0, 3), (4, 5), (5, 6))
        test = _get_first_compatible_tuples(substructsList)
        solution = [((0, 1), (4, 5))]
        assert test == solution, test

    def test_fragmentation_edges(self):
        from Grouper import Group

        smiles = "C(C)C(O)OCOC"
        nodesList = [
            Group("propyl", "[C]C", [0, 0, 0, 1, 1, 1], pattern_type="SMARTS"),
            Group("alc", "[C]O", [0, 0, 0, 1], pattern_type="SMARTS"),
            Group("alkyl", "[C]", [0, 0, 0, 0], pattern_type="SMARTS"),
            Group("ether", "O", [0, 0], pattern_type="SMARTS"),
        ]

        question = fragment(smiles, nodesList, returnHandler="ideal")[0]
        solution = {
            (1, 0, 0, 0, 1),
            (3, 3, 2, 0, 1),
            (2, 3, 1, 1, 1),
        }  # edge form 0 to 1, 2 to 3, 1 to 2
        assert question.edges == solution, question.edges

    def test_hubs_from_smiles_or_smarts_string(self):
        solution = [0, 0, 0, 0]
        test = _get_hubs_from_string("C")
        assert test == solution, test

        solution = [0, 0, 0, 1]
        test = _get_hubs_from_string("CO")
        assert test == solution, test

        solution = [0, 0, 0, 2, 3, 3, 4, 4, 6, 6, 9, 9, 9, 10, 11]
        test = _get_hubs_from_string("CHONPBrBClISiSeSFLi")
        assert test == solution, test

    def test_get_next_available_port_for_groupG(self):
        groupG = GroupGraph()
        for _ in range(3):
            groupG.add_node("a", "N", [0, 0, 0])
        groupG.add_edge((0, 0), (1, 0))
        groupG.add_edge((0, 1), (2, 0))
        solution = 2
        question = _get_next_available_port_from_hub(groupG, 0, 0)
        assert question == solution, question

    def test_smarts_with_ports(self):
        """TODO: Benzene aromatic rings, two character elements, +- charges, and &/, booleans."""
        # group = Group("amine", "CNC", [0,0], pattern_type="SMARTS")
        # solution = "[CH1,CH2,CH3][NH1][CH3]"
        # question = _smarts_with_ports(group.pattern, group.hubs)
        # assert question == solution, question

        # smarts = "[C]"
        # hubs = [0,0]
        # solution = "[CH2,CH3,CH4]"
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "[C]C"
        # hubs = [0,0]
        # solution = "[CH1,CH2,CH3][CH3]"
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "[C]ON"
        # hubs = [0,0,0,2,2]
        # solution = "[CH0,CH1,CH2,CH3][OH0][NH0,NH1,NH2]"
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "[C]CN"
        # hubs= [0,0]
        # solution = '[CH1,CH2,CH3][CH2][NH2]'
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "[C]CCC"
        # hubs = [0,0,1,1,2]
        # solution = '[CH1,CH2,CH3][CH0,CH1,CH2][CH1,CH2][CH3]'
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "[C;R](C)N" # rings and branches
        # hubs= [0,0]
        # solution = '[R;CH0,CH1,CH2]([CH3])[NH2]'
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "[C](C(=O)O)NC" # double bond
        # hubs= [0,0,4,5,5,5]
        # solution = "[CH0,CH1,CH2]([CH0]([OH0])[OH1])[NH0,NH1][CH0,CH1,CH2,CH3]"
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "O=C1CCC(C)CC1" # ring
        # hubs= [5, 7]
        # solution = "[OH0][CH0][CH2][CH2][CH1]([CH2,CH3])[CH2][CH1,CH2]"
        # question = _smarts_with_ports(smarts, hubs)
        # assert question == solution, question

        # smarts = "c1ccccc1" # benzene
        # hubs= [0,1,2]
        # solution = "[c1][c][c][c][c][c1]"
        # question = _smarts_with_ports(smarts, hubs)
        # # No support for aromatics yet
        # with pytest.raises(AssertionError):
        #     assert question == solution, question

        # smarts = "[C&R]" # benzene
        # hubs= [0,0]
        # solution = "[R;CH2,CH3,CH4]"
        # question = _smarts_with_ports(smarts, hubs)
        # # No support boolean &, yet
        # with pytest.raises(AssertionError):
        #     assert question == solution, question

        # smarts = "[C]([Br])[Cl]" # two character elements
        # hubs= [0]
        # solution = "[CH1,CH2]([BrH0])[Cl;H0]"
        # question = _smarts_with_ports(smarts, hubs)
        # # No support two character elements yet
        # with pytest.raises(AssertionError):
        #     assert question == solution, question

        # smarts = "[C][N+](=O)[O-]" # charges
        # hubs= [0]
        # solution = "[CH2,CH3][N+H0](OH0)[O-H0]"
        # question = _smarts_with_ports(smarts, hubs)
        # # No support for charges
        # with pytest.raises(AssertionError):
        #     assert question == solution, question

        smarts = "C(N(C))"  # charges
        hubs = [0, 0]
        solution = "[CH1,CH2,CH3]([NH1]([CH3]))"
        question = _smarts_with_ports(smarts, hubs)
        # No support for charges
        assert question == solution, question


class TestFragmentationOptions(BaseTest):
    ###########################
    # Test fragment rules
    ###########################
    def test_fragmentation_incompleteGraphHandler(self):
        # Test 1.1: incompleteGraphHandler remove
        node_defs = {Group("amine", "N", [0, 0]), Group("alkene", "C", [0, 0, 0, 0])}
        smiles = "COCN"

        question = fragment(smiles, node_defs, incompleteGraphHandler="remove")
        assert not question, question

        # Test 1.2: incompleteGraphHandler remove unconnected graph
        node_defs = {Group("amine", "N", [0]), Group("alkene", "C", [0, 0, 0, 0])}
        smiles = "CNC"

        question = fragment(smiles, node_defs, incompleteGraphHandler="remove")
        assert not question, question

        # Test 1.3: incompleteGraphHandler keep
        node_defs = {Group("alkene", "C", [0, 0, 0, 0])}
        smiles = "CNCO"

        question = fragment(smiles, node_defs, incompleteGraphHandler="keep")[0]
        solution = GroupGraph()
        solution.add_node("alkene", "C", [0, 0, 0, 0])  # only match the carbon atom
        solution.add_node("alkene", "C", [0, 0, 0, 0])  # only match the carbon atom
        assert question == solution, question

        # Test 1.4: incompleteGraphHandler raise error
        node_defs = {Group("amine", "N", [0, 0]), Group("alkene", "C", [0, 0, 0, 0])}
        smiles = "CNO"
        try:
            fragment(smiles, node_defs, incompleteGraphHandler="raise error")
        except ValueError as e:
            assert (
                str(e)
                == "No complete graphs found. Please use a different set of nodDefs."
            ), e

    def test_fragmentation_nodeDefsSorter(self):
        # Test 2.1: nodeDefsSorter "size", "priority" "list"
        #   Will return different groups ordered
        node_defs = [
            Group("amine", "N", [0, 0]),
            Group("alkyl", "C", [0, 0]),
            Group("propyl", "CCC", [0, 2]),
            Group("ether", "COC", [0, 2]),
        ]  # nodeDefsSorter "size" on groups with same number of heavy atoms, then total mass
        smiles = "COCCCCNC"

        solution = [  # oxygen is heavier than carbon, so is nitrogen
            Group("ether", "COC", [0, 2]),
            Group("propyl", "CCC", [0, 2]),
            Group("amine", "N", [0, 0]),
            Group("alkyl", "C", [0, 0]),
        ]
        question = fragment(smiles, node_defs, nodeDefsSorter="size")[0]
        for i in range(len(question.nodes)):
            assert question.nodes[i] == solution[i]  # , (question, i)

        # Test 2.2 nodeDefsSorter "priority"
        # TODO: test for priority once nodeExtensions are implemented

        # Test 2.3  nodeDefsSorter "list"
        node_defs = [
            Group("amine", "N", [0, 0]),
            Group("ether", "COC", [0, 2]),
            Group("alkyl", "C", [0, 0]),
            Group("propyl", "CCC", [0, 2]),
        ]
        smiles = "COCCCCNC"

        solution = [  # keep original order
            Group("amine", "N", [0, 0]),
            Group("ether", "COC", [0, 2]),
            Group("alkyl", "C", [0, 0]),
            Group("alkyl", "C", [0, 0]),
            Group("alkyl", "C", [0, 0]),
            Group("alkyl", "C", [0, 0]),
        ]
        question = fragment(smiles, node_defs, nodeDefsSorter="list")[0]
        for i in range(len(question.nodes)):
            assert question.nodes[i] == solution[i], (question, i)

    def test_fragmentation_returnHandler(self):
        # Test 3.1 returnHandler: "ideal"
        node_defs = [
            Group("ester", "COC", [0, 0, 2, 2]),
            Group("amine", "CN", [0, 0, 1, 1]),
            Group("alkyne", "[C]#[C]", [0, 0, 1, 1], pattern_type="SMARTS"),
            Group("alkane", "C", [0, 0, 0, 0]),
        ]
        smiles = "C#CC(OC)CNCNC"

        solution = GroupGraph()
        solution.add_node("ester", "COC", [0, 0, 2, 2])
        solution.add_node("amine", "CN", [0, 0, 1, 1])
        solution.add_node("amine", "CN", [0, 0, 1, 1])
        solution.add_node("alkyne", "C#C", [0, 0, 1, 1])
        solution.add_node("alkane", "C", [0, 0, 0, 0])
        solution.add_edge((0, 0), (1, 0))
        solution.add_edge((1, 2), (2, 0))
        solution.add_edge((0, 1), (3, 0))
        solution.add_edge((2, 2), (4, 0))
        question = fragment(smiles, node_defs, returnHandler="ideal")[0]
        assert question == solution, question

        # Test 3.2 returnHandler "ideal"
        solution1 = GroupGraph()
        solution1.add_node("ester", "COC", [0, 0, 1, 1])
        solution1.add_node("amine", "CN", [0, 0, 1, 1])
        solution1.add_node("amine", "CN", [0, 0, 1, 1])
        solution1.add_node("alkyne", "C#C", [0, 0, 1, 1])
        solution1.add_node("alkane", "C", [0, 0, 0, 0])
        solution1.add_edge((0, 0), (1, 0))
        solution1.add_edge((1, 2), (2, 0))
        solution1.add_edge((0, 1), (3, 0))
        solution1.add_edge((2, 2), (4, 0))
        question = fragment(smiles, node_defs, returnHandler="ideal")
        assert question[0] == solution  # both are equidealy good interpretations
        assert question[1] == solution1

        # Test 3.3 returnHandler "quick"
        node_defs = [
            Group("ethyl", "CC", [0, 0, 1, 1]),
            Group("amine", "N", [0, 0]),
            Group("alkene", "C", [0, 0, 0, 0]),
        ]
        smiles = "CCNCN"

        question = fragment(smiles, node_defs, returnHandler="quick")[0]
        solution = GroupGraph()
        solution.add_node("ethyl", "CC", [0, 0, 1, 1])
        solution.add_node("amine", "N", [0, 0])
        solution.add_node("amine", "N", [0, 0])
        solution.add_node("alkene", "C", [0, 0, 0, 0])
        solution.add_edge((0, 2), (1, 0))
        solution.add_edge((1, 0), (2, 0))
        solution.add_edge((2, 0), (3, 0))
        assert question == solution, question

        # Test 3.4 returnHandler "exhaustive"
        node_defs = set()
        node_defs.add(Group("oxyl", "O", [0, 0]))  # oxyl group
        node_defs.add(Group("ester", "C(=O)O", [0, 2]))  # Ester group
        node_defs.add(Group("amine", "N", [0, 0, 0]))  # Amine group
        node_defs.add(
            Group("alkene_secondary_amine", "C(N(C))", [0, 0])
        )  # can be made of amine and alkene
        node_defs.add(Group("alkene", "C", [0, 0, 0]))

        solution = GroupGraph()
        solution.add_node("alkene_secondary_amine", "C(N(C))", [0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("amine", "N", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_edge((0, 1), (1, 0))
        solution.add_edge((1, 1), (2, 0))
        solution.add_edge((2, 1), (3, 0))
        solution.add_edge((3, 1), (4, 0))
        solution.add_edge((3, 2), (5, 0))
        solution.add_edge((5, 1), (6, 0))

        question = fragment("CCC(C)NCCNC", node_defs, returnHandler="exhaustive")

        assert solution in question
        assert len(question) == 1

        # Test 3.5 returnHandler "exhaustive" multiple solutions
        node_defs = set()
        node_defs.add(Group("amine", "N", [0, 0, 0]))  # Amine group
        node_defs.add(
            Group("alkene_secondary_amine", "C(N(C))", [0, 0, 0, 1, 2, 2, 2])
        )  # can be made of amine and alkene
        node_defs.add(Group("alkene", "C", [0, 0, 0]))

        solution = GroupGraph()
        solution.add_node("alkene_secondary_amine", "C(N(C))", [0, 0, 0, 1, 2, 2, 2])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("amine", "N", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_edge((0, 1), (1, 0))
        solution.add_edge((1, 1), (2, 0))
        solution.add_edge((2, 1), (3, 0))
        solution.add_edge((3, 1), (4, 0))
        solution.add_edge((3, 2), (5, 0))
        solution.add_edge((5, 1), (6, 0))

        question = fragment("CCC(C)NCCNC", node_defs, returnHandler="exhaustive")

        assert solution in question
        assert len(question) == 9

        # Test 3.6 returnHandler "exhaustive" with smiles
        node_defs = ["O", "C(=O)O", "N", "C(N(C))", "C"]

        solution = GroupGraph()
        solution.add_node("alkene_secondary_amine", "C(N(C))", [0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("amine", "N", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_node("alkene", "C", [0, 0, 0])
        solution.add_edge((0, 1), (1, 0))
        solution.add_edge((1, 1), (2, 0))
        solution.add_edge((2, 1), (3, 0))
        solution.add_edge((3, 1), (4, 0))
        solution.add_edge((3, 2), (5, 0))
        solution.add_edge((5, 1), (6, 0))

        question = fragment("CCC(C)NCCNC", node_defs, returnHandler="exhaustive")

        assert solution in question
        assert len(question) == 9
