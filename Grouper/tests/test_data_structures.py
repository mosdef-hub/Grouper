import pickle
import re

import pytest
from rdkit import Chem

from Grouper import Atom, AtomGraph, Group, GroupGraph
from Grouper._Grouper import GrouperParseException
from Grouper.tests.base_test import BaseTest


class TestAtom(BaseTest):
    def test_atom_initialization(self):
        atom = Atom("C", 4)
        assert atom.type == "C"
        assert atom.valency == 4

    def test_atom_comparison(self):
        atom1 = Atom("C", 4)
        atom2 = Atom("C", 4)
        atom3 = Atom("N", 3)
        assert atom1 == atom2
        assert atom1 != atom3


class TestGroup(BaseTest):
    def test_group_hash(self):
        g1 = Group("C", "[C]", [0], "SMARTS")
        g2 = Group("C", "[C]", [0], "SMARTS")
        g3 = Group("C2", "[C]", [0], "SMARTS")

        assert hash(g1) == hash(g2)
        assert hash(g1) != hash(g3)

    def test_brackets(self):
        gG = GroupGraph()
        g = Group("cl", "[Cl]C=C[Br]", [1, 2], "SMARTS")
        gG.add_node(g.type, g.pattern, g.hubs, g.pattern_type)
        assert (gG.to_smiles() == "[Cl]C=C[Br]") or (gG.to_smiles() == "ClC=CBr")

    def test_charged_species(self):
        # Existing tests
        Group("CH2NO2", "C[N+](=O)[O-]", [0], "SMARTS")
        msg = r"Adding edge from 0 to 1 would exceed the valency for the source node"
        with pytest.raises(ValueError, match=re.escape(msg)):
            Group("cl", "[Cl-]C=C[Br]", [0], "SMARTS")
        Group("cl", "[Cl-]C=C[Br]", [0], "NONATOMIC")
        Group("cl", "[Cl+]C=C[Br]", [0], "NONATOMIC")
        Group("nat", "[Na+]", [0], "SMARTS")

        # New tests for aromatic compounds and different salts
        Group("benz", "c1ccccc1[N+](=O)[O-]", [0], "SMILES")  # Aromatic nitrobenzene
        msg = (
            r"Adding edge from 5 to 6 would exceed the valency for the destination node"
        )
        with pytest.raises(ValueError, match=re.escape(msg)):
            Group("pyrid", "c1ccncc1[Cl-]", [0], "SMARTS")
        Group(
            "pyrid", "c1ccncc1[Cl-]", [0], "NONATOMIC"
        )  # Aromatic pyridinium chloride
        msg = (
            r"Adding edge from 0 to 1 would exceed the valency for the destination node"
        )
        with pytest.raises(ValueError, match=re.escape(msg)):
            Group("ammonium", "[NH4+][Cl-]", [0], "SMARTS")  # Ammonium chloride
        msg = r"Invalid SMARTS: [NH4+].[Cl-] with detached molecules."
        with pytest.raises(GrouperParseException, match=re.escape(msg)):
            Group("ammonium", "[NH4+].[Cl-]", [0], "SMARTS")
        Group("ammonium", "[NH4+][Cl-]", [0], "NONATOMIC")
        Group("sodium_acetate", "[Na+][O-]C=O", [0], "NONATOMIC")  # Sodium acetate
        Group(
            "potassium_permanganate", "[K+][MnO4-]", [0], "NONATOMIC"
        )  # Potassium permanganate

    def test_group_initialization(self):
        group = Group("C", "[C]", [0], "SMARTS")
        assert group.type == "C"
        assert group.pattern == "[C]"
        assert group.hubs == [0]
        assert group.pattern_type == "SMARTS"

    def test_group_to_string(self):
        group = Group("carbon", "C", [0, 0, 0, 0])
        assert (
            str(group)
            == "Group  (carbon) (C) : \n    Ports 0 1 2 3 \n    Hubs  0 0 0 0 "
        )


class TestGroupGraph(BaseTest):
    def to_set_of_sets(self, matches):
        """Helper function for set comparison."""
        return {frozenset(match) for match in matches}

    def test_group_equality(self):
        g1 = Group("C", "[C]", [0], "SMARTS")
        g2 = Group("C", "[C]", [0], "SMARTS")
        g3 = Group("C2", "[C]", [0], "SMARTS")

        assert g1 == g2
        assert g1 != g3

        gg = GroupGraph()
        gg.add_node("O", "[O]", [0], "SMARTS")
        gg.add_node("O", "[O]", [0], "SMARTS")
        gg.add_node("O", "[O]", [0], "SMARTS")
        gg.add_node("O", "[O]", [0], "SMARTS")

        assert all(
            [n1 == n2 for n1, n2 in zip(gg.nodes.values(), list(gg.nodes.values())[1:])]
        )

        gg = GroupGraph()
        gg.add_node("O", "[O]", [0], "SMARTS")
        gg.add_node("C", "[O]", [0], "SMARTS")
        gg.add_node("N", "[O]", [0], "SMARTS")
        gg.add_node("S", "[O]", [0], "SMARTS")

        assert not any([n1 == n2 for n1, n2 in zip(gg.nodes, list(gg.nodes)[1:])])

    def test_add_node(self):
        # Basic node addition
        graph = GroupGraph()

        with pytest.raises(ValueError):  # No smarts
            graph.add_node("type1", "", [0, 0])

        with pytest.raises(ValueError):  # Invalid hubs
            graph.add_node("type1", "C", [0, -1])

        with pytest.raises(ValueError):  # Invalid hubs
            graph.add_node("type1", "C", [0, 1, 2])

        with pytest.raises(GrouperParseException):  # Invalid smarts
            graph.add_node("type1", "asldkfghj", [0])

        with pytest.raises(ValueError):  # No type
            graph.add_node("", "C", [0, 0])

        assert len(graph.nodes) == 0

        graph.add_node("type1", "C", [0, 0])

        assert len(graph.nodes) == 1

        assert set(n.type for n in graph.nodes.values()) == {"type1"}
        assert set(n.pattern for n in graph.nodes.values()) == {"C"}
        assert set(tuple(n.ports) for n in graph.nodes.values()) == {(0, 1)}
        assert set(tuple(n.hubs) for n in graph.nodes.values()) == {(0, 0)}

        graph.add_node("type1")

        # Adding a node with different type and pattern
        graph.add_node("type2", "C", [0])

        assert len(graph.nodes) == 3
        assert set(n.type for n in graph.nodes.values()) == {"type1", "type1", "type2"}
        assert set(n.pattern for n in graph.nodes.values()) == {"C"}
        assert set(tuple(n.ports) for n in graph.nodes.values()) == {(0,), (0, 1)}
        assert set(tuple(n.hubs) for n in graph.nodes.values()) == {(0,), (0, 0)}

        group = Group("alkene", "C=C", [0, 0, 1, 1])
        graph.add_node(group)

        assert len(graph.nodes) == 4
        assert set(n.type for n in graph.nodes.values()) == {
            "type1",
            "type1",
            "type2",
            "alkene",
        }
        assert set(n.pattern for n in graph.nodes.values()) == {"C", "C=C"}
        assert set(tuple(n.ports) for n in graph.nodes.values()) == {
            (0,),
            (0, 1),
            (0, 1, 2, 3),
        }
        assert set(tuple(n.hubs) for n in graph.nodes.values()) == {
            (0,),
            (0, 0),
            (0, 0, 1, 1),
        }

    def test_add_no_port_group(self):
        graph = GroupGraph()
        group = Group(ntype="solo", pattern="c1ccccc1", pattern_type="SMARTS")
        graph.add_node(group)
        group = Group(ntype="solo", pattern="c1ccccc1", hubs=[], pattern_type="SMARTS")
        graph.add_node(group)
        graph.add_node(group.type, group.pattern, [], group.pattern_type)
        assert len(graph.nodes) == 3
        msg = re.escape(r"Source port does not exist")
        with pytest.raises(ValueError, match=msg):
            graph.add_edge((0,0), (1,0))

    def test_add_edge(self):
        graph = GroupGraph()
        graph.add_node("type1", "C", [0, 0])
        graph.add_node("type2", "C", [0, 0])
        graph.add_node("type1", "C", [0, 0])

        graph.add_edge((0, 0), (1, 0))
        assert (0, 0, 1, 0, 1) in graph.edges

        graph.add_edge((2, 1), (1, 1))
        assert (2, 1, 1, 1, 1) in graph.edges

    def test_add_edge_bond_order(self):
        graph = GroupGraph()
        graph.add_node("type1", "C", [0, 0, 0, 0])
        graph.add_node("type2", "C", [0, 0, 0, 0])
        graph.add_node("type3", "C", [0, 0, 0, 0])

        graph.add_edge((0, 0), (1, 0), 2)
        assert (0, 0, 1, 0, 2) in graph.edges

        graph.add_edge((2, 1), (1, 1), 1.5)
        assert (2, 1, 1, 1, 1.5) in graph.edges

        graph.add_edge((2, 2), (1, 2), 1)
        assert (2, 2, 1, 2, 1) in graph.edges

    def test_add_edge_with_invalid_nodes(self):
        graph = GroupGraph()
        graph.add_node("node1", "C", [0, 0])
        graph.add_node("node2", "C", [0])
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (2, 1))

    def test_add_edge_with_invalid_ports(self):
        graph = GroupGraph()
        graph.add_node("node1", "C", [0, 0])
        graph.add_node("node2", "C", [0, 0])
        with pytest.raises(ValueError):
            graph.add_edge((0, 2), (1, 1))

    def test_add_edge_with_occupied_port(self):
        graph = GroupGraph()
        graph.add_node("node1", "C", [0, 0])
        graph.add_node("node2", "C", [0, 0])
        graph.add_edge((0, 1), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (1, 0))

    def test_add_edge_with_same_ports(self):
        graph = GroupGraph()
        graph.add_node("node1", "C", [0, 0])
        graph.add_node("node2", "C", [0, 0])
        graph.add_node("node3", "C", [0, 0])
        graph.add_edge((0, 1), (1, 1))
        graph.add_edge((1, 0), (2, 1))
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (1, 0))
        with pytest.raises(ValueError):
            graph.add_edge((0, 0), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((1, 1), (2, 0))
        with pytest.raises(ValueError):
            graph.add_edge((2, 0), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((1, 1), (2, 1))


        graph = GroupGraph()
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_edge((0, 1), (1, 1))
        graph.add_edge((0, 0), (2, 0))
        with pytest.raises(ValueError):
            graph.add_edge((1, 1), (2, 1))

    def test_equal(self):
        graph1 = GroupGraph()
        graph2 = GroupGraph()
        assert graph1 == graph2
        graph1.add_node("node1", "C", [0, 0])
        assert graph1 != graph2
        graph2.add_node("node1", "C", [0, 0])
        assert graph1 == graph2
        graph1.add_node("node2", "C", [0, 0])
        assert graph1 != graph2
        graph2.add_node("node2", "C", [0, 0])
        assert graph1 == graph2
        graph1.add_edge((0, 0), (1, 0))
        assert graph1 != graph2
        graph2.add_edge((0, 0), (1, 0))
        assert graph1 == graph2

    def test_in(self):
        graph1 = GroupGraph()
        graph2 = GroupGraph()

        assert graph1 in [graph1]

        graph1.add_node("node1", "C", [0, 0])
        assert graph1 not in [graph2]

        graph2.add_node("node1", "C", [0, 0])
        assert graph1 in [graph2]

        graph1.add_node("node2", "C", [0, 0])
        assert graph1 not in [graph2]

        graph2.add_node("node2", "C", [0, 0])
        assert graph1 in [graph2]

        graph1.add_edge((0, 0), (1, 0))
        assert graph1 not in [graph2]

        graph2.add_edge((0, 0), (1, 0))
        assert graph1 in [graph2]

        assert graph1 in [graph1, graph2]

    def test_rings(self):
        graph = GroupGraph()
        graph.add_node("ring", "c1ccccc1", [0, 1, 2, 3, 4, 5])
        graph.add_node("ring")
        graph.add_node("ring")

        graph.add_edge((0, 0), (1, 0))
        with pytest.raises(ValueError):
            graph.add_edge((0, 0), (2, 0))

    @pytest.mark.parametrize(
        "graph_fixture", ["empty_graph", "basic_graph", "single_node_graph"]
    )
    def test_to_atom_graph(self, request, graph_fixture):
        # Access the graph using request.getfixturevalue
        graph = request.getfixturevalue(graph_fixture)

        # Example assertions for different graphs
        if graph_fixture == "empty_graph":
            with pytest.raises(ValueError):
                graph.to_atom_graph()
        elif graph_fixture == "basic_graph":
            atomic_graph = graph.to_atom_graph()
            truth = AtomGraph()
            truth.add_node("C", 4)
            truth.add_node("N", 3)
            truth.add_edge(0, 1)
            assert atomic_graph == truth
        elif graph_fixture == "single_node_graph":
            atomic_graph = graph.to_atom_graph()
            truth = AtomGraph()
            truth.add_node("C", 1)
            assert atomic_graph == truth

    def test_n_free_ports(self):
        graph = GroupGraph()
        graph.add_node("node1", "[C]", [0, 0], "SMARTS")
        graph.add_node("node2", "[C]", [0, 0, 0], "SMARTS")
        assert graph.n_free_ports(0) == 2

        # Connect a edge and recheck
        graph.add_edge((0, 0), (1, 0))
        assert graph.n_free_ports(0) == 1
        assert graph.n_free_ports(1) == 2

    @pytest.mark.parametrize(
        "graph_fixture",
        [
            "empty_graph",
            "basic_graph",
            "single_node_graph",
            "single_edge_graph",
            "five_member_ring_graph",
        ],
    )
    def test_to_smiles(self, request, graph_fixture):
        graph = request.getfixturevalue(graph_fixture)
        if graph_fixture == "empty_graph":
            assert graph.to_smiles() == ""
        if graph_fixture == "basic_graph":
            assert graph.to_smiles() == "CN"
        if graph_fixture == "single_node_graph":
            assert graph.to_smiles() == "C"
        if graph_fixture == "single_edge_graph":
            assert graph.to_smiles() == "CC"
        if graph_fixture == "five_member_ring_graph":
            assert graph.to_smiles() == "C1CCCC1"

    def test_add_node_performance(self, benchmark):
        graph = GroupGraph()
        for i in range(100):
            graph.add_node(f"type{i}", "C", [0, 0])

        def benchmark_add_node():
            graph.add_node("type100", "C", [0, 0])

        # Benchmark the add_node method
        benchmark(benchmark_add_node)

    # def test_canonize(self):
    #     # Test simple linear structure
    #     graph1 = GroupGraph()
    #     graph1.add_node("node1", "C", [0, 0])
    #     graph1.add_node("node2", "C", [0, 0])
    #     graph1.add_node("node3", "C", [0, 0])
    #     graph1.add_edge((0, 0), (1, 0))
    #     graph1.add_edge((1, 0), (2, 0))

    #     graph2 = GroupGraph()
    #     graph2.add_node("node3", "C", [0, 0])
    #     graph2.add_node("node2", "C", [0, 0])
    #     graph2.add_node("node1", "C", [0, 0])
    #     graph2.add_edge((2, 0), (1, 0))
    #     graph2.add_edge((1, 0), (0, 0))

    #     assert graph1.to_canonical() == graph2.to_canonical()

    #     # Test with different node types
    #     graph3 = GroupGraph()
    #     graph3.add_node("node1", "CO", [0, 0, 0, 1])
    #     graph3.add_node("node2", "C", [0, 0])
    #     graph3.add_node("node3", "C", [0, 0])
    #     graph3.add_edge((0, 0), (1, 0))
    #     graph3.add_edge((1, 0), (2, 0))

    #     graph4 = GroupGraph()
    #     graph4.add_node("node1", "CO", [0, 0, 0, 1])
    #     graph4.add_node("node3", "C", [0, 0])
    #     graph4.add_node("node2", "C", [0, 0])
    #     graph4.add_edge((0, 1), (1, 0))
    #     graph4.add_edge((1, 0), (2, 0))

    #     assert graph3.to_canonical() == graph4.to_canonical()

    #     # Test cyclic structure
    #     graph5 = GroupGraph()
    #     graph5.add_node("A", "C", [0, 0])
    #     graph5.add_node("B", "C", [0, 0])
    #     graph5.add_node("C", "C", [0, 0])
    #     graph5.add_edge((0, 0), (1, 0))
    #     graph5.add_edge((1, 0), (2, 0))
    #     graph5.add_edge((2, 0), (0, 0))  # Closing the cycle

    #     graph6 = GroupGraph()
    #     graph6.add_node("B", "C", [0, 0])
    #     graph6.add_node("C", "C", [0, 0])
    #     graph6.add_node("A", "C", [0, 0])
    #     graph6.add_edge((1, 0), (2, 0))
    #     graph6.add_edge((2, 0), (0, 0))
    #     graph6.add_edge((0, 0), (1, 0))  # Different order but same cycle

    #     assert graph5.to_canonical() == graph6.to_canonical()

    #     # Test disconnected graphs with identical structure
    #     graph7 = GroupGraph()
    #     graph7.add_node("X", "C", [0, 0])
    #     graph7.add_node("Y", "C", [0, 0])

    #     graph8 = GroupGraph()
    #     graph8.add_node("Y", "C", [0, 0])
    #     graph8.add_node("X", "C", [0, 0])

    #     assert graph7.to_canonical() == graph8.to_canonical()

    #     # Test graphs with different structures
    #     graph9 = GroupGraph()
    #     graph9.add_node("X", "C", [0, 0])
    #     graph9.add_node("Y", "C", [0, 0])
    #     graph9.add_edge((0, 0), (1, 0))

    #     graph10 = GroupGraph()
    #     graph10.add_node("X", "C", [0, 0])
    #     graph10.add_node("Y", "C", [0, 0])

    #     assert (
    #         graph9.to_canonical() != graph10.to_canonical()
    #     )  # One has an edge, the other does not

    #     # Test automorphic graph
    #     graph11 = GroupGraph()
    #     graph11.add_node("X", "C", [0, 0])
    #     graph11.add_node("Y", "C", [0, 0])
    #     graph11.add_edge((0, 0), (1, 0))

    #     graph12 = GroupGraph()
    #     graph12.add_node("Y", "C", [0, 0])
    #     graph12.add_node("X", "C", [0, 0])
    #     graph12.add_edge((1, 0), (0, 0))

    def test_calculate_orbits(self):
        graph = GroupGraph()
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("ester", "C(O)=O", [0, 1])

        edge_list = [ (0, 2),(0, 3),(1, 2),(1, 3) ]
        node_colors = [0,0,0,0]
        node_orbits, edge_orbits, _ = graph.compute_orbits(edge_list, node_colors)
        assert node_orbits == [0, 0, 0, 0]
        assert edge_orbits == [0, 0, 0, 0]

        graph = GroupGraph()
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("ester", "C(O)=O", [0, 1])
        graph.add_node("not-ester", "C", [0, 0])

        edge_list = [ (0, 2),(0, 3),(1, 2),(1, 3) ]
        node_colors = [0,0,0,1]
        node_orbits, edge_orbits, _ = graph.compute_orbits(edge_list, node_colors)
        assert node_orbits == [0, 0, 2, 3]
        assert edge_orbits == [0, 1, 0, 1]


class TestAtomGraph(BaseTest):
    def to_set_of_sets(self, matches):
        """Helper function for set comparison."""
        return {frozenset(match) for match in matches}

    def test_atom_equality(self):
        a1 = Atom("C")
        a2 = Atom("C", 4)
        a3 = Atom("C", 3)
        a4 = Atom("N")

        assert a1 == a2
        assert a1 != a3
        assert a1 != a4
        assert a3 != a4

        ag = AtomGraph()
        ag.add_node("C")
        ag.add_node("C", 4)
        nodeList = list(ag.nodes.values())

        assert nodeList[0] == nodeList[1]

        ag = AtomGraph()
        ag.add_node("C")
        ag.add_node("C", 3)
        ag.add_node("N")
        ag.add_node("N", 2)
        nodeList = list(ag.nodes.values())

        assert not any([n1 == n2 for n1, n2 in zip(nodeList, nodeList[1:])])

    def test_add_node(self):
        # Basic node addition for AtomGraph
        agraph = AtomGraph()
        agraph.add_node("C", 4)
        assert set(n.type for n in agraph.nodes.values()) == {"C"}
        assert set(n.valency for n in agraph.nodes.values()) == {4}

        atom = Atom("C", 4)
        agraph.add_node(atom)
        assert set(n.type for n in agraph.nodes.values()) == {"C", "C"}
        assert set(n.valency for n in agraph.nodes.values()) == {4, 4}
        assert len(agraph.nodes) == 2
        assert len(agraph.edges) == 0

    def test_from_smiles(self):
        graph = AtomGraph()
        graph.from_smiles("CN")
        truth = AtomGraph()
        truth.add_node("C", 4)
        truth.add_node("N", 3)
        truth.add_edge(0, 1, 1)
        assert graph == truth

        graph = AtomGraph()
        graph.from_smiles("C")
        truth = AtomGraph()
        truth.add_node("C", 4)
        assert graph == truth

        graph = AtomGraph()
        graph.from_smiles("CC")
        truth = AtomGraph()
        truth.add_node("C", 4)
        truth.add_node("C", 4)
        truth.add_edge(0, 1)
        assert graph == truth

        graph = AtomGraph()
        graph.from_smiles("O=COOC(=O)O")
        truth = AtomGraph()
        truth.add_node("O", 2)
        truth.add_node("C", 4)
        truth.add_node("O", 2)
        truth.add_node("O", 2)
        truth.add_node("C", 4)
        truth.add_node("O", 2)
        truth.add_node("O", 2)
        truth.add_edge(0, 1, 2)
        truth.add_edge(1, 2, 1)
        truth.add_edge(2, 3, 1)
        truth.add_edge(3, 4, 1)
        truth.add_edge(4, 5, 2)
        truth.add_edge(4, 6, 1)
        assert graph == truth

    @pytest.mark.parametrize(
        "smarts",
        [
            "CCCC",
            "[C](N)O",
            "[C](N)O(C)",
            "[C](N)OC",
            "[C](O(N))CC",
            "[C](CC(N(C))C)O",
            "[C](C1C(N(C))CCCC1)O",
        ],
    )
    def test_from_smarts(self, smarts):
        # missing (8,1,1) Extra (0,0,1)
        ag = AtomGraph()
        ag.from_smarts(smarts)
        mol = Chem.MolFromSmarts(smarts)
        # check nodes
        assert len(mol.GetAtoms()) == len(ag.nodes)
        # check edges
        assert int(len(mol.GetBonds())) * 2 == len(ag.edges)

        # check connectivity
        agBonds = [
            {ag.nodes[edge1].type, ag.nodes[edge2].type} for edge1, edge2, _ in ag.edges
        ]
        for bond in mol.GetBonds():
            assert {
                bond.GetBeginAtom().GetSymbol(),
                bond.GetEndAtom().GetSymbol(),
            } in agBonds

    def test_disconnected_graph(self):
        msg = re.escape(r"Invalid SMARTS: [C].[C] with detached molecules.")
        with pytest.raises(GrouperParseException, match=msg):
            Group("a", "[C].[C]", [0, 1], "SMARTS")
        gG = GroupGraph()
        msg = re.escape(r"Invalid SMARTS: OCC[N].[N]CCO with detached molecules.")
        with pytest.raises(GrouperParseException, match=msg):
            gG.add_node("a", "OCC[N].[N]CCO", [0], "SMARTS")

    def test_substructure_search(self):
        graph = AtomGraph()
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        sub = AtomGraph()
        sub.add_node("C", 4)
        matches = graph.substructure_search(sub, [0])
        assert self.to_set_of_sets(matches) == {
            frozenset({(0, 0)}),
            frozenset({(0, 2)}),
        }
        matches = graph.substructure_search(sub, [0, 0, 0])
        assert self.to_set_of_sets(matches) == set()

        matches = graph.substructure_search(sub, [0, 0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 1)})}

        graph = AtomGraph()  # "CCOCO"
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("O", 2)
        graph.add_node("C", 4)
        graph.add_node("O", 2)
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        graph.add_edge(2, 3)
        graph.add_edge(3, 4)
        methanol = AtomGraph()  # "CO"
        methanol.add_node("C", 4)
        methanol.add_node("O", 2)
        methanol.add_edge(0, 1)
        matches = graph.substructure_search(methanol, [0, 0])
        assert self.to_set_of_sets(matches) == set()

        graph = AtomGraph()  # C=CO
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("O", 2)
        graph.add_edge(0, 1, 2)
        graph.add_edge(1, 2)
        alkene = AtomGraph()  # C=C
        alkene.add_node("C", 4)
        alkene.add_node("C", 4)
        alkene.add_edge(0, 1, 2)
        matches = graph.substructure_search(alkene, [0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 1), (1, 0)})}
        
        matches = graph.substructure_search(alkene, [0, 0, 1, 1])
        assert self.to_set_of_sets(matches) == set()
        
        matches = graph.substructure_search(alkene, [0, 0, 1])
        assert self.to_set_of_sets(matches) == set()
        
        matches = graph.substructure_search(alkene, [1, 1, 0])
        assert self.to_set_of_sets(matches) == set()
        
        matches = graph.substructure_search(alkene, [0, 0])
        assert self.to_set_of_sets(matches) == set()

        graph = AtomGraph()  # C=CO
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("O", 2)
        graph.add_edge(0, 1, 2)
        graph.add_edge(1, 2)
        oxyl = AtomGraph()
        oxyl.add_node("O", 2)
        matches = graph.substructure_search(oxyl, [0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 2)})}

    def test_substructure_search_2(self):
        truth = AtomGraph()
        truth.from_smiles("CNCNOC=O")

        sub = AtomGraph()
        sub.add_node("C", 4)
        sub.add_node("O", 2)
        sub.add_node("O", 2)
        sub.add_edge(0, 1)
        sub.add_edge(0, 2, 2)

        matches = truth.substructure_search(sub, [0])
        assert len(matches) == 0

        matches = truth.substructure_search(sub, [1])
        assert len(matches) == 1

    def test_substructure_search_3(self):
        truth = AtomGraph()
        truth.from_smiles("CNCNC(=O)O")

        sub = AtomGraph()
        sub.add_node("C", 4)
        sub.add_node("N", 3)
        sub.add_node("C", 4)
        sub.add_edge(0, 1)
        sub.add_edge(1, 2)

        matches = truth.substructure_search(sub, [2])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 0), (1, 1), (2, 2)})}

        matches = truth.substructure_search(sub, [0, 2, 2])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 2), (1, 3), (2, 4)})}

    def test_json(self):
        graph = GroupGraph()
        graph.add_node("type1", "C", [0, 0])
        graph.add_node("type1", "C", [0, 0])
        graph.add_node("type2", "C", [0, 0])
        graph.add_edge((0, 0), (1, 0))
        graph.add_edge((1, 1), (2, 0))

        json = graph.to_json()
        graph2 = GroupGraph()
        graph2.from_json(json)
        assert graph == graph2

    def test_pickle(self):
        graph = GroupGraph()
        graph.add_node("type1", "C", [0, 0])
        graph.add_node("type1", "C", [0, 0])
        graph.add_node("type2", "C", [0, 0])
        graph.add_edge((0, 0), (1, 0))
        graph.add_edge((1, 1), (2, 0))

        with open("test.pkl", "wb") as f:
            pickle.dump(graph, f)

        with open("test.pkl", "rb") as f:
            graph2 = pickle.load(f)

        assert graph == graph2

    # def test_canonize(self):
    #     # Test case 1: Simple linear chain with single bonds
    #     aGraph = AtomGraph()
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_edge(0, 1, 1)  # Single bond
    #     aGraph.add_edge(1, 2, 1)  # Single bond

    #     bGraph = AtomGraph()
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_edge(0, 1, 1)  # Single bond
    #     bGraph.add_edge(1, 2, 1)  # Single bond

    #     assert aGraph.to_canonical() == bGraph.to_canonical()

    #     # Test case 2: Linear chain with a double bond
    #     aGraph = AtomGraph()
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_edge(0, 1, 2)  # Double bond
    #     aGraph.add_edge(1, 2, 1)  # Single bond

    #     bGraph = AtomGraph()
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_edge(0, 1, 2)  # Double bond
    #     bGraph.add_edge(1, 2, 1)  # Single bond

    #     assert aGraph.to_canonical() == bGraph.to_canonical()

    #     # Test case 3: Cyclic structure with single bonds
    #     aGraph = AtomGraph()
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_edge(0, 1, 1)  # Single bond
    #     aGraph.add_edge(1, 2, 1)  # Single bond
    #     aGraph.add_edge(2, 0, 1)  # Single bond

    #     bGraph = AtomGraph()
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_edge(0, 1, 1)  # Single bond
    #     bGraph.add_edge(1, 2, 1)  # Single bond
    #     bGraph.add_edge(2, 0, 1)  # Single bond

    #     assert aGraph.to_canonical() == bGraph.to_canonical()

    #     # Test case 4: Cyclic structure with a double bond
    #     aGraph = AtomGraph()
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_edge(0, 1, 2)  # Double bond
    #     aGraph.add_edge(1, 2, 1)  # Single bond
    #     aGraph.add_edge(2, 0, 1)  # Single bond

    #     bGraph = AtomGraph()
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_edge(0, 1, 2)  # Double bond
    #     bGraph.add_edge(1, 2, 1)  # Single bond
    #     bGraph.add_edge(2, 0, 1)  # Single bond

    #     assert aGraph.to_canonical() == bGraph.to_canonical()

    #     # Test case 5: Automorphic graphs with different bond orders
    #     aGraph = AtomGraph()
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_node("C", 4)
    #     aGraph.add_edge(0, 1, 1)  # Single bond
    #     aGraph.add_edge(1, 2, 2)  # Double bond

    #     bGraph = AtomGraph()
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_node("C", 4)
    #     bGraph.add_edge(0, 2, 2)  # Double bond
    #     bGraph.add_edge(2, 1, 1)  # Single bond

    #     assert aGraph.to_canonical() == bGraph.to_canonical()



    #     molecules = [
    #         "C1CC1",
    #         "C1OO1",
    #         "O1OO1",
    #         "N1OO1",
    #         "N1NO1",
    #         "N1NN1",
    #         "C1NO1",
    #         "CCC",
    #         "CCO",
    #         "C1NN1",
    #         "COO",
    #         "C1CO1",
    #         "NNN",
    #         "NON",
    #         "NNO",
    #         "C1CN1",
    #         "NOO",
    #         "COC",
    #         "OCO",
    #         "OOO",
    #         "CNC",
    #         "NCN",
    #         "CCN",
    #         "ONO",
    #         "NCO",
    #         "CNN",
    #         "CON",
    #         "CNO"
    #     ]

    #     basis = set()

    #     for mol in molecules:
    #         graph = AtomGraph()
    #         graph.from_smiles(mol)
    #         canon = graph.to_canonical()
    #         canon = tuple(sorted(canon))
    #         basis.add(canon)
    #     assert len(basis) == len(molecules), "Not all molecules are unique after canonicalization"

    # def test_add_edge_performance(self, benchmark):
    #     graph = GroupGraph()
    #     for i in range(100):
    #         graph.add_node(f'type{i}', '', [0, 1], [0, 0])
    #     def benchmark_add_edge():
    #         graph.add_edge((0, 0), (1, 1))
    #         graph.remove_edge((0, 0), (1, 1))

    #     # Benchmark the add_edge method
    #     benchmark(benchmark_add_edge)


def TestComplexPatternedGroups(BaseTest):
    def test_error_messages(self):
        msg = r"Adding edge from 0 to 1 would exceed the valency for the source node"
        with pytest.raises(ValueError, match=re.escpape(msg)):
            Group("cl", "[Cl]C=C[Br]", [1, 2], "SMARTS")

    def test_smiles_brackets(self):
        assert Group("type1", "[C]", [0])
        assert Group("type1", "[C][C]", [0])
        assert Group("type1", "[Cl][N](C)[Li+]", [0])
        gG = GroupGraph()
        gG.add_node("type1", "C([Cl])([Br])[Li+]", [0])
        assert gG.to_smiles() == "[Li+]C(Cl)Br"

    def test_non_elemental_graphs(self):
        graph = GroupGraph()
        graph.add_node("type1", "[Cl-]C", [1, 1], "NONATOMIC")  # negative charge
        graph.add_node("type2", "[Li+2]C", [1], "NONATOMIC")  # positive charge
        graph.add_edge((0, 0), (1, 0))
        assert graph.to_smiles() == "[Li+2]CC[Cl-]"

        # Add too many ports
        matchMsg = """
        This may be a non-atomistic representation. Check your Graph using
        Graph.visualize(), to see where it fails to represent a single molecule
        using classical valencies.
        """
        with pytest.raises(GrouperParseException, match=re.escape(matchMsg)):
            graph.add_node("type3", "C", [0, 0, 0, 0, 0])
        graph.add_node("type3", "C", [0, 0, 0, 0, 0], "NONATOMIC")  # too many ports
        graph.add_edge((2, 0), (3, 0))
        graph.add_edge((0, 1), (2, 4))
        matchMsg = """
        This is non-atomistic representation. Check your Graph using
        Graph.visualize(), to see where it fails to represent a single molecule
        using classical valencies.
        """
        with pytest.raises(GrouperParseException, match=matchMsg):
            graph.to_smiles()

    def test_coarse_grained(self):
        graph = GroupGraph()
        # Try non elemental groups
        with pytest.raises(
            GrouperParseException,
            match=re.escape("Invalid SMARTS or SMILES: [LP] provided"),
        ):
            graph.add_node("type1", "[LP]", [0, 0], "SMARTS")
        with pytest.raises(
            GrouperParseException,
            match=re.escape("Invalid SMARTS or SMILES: [LP] provided"),
        ):
            graph.add_node("type1", "[LP]", [0, 0])
        msg = "SMILES character `L` not in standard element map. Try to add brackets around each element for clarity: i.e. 'Li'->'[Li]', or try setting `Grouper.Group` patternType argument to `NONATOMIC`"
        with pytest.raises(
            GrouperParseException,
            match=re.escape(msg),
        ):
            graph.add_node("type1", "LP", [0, 0])

        graph.is_coarse_grained = True
        graph.add_node("type1", "[LP]", [0], "NONATOMIC")  # skip AtomGraph stuff
        graph.add_node("type1", "C", [0, 0, 0, 0, 0], "NONATOMIC")
        graph.add_edge((0, 0), (1, 0), 2)
        # Methods that only work on atomistic graphs
        with pytest.raises(Exception):
            graph.to_smiles()
        with pytest.raises(Exception):
            graph.to_atom_graph()

        # Methods that work on both
        assert graph.to_vector()
        obj = graph.to_json()
        assert obj
        new_graph = GroupGraph()
        assert new_graph.from_json()
        assert graph.to_canonical()
