import pytest
from rdkit import Chem

from Grouper import AtomGraph, GroupGraph
from Grouper.tests.base_test import BaseTest


class TestGroupGraph(BaseTest):
    def to_set_of_sets(self, matches):
        """Helper function for set comparison."""
        return {frozenset(match) for match in matches}

    def test_add_node(self):
        # Basic node addition
        graph = GroupGraph()
        graph.add_node("type1", "", [0, 0])
        assert set(n.type for n in graph.nodes.values()) == set(["type1"])
        assert set(n.smarts for n in graph.nodes.values()) == set([""])
        assert [n.ports for n in graph.nodes.values()] == [[0, 1]]

        # Basic node addition for AtomGraph
        agraph = AtomGraph()
        agraph.add_node("C", 4)
        assert set(n.type for n in agraph.nodes.values()) == set(["C"])
        assert set(n.valency for n in agraph.nodes.values()) == set([4])

        # Adding a node with different type and smarts
        graph.add_node("", "C", [0, 0])
        assert len(graph.nodes) == 2
        assert set(n.type for n in graph.nodes.values()) == set(["type1", "C"])
        assert set(n.smarts for n in graph.nodes.values()) == set(["", "C"])
        assert [n.ports for n in graph.nodes.values()] == [[0, 1], [0, 1]]

        # Adding a node with only a type
        graph.add_node("type1")
        assert len(graph.nodes) == 3
        assert set(n.type for n in graph.nodes.values()) == set(["type1", "C", "type1"])
        assert set(n.smarts for n in graph.nodes.values()) == set(["", "C", ""])
        assert [n.ports for n in graph.nodes.values()] == [[0, 1], [0, 1], [0, 1]]

    def test_add_edge(self):
        graph = GroupGraph()
        graph.add_node("type1", "", [0, 0])
        graph.add_node("", "C", [0, 0])
        graph.add_node("type1", "", [0, 0])

        graph.add_edge((0, 0), (1, 0))
        assert (0, 0, 1, 0, 1) in graph.edges

        graph.add_edge((2, 1), (1, 1))
        assert (2, 1, 1, 1, 1) in graph.edges

    def test_add_edge_with_invalid_nodes(self):
        graph = GroupGraph()
        graph.add_node("node1", "", [0, 0])
        graph.add_node("node2", "", [0])
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (2, 1))

    def test_add_edge_with_invalid_ports(self):
        graph = GroupGraph()
        graph.add_node("node1", "", [0, 0])
        graph.add_node("node2", "", [0, 0])
        with pytest.raises(ValueError):
            graph.add_edge((0, 2), (1, 1))

    def test_add_edge_with_occupied_port(self):
        graph = GroupGraph()
        graph.add_node("node1", "", [0, 0])
        graph.add_node("node2", "", [0, 0])
        graph.add_edge((0, 1), (1, 1))
        with pytest.raises(ValueError):
            graph.add_edge((0, 1), (1, 0))

    def test_add_edge_with_same_ports(self):
        graph = GroupGraph()
        graph.add_node("node1", "", [0, 0])
        graph.add_node("node2", "", [0, 0])
        graph.add_node("node3", "", [0, 0])
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
        graph.add_node("node1", "[C]", [0, 0])
        graph.add_node("node2", "[CX4]", [0, 0, 1])
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

    def test_from_smiles(self):
        graph = AtomGraph()
        graph.from_smiles("CN")
        truth = AtomGraph()
        truth.add_node("C", 4)
        truth.add_node("N", 3)
        truth.add_edge(0, 1, 1)
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
        ag = AtomGraph()
        ag.from_smarts(smarts)
        mol = Chem.MolFromSmarts(smarts)
        # check nodes
        assert len(mol.GetAtoms()) == len(ag.nodes)
        # check edges
        assert len(mol.GetBonds()) == len(ag.edges)

        # check connectivity
        agBonds = [{edge1.type, edge2.type} for edge1, edge2 in ag.edges]
        for bond in mol.GetBonds():
            assert {
                bond.GetBeginAtom().GetSymbol(),
                bond.GetEndAtom().GetSymbol(),
            } in agBonds

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
        ag = AtomGraph()
        ag.from_smarts(smarts)
        mol = Chem.MolFromSmarts(smarts)
        # check nodes
        assert len(mol.GetAtoms()) == len(ag.nodes)
        # check edges
        assert len(mol.GetBonds()) == len(ag.edges)

        # check connectivity
        agBonds = [{edge1.type, edge2.type} for edge1, edge2 in ag.edges]
        for bond in mol.GetBonds():
            assert {
                bond.GetBeginAtom().GetSymbol(),
                bond.GetEndAtom().GetSymbol(),
            } in agBonds

    def test_add_node_performance(self, benchmark):
        graph = GroupGraph()
        for i in range(100):
            graph.add_node(f"type{i}", "", [0, 0])

        def benchmark_add_node():
            graph.add_node("type100", "", [0, 0])

        # Benchmark the add_node method
        benchmark(benchmark_add_node)

    def test_substructure_search(self):
        graph = AtomGraph()
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        sub = AtomGraph()
        sub.add_node("C", 4)
        matches = graph.substructure_search(sub, [0, 0, 0])
        assert self.to_set_of_sets(matches) == {
            frozenset({(0, 0)}),
            frozenset({(0, 2)}),
        }
        matches = graph.substructure_search(sub, [0, 0, 0])
        assert self.to_set_of_sets(matches) == {
            frozenset({(0, 0)}),
            frozenset({(0, 2)}),
        }

        matches = graph.substructure_search(sub, [0, 0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 1)})}
        matches = graph.substructure_search(sub, [0, 0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 1)})}

        matches = graph.substructure_search(sub, [0])
        assert self.to_set_of_sets(matches) == set()
        assert self.to_set_of_sets(matches) == set()

        graph = AtomGraph()  # "CCOCO"
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
        methanol = AtomGraph()  # "CO"
        methanol.add_node("C", 4)
        methanol.add_node("O", 2)
        methanol.add_edge(0, 1)
        matches = graph.substructure_search(methanol, [0, 0])
        assert self.to_set_of_sets(matches) == {
            frozenset({(0, 1), (1, 2)}),
            frozenset({(0, 3), (1, 2)}),
        }
        matches = graph.substructure_search(methanol, [0, 0])
        assert self.to_set_of_sets(matches) == {
            frozenset({(0, 1), (1, 2)}),
            frozenset({(0, 3), (1, 2)}),
        }

        graph = AtomGraph()  # C=CO
        graph = AtomGraph()  # C=CO
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("O", 2)
        graph.add_edge(0, 1, 2)
        graph.add_edge(1, 2)
        alkene = AtomGraph()  # C=C
        alkene = AtomGraph()  # C=C
        alkene.add_node("C", 4)
        alkene.add_node("C", 4)
        alkene.add_edge(0, 1, 2)
        matches = graph.substructure_search(alkene, [0])
        assert self.to_set_of_sets(matches) == set()
        matches = graph.substructure_search(alkene, [0, 0, 1, 1])
        assert self.to_set_of_sets(matches) == set()
        matches = graph.substructure_search(alkene, [0, 0, 1])
        assert self.to_set_of_sets(matches) == {frozenset({(1, 1), (0, 0)})}
        matches = graph.substructure_search(alkene, [1, 1, 0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 1), (1, 0)})}
        matches = graph.substructure_search(alkene, [0, 0])
        assert self.to_set_of_sets(matches) == set()
        assert self.to_set_of_sets(matches) == set()
        matches = graph.substructure_search(alkene, [0, 0, 1, 1])
        assert self.to_set_of_sets(matches) == set()
        matches = graph.substructure_search(alkene, [0, 0, 1])
        assert self.to_set_of_sets(matches) == {frozenset({(1, 1), (0, 0)})}
        matches = graph.substructure_search(alkene, [1, 1, 0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 1), (1, 0)})}
        matches = graph.substructure_search(alkene, [0, 0])
        assert self.to_set_of_sets(matches) == set()

        graph = AtomGraph()  # C=CO
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

        matches = truth.substructure_search(sub, [0, 0, 0, 1, 2, 2])
        matches = truth.substructure_search(sub, [0, 0, 0, 1, 2, 2])
        # matches = truth.substructure_search(sub, [2,2,2, 1, 0, 0])
        assert self.to_set_of_sets(matches) == {frozenset({(0, 0), (1, 1), (2, 2)})}
        assert self.to_set_of_sets(matches) == {frozenset({(0, 0), (1, 1), (2, 2)})}

    # def test_add_edge_performance(self, benchmark):
    #     graph = GroupGraph()
    #     for i in range(100):
    #         graph.add_node(f'type{i}', '', [0, 1], [0, 0])
    #     def benchmark_add_edge():
    #         graph.add_edge((0, 0), (1, 1))
    #         graph.remove_edge((0, 0), (1, 1))

    #     # Benchmark the add_edge method
    #     benchmark(benchmark_add_edge)
