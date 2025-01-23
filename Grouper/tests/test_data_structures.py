import pytest

from Grouper import AtomGraph, GroupGraph
from Grouper.tests.base_test import BaseTest


class TestGroupGraph(BaseTest):
    def test_add_node(self):
        # Basic node addition
        graph = GroupGraph()
        graph.add_node("type1", "", [0, 0])
        assert set(n.type for n in graph.nodes.values()) == set(["type1"])
        assert set(n.smarts for n in graph.nodes.values()) == set([""])
        assert [n.ports for n in graph.nodes.values()] == [[0, 1]]

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
        assert (0, 0, 1, 0) in graph.edges

        graph.add_edge((2, 1), (1, 1))
        assert (2, 1, 1, 1) in graph.edges

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

    def test_add_node_performance(self, benchmark):
        graph = GroupGraph()
        for i in range(100):
            graph.add_node(f"type{i}", "", [0, 0])

        def benchmark_add_node():
            graph.add_node("type100", "", [0, 0])

        # Benchmark the add_node method
        benchmark(benchmark_add_node)

    def test_substructure_search(self):
        from Grouper import AtomGraph
        graph = AtomGraph()
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        sub = AtomGraph()
        sub.add_node("C", 4)
        matches = graph.substructure_search(sub, [0])
        set_matches = set(tuple(m) for m in matches)
        assert set_matches == {(0,), (2,)}

        matches = graph.substructure_search(sub, [0,0])
        set_matches = set(tuple(m) for m in matches)
        assert set_matches == {(1,)}

        matches = graph.substructure_search(sub, [0,0,0])
        set_matches = set(tuple(m) for m in matches)
        assert set_matches == set()

        graph = AtomGraph() # "CCOCO"
        graph.add_node("C", 4)
        graph.add_node("C", 4)
        graph.add_node("O", 2)
        graph.add_node("C", 4)
        graph.add_node("O", 2)
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        graph.add_edge(2, 3)
        graph.add_edge(3, 4)
        methanol = AtomGraph()
        methanol.add_node("C", 4)
        methanol.add_node("O", 2)
        methanol.add_edge(0, 1)
        matches = graph.substructure_search(methanol, [0])
        set_matches = set(tuple(m) for m in matches)
        assert set_matches == {(3,4)}
        ether = AtomGraph()
        ether.add_node("C", 4)
        ether.add_node("O", 2)
        ether.add_edge(0, 1)
        matches = graph.substructure_search(ether, [0,1])
        set_matches = set(tuple(m) for m in matches)




    # def test_add_edge_performance(self, benchmark):
    #     graph = GroupGraph()
    #     for i in range(100):
    #         graph.add_node(f'type{i}', '', [0, 1], [0, 0])
    #     def benchmark_add_edge():
    #         graph.add_edge((0, 0), (1, 1))
    #         graph.remove_edge((0, 0), (1, 1))

    #     # Benchmark the add_edge method
    #     benchmark(benchmark_add_edge)
