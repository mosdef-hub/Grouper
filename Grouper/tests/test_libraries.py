import pytest

from Grouper import Node
from Grouper.libraries.Libraries import BasisSet, Libraries, NodeTrace
from Grouper.tests.base_test import BaseTest


class TestLibraries(BaseTest):
    @pytest.mark.parametrize(
        "library,n_graphs", [("saftgm", 1), ("joback", 41), ("UNIFAC", 89), ("base", 0)]
    )
    def test_build_library(self, library, n_graphs):
        if library == "base":
            library = BasisSet()
        else:
            library = Libraries[library]()
        assert library.n_nodes == n_graphs

    def test_node_trace(self):
        library = Libraries["saftgm"]()
        node = Node("-CH3", "CH3", [0])
        nt = library.query_nodes({"node": node})[0]

        assert nt.node == node
        assert nt.smarts == "[CX4H3]"
        assert nt.doi == ""
        assert nt.priority is None

        nt = NodeTrace(node, "", "[CX4H3]", 1)
        assert nt.priority == 1

    def test_add_node(self):
        library = BasisSet()
        library.add_node(Node("-CH3", "CH3", [0]), "", "[CX4H3]", None)
        assert library.n_nodes == 1

    def test_query_node(self):
        library = Libraries["joback"]()
        nt = library.query_nodes({"smarts": "[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]"})[0]
        assert nt.node == Node("=C<", "C", [0, 0, 0])

    def test_list_nodes(self):
        library = Libraries["joback"]()
        nodes = list(library.get_nodes())
        assert len(nodes) == 41 == library.n_nodes
