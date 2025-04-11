import pytest

from Grouper import Group
from Grouper.libraries.Libraries import (
    BasisSet,
    GroupExtension,
    Joback,
    SaftGammaMie,
    Unifac,
)
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
        node = Group("-CH3", "[CH3]", [0], True)
        nt = library.query_nodes({"node": node})[0]

        assert nt.node == node
        assert nt.extended_smarts == "[CX4H3]"
        assert nt.doi == ""
        assert nt.priority is None

        nt = GroupExtension(node, "", "[CX4H3]", 1)
        assert nt.priority == 1

    def test_add_node(self):
        library = BasisSet()
        library.add_node(Group("-CH3", "[CH3]", [0], True), "", "[CX4H3]", None)
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
