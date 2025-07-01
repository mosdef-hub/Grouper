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
        [("saftgm", 1), ("joback", 41), ("UNIFAC", 72), ("base", 0)],
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
        nodes = list(library.get_nodes())
        assert len(nodes) == 41 == library.n_nodes
