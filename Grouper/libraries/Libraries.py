"""Module for pre-loading libraries of GroupGraphs.

This module defines classes for managing and querying libraries commonly used in group contribution methods.
It includes predefined libraries for Joback, UNIFAC, and SaftGammaMie methods.
"""

from collections import namedtuple
from typing import Dict, Union

import rdkit.Chem
import rdkit.Chem.Draw

from Grouper import Group

# Declaring namedtuple for defining adding nodes and their SMARTS to a library.
GroupExtension = namedtuple(
    "GroupExtension", ["node", "doi", "extended_smarts", "priority"]
)


class BasisSet(object):
    """Base class for managing a library of chemical groups.

    Attributes:
        node_traces_ (list): A list of GroupExtension objects representing the nodes in the library.
    """

    def __init__(self):
        """Initialize an empty BasisSet."""
        self.node_traces_ = []

    def __repr__(self):
        """Return a string representation of the BasisSet."""
        init_message = (
            f"{self.__class__.__name__} has {len(self.node_traces_)} unique nodes."
        )
        if self.doi and self.ref:
            return (
                init_message
                + f"\nSee Groups are defined from {self.ref}, and more information about the Fragmentation System can be found at {self.doi}"
            )
        return init_message

    def add_node(
        self, node: Group, doi: str = None, smarts: str = None, priority: int = None
    ):
        """Add a node to the library.

        Args:
            node (Group): The chemical group to add.
            doi (str, optional): DOI reference for the group. Defaults to None.
            smarts (str, optional): SMARTS pattern for the group. Defaults to None.
            priority (int, optional): Priority of the group. Defaults to None.
        """
        trace = GroupExtension(node, doi, smarts, priority)
        self.node_traces_.append(trace)

    def query_nodes(self, query: Dict[str, Union[str, int]]):
        """Query the library for nodes that match the query.

        Args:
            query (dict): A dictionary of attributes and their values to match.

        Returns:
            list: A list of nodes that match the query.
        """
        matched_nodes = []
        for node in self.node_traces:
            if any([getattr(node, var) == value for var, value in query.items()]):
                matched_nodes.append(node)
        return matched_nodes

    def get_nodes(self):
        """Return an iterator of the nodes in the library.

        Yields:
            Group: The next node in the library.
        """
        for trace in self.node_traces:
            yield trace.node

    def visualize_library(self):
        """Visualize the library using RDKit.

        Note:
            This method is not yet implemented.
        """
        mols = []
        for node in self.get_nodes():
            mol = rdkit.Chem.MolFromSmarts(node.pattern)
            if mol:
                mols.append(mol)
        return rdkit.Chem.Draw.MolsToGridImage(
            mols,
            molsPerRow=10,
            subImgSize=(200, 200),
            legends=[node.type for node in self.get_nodes()],
        )

    @property
    def node_traces(self):
        """Get the list of node traces."""
        return self.node_traces_

    @node_traces.setter
    def node_traces(self, value):
        """Set the list of node traces.

        Args:
            value (list): A list of GroupExtension objects.
        """
        self.node_traces_ = value

    @property
    def n_nodes(self):
        """Get the number of nodes in the library.

        Returns:
            int: The number of nodes.
        """
        return len(self.node_traces_)


class Joback(BasisSet):
    """Joback group contribution method for estimating physical properties of organic compounds.

    This class predefines a library of chemical groups based on the Joback method.
    """

    def __init__(self):
        """Initialize the Joback library with predefined nodes."""
        super().__init__()
        doi = "https://doi.org/10.1021/i200013a024"
        ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/Joback.jl"
        self.node_traces = [
            GroupExtension(Group("-CH3", "C", [0]), "", "[CX4H3]", None),
            GroupExtension(Group("-CH2-", "C", [0, 0]), "", "[!R;CX4H2]", None),
            GroupExtension(Group(">CH-", "C", [0, 0, 0]), "", "[!R;CX4H]", None),
            GroupExtension(Group(">C<", "C", [0, 0, 0, 0]), "", "[!R;CX4H0]", None),
            GroupExtension(Group("CH2=CH-", "C=C", [0]), "", "[CX3H2][CX3H1]", None),
            GroupExtension(Group("-CH=CH-", "C=C", [0, 1]), "", "[CX3H1][CX3H1]", None),
            GroupExtension(
                Group("=C<", "C", [0, 0, 0]),
                "",
                "[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]",
                None,
            ),
            GroupExtension(Group("=C=", "C", [0, 0]), "", "[$([CX2H0](=*)=*)]", None),
            GroupExtension(Group("CH", "C", [0]), "", "[$([CX2H1]#[!#7])]", None),
            GroupExtension(Group("C", "C", [0, 0]), "", "[$([CX2H0]#[!#7])]", None),
            GroupExtension(Group("ring-CH2-", "C", [0, 0]), "", "[R;CX4H2]", None),
            GroupExtension(Group("ring>CH-", "C", [0, 0, 0]), "", "[R;CX4H]", None),
            GroupExtension(Group("ring>C<", "C", [0, 0, 0, 0]), "", "[R;CX4H0]", None),
            GroupExtension(Group("ring=CH-", "C", [0, 0]), "", "[R;CX3H1,cX3H1]", None),
            GroupExtension(
                Group("ring=C<", "C", [0, 0, 0]),
                "",
                "[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]",
                None,
            ),
            GroupExtension(Group("-F", "F", [0]), "", "[F;X1]", None),
            GroupExtension(Group("-Cl", "[Cl]", [0], True), "", "[Cl;X1]", None),
            GroupExtension(Group("-Br", "[Br]", [0], True), "", "[Br;X1]", None),
            GroupExtension(Group("-I", "I", [0]), "", "[I;X1]", None),
            GroupExtension(
                Group("-OH (alcohol)", "O", [0]),
                "",
                "[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]",
                None,
            ),
            GroupExtension(
                Group("-OH (phenol)", "O", [0]), "", "[O;H1;$(O-!@c)]", None
            ),
            GroupExtension(
                Group("-O- (non-ring)", "O", [0, 0]),
                "",
                "[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]",
                None,
            ),
            GroupExtension(
                Group("-O- (ring)", "O", [0, 0]),
                "",
                "[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]",
                None,
            ),
            GroupExtension(
                Group(">C=O (non-ring)", "C=O", [0, 0]),
                "",
                "[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O",
                None,
            ),
            GroupExtension(
                Group(">C=O (ring)", "C=O", [0, 0]),
                "",
                "[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O",
                None,
            ),
            GroupExtension(
                Group("O=CH- (aldehyde)", "O=C", [1]), "", "[CH;D2](=O)", None
            ),
            GroupExtension(
                Group("-COOH (acid)", "C(=O)O", [0]), "", "[OX2H]-[C]=O", None
            ),
            GroupExtension(
                Group("-COO- (ester)", "C(=O)O", [0, 2]),
                "",
                "[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]",
                None,
            ),
            GroupExtension(
                Group("=O (other than above)", "O", [0]),
                "",
                "[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]",
                None,
            ),
            GroupExtension(Group("-NH2", "N", [0]), "", "[NX3H2]", None),
            GroupExtension(
                Group(">NH (non-ring)", "N", [0, 0]), "", "[NX3H1;!R]", None
            ),
            GroupExtension(Group(">NH (ring)", "N", [0, 0]), "", "[#7X3H1;R]", None),
            GroupExtension(
                Group(">N- (non-ring)", "N", [0, 0, 0]),
                "",
                "[#7X3H0;!$([#7](~O)~O)]",
                None,
            ),
            GroupExtension(
                Group("-N= (non-ring)", "N", [0, 0]), "", "[#7X2H0;!R]", None
            ),
            GroupExtension(Group("-N= (ring)", "N", [0, 0]), "", "[#7X2H0;R]", None),
            GroupExtension(Group("=NH", "N", [0]), "", "[#7X2H1]", None),
            GroupExtension(Group("-CN", "CN", [0]), "", "[#6X2]#[#7X1H0]", None),
            GroupExtension(
                Group("-NO2", "[N+]([O-])[O]", [0], True),
                "",
                "[$([#7X3,#7X3+][!#8])](=[O])~[O-]",
                None,
            ),
            GroupExtension(Group("-SH", "S", [0]), "", "[SX2H]", None),
            GroupExtension(
                Group("-S- (non-ring)", "S", [0, 0]), "", "[#16X2H0;!R]", None
            ),
            GroupExtension(Group("-S- (ring)", "S", [0, 0]), "", "[#16X2H0;R]", None),
        ]


class Unifac(BasisSet):
    """UNIFAC group contribution method for estimating activity coefficients.

    This class predefines a library of chemical groups based on the UNIFAC method.
    """

    def __init__(self):
        """Initialize the UNIFAC library with predefined nodes."""
        super().__init__()
        doi = "https://doi.org/10.1021/i200013a024"
        ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/UNIFAC.jl"
        self.node_traces = [
            GroupExtension(Group("CH3", "C", [0]), "", "[CX4;H3;!R]", None),
            GroupExtension(Group("CH2", "C", [0, 0]), "", "[CX4;H2;!R]", None),
            GroupExtension(Group("CH", "C", [0, 0, 0]), "", "[CX4;H1;!R]", None),
            GroupExtension(Group("C", "C", [0, 0, 0, 0]), "", "[CX4;H0;!R]", None),
            GroupExtension(Group("CH2=CH", "C=C", [0]), "", "[CX3;H2]=[CX3;H1]", None),
            GroupExtension(
                Group("CH=CH", "C=C", [0, 1]), "", "[CX3;H1]=[CX3;H1]", None
            ),
            GroupExtension(
                Group("CH2=C", "C=C", [0, 0]), "", "[CX3;H2]=[CX3;H0]", None
            ),
            GroupExtension(
                Group("CH=C", "C=C", [0, 0, 1]), "", "[CX3;H1]=[CX3;H0]", None
            ),
            GroupExtension(
                Group("OH (P)", "O", [0]), "", "[OH1;$([OH1][CX4H2])]", None
            ),
            GroupExtension(Group("CH3OH", "CO", [0]), "", "[CX4;H3][OX2;H1]", None),
            GroupExtension(
                Group("ACOH", "C=O", [0, 0]), "", "[cX3;H0;R][OX2;H1]", None
            ),
            GroupExtension(
                Group("CH2CO", "CC=O", [0]), "", "[CX4;H2][CX3;!H1](=O)", None
            ),
            GroupExtension(Group("CHO", "C=O", [0, 0]), "", "[CX3H1](=O)", None),
            GroupExtension(
                Group("CH3COO", "CC(=O)O", [3]), "", "[CH3][CX3;H0](=[O])[OH0]", None
            ),
            GroupExtension(
                Group("CH2COO", "CC(=O)O", [0, 3]),
                "",
                "[CX4;H2][CX3](=[OX1])[OX2]",
                None,
            ),
            GroupExtension(
                Group("HCOO", "C(=O)O", [0, 2]), "", "[CX3;H1](=[OX1])[OX2]", None
            ),
            GroupExtension(Group("CH3O", "CO", [0, 1]), "", "[CH3;!R][OH0;!R]", None),
            GroupExtension(Group("CH2O", "CO", [0, 1]), "", "[CH2;!R][OH0;!R]", None),
            GroupExtension(Group("CHO", "CO", [0, 1]), "", "[C;H1;!R][OH0;!R]", None),
            GroupExtension(Group("ACH", "C=O", [0]), "", "[cX3;H1]", None),
            GroupExtension(Group("AC", "C=O", [0]), "", "[cX3;H0]", None),
            GroupExtension(Group("ACCH3", "C(=O)C", [0]), "", "[cX3;H0][CX4;H3]", None),
            GroupExtension(
                Group("ACCH2", "C(=O)C", [0, 2]), "", "[cX3;H0][CX4;H2]", None
            ),
            GroupExtension(
                Group("ACCH", "C(=O)C", [0, 2, 2]), "", "[cX3;H0][CX4;H1]", None
            ),
            GroupExtension(Group("CH2NH2", "CN", [0]), "", "[CX4;H2][NX3;H2]", None),
            GroupExtension(Group("CHNH2", "CN", [0, 0]), "", "[CX4;H1][NX3;H2]", None),
            GroupExtension(Group("CH3NH", "CN", [1]), "", "[CX4;H3][NX3;H1]", None),
            GroupExtension(Group("CH2NH", "CN", [0, 1]), "", "[CX4;H2][NX3;H1]", None),
            GroupExtension(
                Group("CHNH", "CN", [0, 0, 1]), "", "[CX4;H1][NX3;H1]", None
            ),
            GroupExtension(Group("CH3N", "CN", [1, 1]), "", "[CX4;H3][NX3;H0]", None),
            GroupExtension(
                Group("CH2N", "CN", [0, 1, 1]), "", "[CX4;H2][NX3;H0]", None
            ),
            GroupExtension(Group("ACNH2", "C(=O)N", [0]), "", "[c][NX3;H2]", None),
            GroupExtension(
                Group("AC2H2N", "C(=O)NC=O", [2]), "", "[cX3H1][n][cX3H1]", None
            ),
            GroupExtension(
                Group("AC2HN", "C(=O)NC=O", [0, 2]), "", "[cX3H0][n][cX3H1]", None
            ),
            GroupExtension(
                Group("AC2N", "C(=O)NC=O", [0, 2, 3]), "", "[cX3H0][n][cX3H0]", None
            ),
            GroupExtension(
                Group("CH2CN", "CC#N", [0]), "", "[CX4;H2][CX2]#[NX1]", None
            ),
            GroupExtension(
                Group("COO", "C(=O)O", [0, 1]),
                "",
                "[CX3,cX3](=[OX1])[OX2H0,oX2H0]",
                None,
            ),
            GroupExtension(
                Group("COOH", "C(=O)O", [0]), "", "[CX3](=[OX1])[O;H1]", None
            ),
            GroupExtension(
                Group("CH2CL", "C([Cl])", [0], True),
                "",
                "[CX4;H2;!$(C(Cl)(Cl))](Cl)",
                None,
            ),
            GroupExtension(
                Group("CHCL", "C([Cl])", [0, 0], True),
                "",
                "[CX4;H1;!$(C(Cl)(Cl))](Cl)",
                None,
            ),
            GroupExtension(
                Group("CCL", "C([Cl])", [0, 0, 0], True), "", "[CX4;H0](Cl)", None
            ),
            GroupExtension(
                Group("CH2CL", "C([Cl])", [0], True),
                "",
                "[CX4;H2;!$(C(Cl)(Cl))](Cl)",
                None,
            ),
            GroupExtension(
                Group("CHCL", "C([Cl])", [0, 0], True),
                "",
                "[CX4;H1;!$(C(Cl)(Cl))](Cl)",
                None,
            ),
            GroupExtension(
                Group("CCL", "C([Cl])", [0, 0, 0], True), "", "[CX4;H0](Cl)", None
            ),
            GroupExtension(
                Group("CHCL2", "C([Cl])([Cl])", [0], True),
                "",
                "[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)",
                None,
            ),
            GroupExtension(
                Group("CCL2", "C([Cl])([Cl])", [0, 0], True),
                "",
                "[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)",
                None,
            ),
            GroupExtension(
                Group("CCL3", "C([Cl])([Cl])([Cl])", [0], True),
                "",
                "[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)",
                None,
            ),
            GroupExtension(
                Group("ACCL", "C(=O)([Cl])", [0], True), "", "[c][Cl]", None
            ),
            GroupExtension(
                Group("CH2NO2", "C[N+](=O)[O-]", [0], True),
                "",
                "[CX4;H2][NX3](=[OX1])([OX1])",
                None,
            ),
            GroupExtension(
                Group("CHNO2", "C[N+](=O)[O-]", [0, 0], True),
                "",
                "[CX4;H1][NX3](=[OX1])([OX1])",
                None,
            ),
            GroupExtension(
                Group("ACNO2", "C(=O)[N+](=O)[O-]", [0], True),
                "",
                "[cX3][NX3](=[OX1])([OX1])",
                None,
            ),
            GroupExtension(Group("CH2SH", "CS", [0]), "", "[SX2H][CX4;H2]", None),
            GroupExtension(Group("I", "I", [0]), "", "[I]", None),
            GroupExtension(Group("BR", "[Br]", [0], True), "", "[Br]", None),
            GroupExtension(Group("CH=-C", "C#C", [0]), "", "[CX2;H1]#[CX2;H0]", None),
            GroupExtension(Group("C=-C", "C#C", [0, 1]), "", "[CX2;H0]#[CX2;H0]", None),
            GroupExtension(
                Group("CL-(C=C)", "[Cl][C]=[C]", [1, 2, 2], True),
                "",
                "[$([Cl;H0]([C]=[C]))]",
                None,
            ),
            GroupExtension(
                Group("C=C", "C=C", [0, 0, 1, 1]), "", "[CX3;H0]=[CX3;H0]", None
            ),
            GroupExtension(Group("ACF", "C(=O)F", [0]), "", "[cX3][F]", None),
            GroupExtension(Group("CF3", "C(F)(F)F", [0]), "", "C(F)(F)F", None),
            GroupExtension(Group("CF2", "C(F)F", [0, 0]), "", "C(F)F", None),
            GroupExtension(Group("CF", "CF", [0, 0, 0]), "", "C(F)", None),
            GroupExtension(Group("CY-CH2", "C", [0, 0]), "", "[CH2;R]", None),
            GroupExtension(Group("CY-CH", "C", [0, 0, 0]), "", "[CH1;R]", None),
            GroupExtension(Group("CY-C", "C", [0, 0, 0, 0]), "", "[CH0;R]", None),
            GroupExtension(
                Group("OH (S)", "CO", [0, 0]), "", "[OH1;$([OH1][CX4H1])]", None
            ),
            GroupExtension(
                Group("OH (T)", "CO", [0, 0, 0]), "", "[OH1;$([OH1][CX4H0])]", None
            ),
            GroupExtension(
                Group("CY-CH2O", "CO", [0, 1]),
                "",
                "[CX4H2;R][OX2;R;$(O(CC)C)][CX4H2;R][OX2;R][CX4H2;R]",
                None,
            ),
            GroupExtension(Group("CNH2", "CN", [0, 0, 0]), "", "[CX4H0][NH2]", None),
            GroupExtension(
                Group("CONH2", "C(=O)N", [0]), "", "[CX3H0](=[OX1H0])[NX3H2]", None
            ),
            GroupExtension(
                Group("CONHCH3", "C(=O)NC", [0]),
                "",
                "[OX1H0;!R]=[CX3H0;!R][NH1X3;!R][CH3;!R]",
                None,
            ),
            GroupExtension(
                Group("CONHCH2", "C(=O)NC", [0, 3]),
                "",
                "[CH2X4;!R][NH1X3;!R][CX3H0;!R]=[OX1H0;!R]",
                None,
            ),
        ]


class SaftGammaMie(BasisSet):
    """SAFT-gamma-Mie group contribution method for thermodynamic modeling.

    This class is a placeholder for the SAFT-gamma-Mie library.
    """

    def __init__(self):
        """Initialize the SAFT-gamma-Mie library.

        Note:
            This library is not fully implemented yet.
        """
        super().__init__()
        self.doi = "https://doi.org/10.1080/00268978800101601"
        self.ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/SAFTgammaMie.jl"
        self.node_traces = [
            GroupExtension(Group("-CH3", "[CH3]", [0], True), "", "[CX4H3]", None),
        ]
        print(
            "This library is not fully implemented yet... please raise a GitHub issue if you wish to use this."
        )
