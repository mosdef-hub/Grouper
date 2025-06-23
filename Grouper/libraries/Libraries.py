"""Module for pre-loading libraries of GroupGraphs.

This module defines classes for managing and querying libraries commonly used in group contribution methods.
It includes predefined libraries for Joback, UNIFAC, and SaftGammaMie methods.
"""

from collections import namedtuple
from typing import Dict, Union

import rdkit.Chem
import rdkit.Chem.Draw

from Grouper import Group

# Declaring namedtuple for defining adding groups and their SMARTS to a library.
GroupExtension = namedtuple(
    "GroupExtension", ["group", "doi", "extended_smarts", "priority"]
)


class BasisSet(object):
    """Base class for managing a library of chemical groups.

    Attributes:
        group_traces_ (list): A list of GroupExtension objects representing the groups in the library.
    """

    def __init__(self):
        """Initialize an empty BasisSet."""
        self.group_traces_ = []

    def __repr__(self):
        """Return a string representation of the BasisSet."""
        init_message = (
            f"{self.__class__.__name__} has {len(self.group_traces_)} unique groups."
        )
        if self.doi and self.ref:
            return (
                init_message
                + f"\nSee Groups are defined from {self.ref}, and more information about the Fragmentation System can be found at {self.doi}"
            )
        return init_message

    def add_group(
        self, group: Group, doi: str = None, smarts: str = None, priority: int = None
    ):
        """Add a group to the library.

        Args:
            group (Group): The chemical group to add.
            doi (str, optional): DOI reference for the group. Defaults to None.
            smarts (str, optional): SMARTS pattern for the group. Defaults to None.
            priority (int, optional): Priority of the group. Defaults to None.
        """
        trace = GroupExtension(group, doi, smarts, priority)
        self.group_traces_.append(trace)

    def query_groups(self, query: Dict[str, Union[str, int]]):
        """Query the library for groups that match the query.

        Args:
            query (dict): A dictionary of attributes and their values to match.

        Returns:
            list: A list of groups that match the query.
        """
        matched_groups = []
        for group in self.group_traces:
            if any([getattr(group, var) == value for var, value in query.items()]):
                matched_groups.append(group)
        return matched_groups

    def get_groups(self):
        """Return an iterator of the groups in the library.

        Yields:
            Group: The next group in the library.
        """
        for trace in self.group_traces:
            yield trace.group

    def visualize_library(self):
        """Visualize the library using RDKit.

        Note:
            This method is not yet implemented.
        """
        mols = []
        for group in self.get_groups():
            mol = rdkit.Chem.MolFromSmarts(group.pattern)
            if mol:
                mols.append(mol)
        return rdkit.Chem.Draw.MolsToGridImage(
            mols,
            molsPerRow=10,
            subImgSize=(200, 200),
            legends=[group.type for group in self.get_groups()],
        )

    @property
    def group_traces(self):
        """Get the list of group traces."""
        return self.group_traces_

    @group_traces.setter
    def group_traces(self, value):
        """Set the list of group traces.

        Args:
            value (list): A list of GroupExtension objects.
        """
        self.group_traces_ = value

    @property
    def n_groups(self):
        """Get the number of groups in the library.

        Returns:
            int: The number of groups.
        """
        return len(self.group_traces_)


class Joback(BasisSet):
    """Joback group contribution method for estimating physical properties of organic compounds.

    This class predefines a library of chemical groups based on the Joback method.
    """

    def __init__(self):
        """Initialize the Joback library with predefined groups."""
        super().__init__()
        doi = "https://doi.org/10.1021/i200013a024"
        ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/Joback.jl"
        self.group_traces = [
            GroupExtension(Group("-CH3", "C", [0]), doi, "[CX4H3]", None),
            GroupExtension(Group("-CH2-", "C", [0, 0]), doi, "[!R;CX4H2]", None),
            GroupExtension(Group(">CH-", "C", [0, 0, 0]), doi, "[!R;CX4H]", None),
            GroupExtension(Group(">C<", "C", [0, 0, 0, 0]), doi, "[!R;CX4H0]", None),
            GroupExtension(Group("CH2=CH-", "C=C", [0]), doi, "[CX3H2][CX3H1]", None),
            GroupExtension(
                Group("-CH=CH-", "C=C", [0, 1]), doi, "[CX3H1][CX3H1]", None
            ),
            GroupExtension(
                Group("=C<", "C", [0, 0, 0]),
                doi,
                "[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]",
                None,
            ),
            GroupExtension(Group("=C=", "C", [0, 0]), doi, "[$([CX2H0](=*)=*)]", None),
            GroupExtension(Group("CH", "C", [0]), doi, "[$([CX2H1]#[!#7])]", None),
            GroupExtension(Group("C", "C", [0, 0]), doi, "[$([CX2H0]#[!#7])]", None),
            GroupExtension(Group("ring-CH2-", "C", [0, 0]), doi, "[R;CX4H2]", None),
            GroupExtension(Group("ring>CH-", "C", [0, 0, 0]), doi, "[R;CX4H]", None),
            GroupExtension(Group("ring>C<", "C", [0, 0, 0, 0]), doi, "[R;CX4H0]", None),
            GroupExtension(
                Group("ring=CH-", "C", [0, 0]), doi, "[R;CX3H1,cX3H1]", None
            ),
            GroupExtension(
                Group("ring=C<", "C", [0, 0, 0]),
                doi,
                "[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]",
                None,
            ),
            GroupExtension(Group("-F", "F", [0]), doi, "[F;X1]", None),
            GroupExtension(Group("-Cl", "[Cl]", [0], "SMARTS"), doi, "[Cl;X1]", None),
            GroupExtension(Group("-Br", "[Br]", [0], "SMARTS"), doi, "[Br;X1]", None),
            GroupExtension(Group("-I", "I", [0]), doi, "[I;X1]", None),
            GroupExtension(
                Group("-OH (alcohol)", "O", [0]),
                doi,
                "[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]",
                None,
            ),
            GroupExtension(
                Group("-OH (phenol)", "O", [0]), doi, "[O;H1;$(O-!@c)]", None
            ),
            GroupExtension(
                Group("-O- (non-ring)", "O", [0, 0]),
                doi,
                "[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]",
                None,
            ),
            GroupExtension(
                Group("-O- (ring)", "O", [0, 0]),
                doi,
                "[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]",
                None,
            ),
            GroupExtension(
                Group(">C=O (non-ring)", "C=O", [0, 0]),
                doi,
                "[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O",
                None,
            ),
            GroupExtension(
                Group(">C=O (ring)", "C=O", [0, 0]),
                doi,
                "[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O",
                None,
            ),
            GroupExtension(
                Group("O=CH- (aldehyde)", "O=C", [1]), doi, "[CH;D2](=O)", None
            ),
            GroupExtension(
                Group("-COOH (acid)", "C(=O)O", [0]), doi, "[OX2H]-[C]=O", None
            ),
            GroupExtension(
                Group("-COO- (ester)", "C(=O)O", [0, 2]),
                doi,
                "[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]",
                None,
            ),
            GroupExtension(
                Group("=O (other than above)", "O", [0]),
                doi,
                "[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]",
                None,
            ),
            GroupExtension(Group("-NH2", "N", [0]), doi, "[NX3H2]", None),
            GroupExtension(
                Group(">NH (non-ring)", "N", [0, 0]), doi, "[NX3H1;!R]", None
            ),
            GroupExtension(Group(">NH (ring)", "N", [0, 0]), doi, "[#7X3H1;R]", None),
            GroupExtension(
                Group(">N- (non-ring)", "N", [0, 0, 0]),
                doi,
                "[#7X3H0;!$([#7](~O)~O)]",
                None,
            ),
            GroupExtension(
                Group("-N= (non-ring)", "N", [0, 0]), doi, "[#7X2H0;!R]", None
            ),
            GroupExtension(Group("-N= (ring)", "N", [0, 0]), doi, "[#7X2H0;R]", None),
            GroupExtension(Group("=NH", "N", [0]), doi, "[#7X2H1]", None),
            GroupExtension(Group("-CN", "CN", [0]), doi, "[#6X2]#[#7X1H0]", None),
            GroupExtension(
                Group("-NO2", "[N+]([O-])[O]", [0], "SMARTS"),
                doi,
                "[$([#7X3,#7X3+][!#8])](=[O])~[O-]",
                None,
            ),
            GroupExtension(Group("-SH", "S", [0]), doi, "[SX2H]", None),
            GroupExtension(
                Group("-S- (non-ring)", "S", [0, 0]), doi, "[#16X2H0;!R]", None
            ),
            GroupExtension(Group("-S- (ring)", "S", [0, 0]), doi, "[#16X2H0;R]", None),
        ]


class Unifac(BasisSet):
    """UNIFAC group contribution method for estimating activity coefficients.

    This class predefines a library of chemical groups based on the UNIFAC method.
    """

    def __init__(self):
        """Initialize the UNIFAC library with predefined groups."""
        super().__init__()
        self.doi = "https://doi.org/10.1021/i200013a024"
        self.ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/UNIFAC.jl"
        self.group_traces = [
            GroupExtension(Group('CH3', '[CX4;H3;!R]', [0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH2', '[CX4;H2;!R]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH', '[CX4;H1;!R]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('C', '[CX4;H0;!R]', [0,0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH2=CH', '[CX3;H2]=[CX3;H1]', [0], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('CH=CH', '[CX3;H1]=[CX3;H1]', [0,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('CH2=C', '[CX3;H2]=[CX3;H0]', [0,0], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('CH=C', '[CX3;H1]=[CX3;H0]', [0,0,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('OH (P)', '[OH1;$([OH1][CX4H2])]', [0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('ACOH', '[cX3;H0;R][OX2;H1]', [0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('CH2CO', '[CX4;H2][CX3;!H1](=O)', [0], 'SMARTS'), self.doi, 'CC=O', None),
            GroupExtension(Group('CHO', '[CX3H1](=O)', [0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('CH3COO', '[CH3][CX3;H0](=[O])[OH0]', [3], 'SMARTS'), self.doi, 'CC(=O)O', None),
            GroupExtension(Group('CH2COO', '[CX4;H2][CX3](=[OX1])[OX2]', [0,3], 'SMARTS'), self.doi, 'CC(=O)O', None),
            GroupExtension(Group('HCOO', '[CX3;H1](=[OX1])[OX2]', [2], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('CH3O', '[CH3;!R][OH0;!R]', [0,1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CH2O', '[CH2;!R][OH0;!R]', [0,1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CHO', '[C;H1;!R][OH0;!R]', [0,1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('ACH', '[cX3;H1]', [0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('AC', '[cX3;H0]', [0,0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('ACCH3', '[cX3;H0][CX4;H3]', [0,0], 'SMARTS'), self.doi, 'CringCH3', None),
            GroupExtension(Group('ACCH2', '[cX3;H0][CX4;H2]', [0,0,1], 'SMARTS'), self.doi, 'CringCH2', None),
            GroupExtension(Group('ACCH', '[cX3;H0][CX4;H1]', [0,0,1,1], 'SMARTS'), self.doi, 'CringCH', None),
            GroupExtension(Group('CH2NH2', '[CX4;H2][NX3;H2]', [0], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CHNH2', '[CX4;H1][NX3;H2]', [0,0], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH3NH', '[CX4;H3][NX3;H1]', [1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH2NH', '[CX4;H2][NX3;H1]', [0,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CHNH', '[CX4;H1][NX3;H1]', [0,0,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH3N', '[CX4;H3][NX3;H0]', [1,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH2N', '[CX4;H2][NX3;H0]', [0,1,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('ACNH2', '[c][NX3;H2]', [0,0], 'SMARTS'), self.doi, 'C(=O)N', None),
            GroupExtension(Group('AC2H2N', '[cX3H1][n][cX3H1]', [0,1,2], 'SMARTS'), self.doi, 'C(=O)NC=O', None),
            GroupExtension(Group('AC2HN', '[cX3H0][n][cX3H1]', [0,0,1,2], 'SMARTS'), self.doi, 'C(=O)NC=O', None),
            GroupExtension(Group('AC2N', '[cX3H0][n][cX3H0]', [0,0,1,2,2], 'SMARTS'), self.doi, 'CringNringCring', None),
            GroupExtension(Group('CH2CN', '[CX4;H2][CX2]#[NX1]', [0], 'SMARTS'), self.doi, 'CC#N', None),
            GroupExtension(Group('COO', '[CX3,cX3](=[OX1])[OX2H0,oX2H0]', [0,1], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('COOH', '[CX3](=[OX1])[O;H1]', [0], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('CH2CL', '[CX4;H2;!$(C(Cl)(Cl))](Cl)', [0], 'SMARTS'), self.doi, 'C(Cl)', None),
            GroupExtension(Group('CHCL', '[CX4;H1;!$(C(Cl)(Cl))](Cl)', [0,0], 'SMARTS'), self.doi, 'C(Cl)', None),
            GroupExtension(Group('CCL', '[CX4;H0](Cl)', [0,0,0], 'SMARTS'), self.doi, 'C(Cl)', None),
            GroupExtension(Group('CHCL2', '[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', [0], 'SMARTS'), self.doi, 'C(Cl)(Cl)', None),
            GroupExtension(Group('CCL2', '[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', [0,0], 'SMARTS'), self.doi, 'C(Cl)(Cl)', None),
            GroupExtension(Group('CCL3', '[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)', [0], 'SMARTS'), self.doi, 'C(Cl)(Cl)(Cl)', None),
            GroupExtension(Group('ACCL', '[c][Cl]', [0,0], 'SMARTS'), self.doi, 'C(=O)(Cl)', None),
            GroupExtension(Group('CH2NO2', '[CX4;H2][N+X3](=[OX1])([O-X1])', [0], 'SMARTS'), self.doi, 'CN(=O)O', None),
            GroupExtension(Group('CHNO2', '[CX4;H1][N+X3](=[OX1])([O-X1])', [0,0], 'SMARTS'), self.doi, 'CN(=O)O', None),
            GroupExtension(Group('ACNO2', '[cX3][N+X3](=[OX1])([O-X1])', [0,0], 'SMARTS'), self.doi, 'C(=O)N(=O)O', None),
            GroupExtension(Group('CH2SH', '[SX2H][CX4;H2]', [0], 'SMARTS'), self.doi, 'CS', None),
            GroupExtension(Group('I', '[I]', [0], 'SMARTS'), self.doi, 'I', None),
            GroupExtension(Group('BR', '[Br]', [0], 'SMARTS'), self.doi, '(Br)', None),
            GroupExtension(Group('CH=-C', '[CX2;H1]#[CX2;H0]', [0], 'SMARTS'), self.doi, 'C#C', None),
            GroupExtension(Group('C=-C', '[CX2;H0]#[CX2;H0]', [0,1], 'SMARTS'), self.doi, 'C#C', None),
            GroupExtension(Group('CL-(C=C)', '[$([Cl;H0]([C]=[C]))]', [0], 'SMARTS'), self.doi, '(Cl)C=C', None),
            GroupExtension(Group('C=C', '[CX3;H0]=[CX3;H0]', [0,0,1,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('ACF', '[cX3][F]', [0,0], 'SMARTS'), self.doi, 'C(=O)F', None),
            GroupExtension(Group('HCON(CH2)2', '[NX3]([CX4;H2])([CX4;H2])[CX3;H1](=[OX1])', [1,1], 'SMARTS'), self.doi, 'C(=O)N(C)C', None),
            GroupExtension(Group('CF3', 'C(F)(F)F', [0], 'SMARTS'), self.doi, 'C(F)(F)F', None),
            GroupExtension(Group('CF2', 'C(F)F', [0,0], 'SMARTS'), self.doi, 'C(F)F', None),
            GroupExtension(Group('CF', 'C(F)', [0,0,0], 'SMARTS'), self.doi, 'CF', None),
            GroupExtension(Group('CY-CH2', '[CH2;R]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CY-CH', '[CH1;R]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CY-C', '[CH0;R]', [0,0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('OH (S)', '[OH1;$([OH1][CX4H1])]', [0], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('OH (T)', '[OH1;$([OH1][CX4H0])]', [0], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CY-CH2O', '[CX4H2;R][OX2;R;$(O(CC)C)][CX4H2;R][OX2;R][CX4H2;R]', [0,1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CNH2', '[CX4H0][NH2]', [0,0,0], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CONH2', '[CX3H0](=[OX1H0])[NX3H2]', [0], 'SMARTS'), self.doi, 'C(=O)N', None),
            GroupExtension(Group('CONHCH3', '[OX1H0;!R]=[CX3H0;!R][NH1X3;!R][CH3;!R]', [1], 'SMARTS'), self.doi, 'C(=O)NC', None),
            GroupExtension(Group('CONHCH2', '[CH2X4;!R][NH1X3;!R][CX3H0;!R]=[OX1H0;!R]', [0,3], 'SMARTS'), self.doi, 'C(=O)NC', None),
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
        doi = "https://doi.org/10.1080/00268978800101601"
        ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/SAFTgammaMie.jl"
        self.group_traces = [
            GroupExtension(Group("-CH3", "[CH3]", [0], "SMARTS"), doi, "[CX4H3]", None),
        ]
        print(
            "This library is not fully implemented yet... please raise a GitHub issue if you wish to use this."
        )
