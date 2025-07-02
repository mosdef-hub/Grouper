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
# Declaring namedtuple for defining adding groups and their SMARTS to a library.
GroupExtension = namedtuple(
    "GroupExtension", ["group", "doi", "extended_smarts", "priority"]
)


class BasisSet(object):
    """Base class for managing a library of chemical groups.

    Attributes:
        group_traces_ (list): A list of GroupExtension objects representing the groups in the library.
        group_traces_ (list): A list of GroupExtension objects representing the groups in the library.
    """

    def __init__(self):
        """Initialize an empty BasisSet."""
        self.group_traces_ = []
        self.group_traces_ = []

    def __repr__(self):
        """Return a string representation of the BasisSet."""
        init_message = (
            f"{self.__class__.__name__} has {len(self.group_traces_)} unique groups."
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
            group (Group): The chemical group to add.
            doi (str, optional): doi reference for the group. Defaults to None.
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
        self.doi = "https://doi.org/10.1021/i200013a024"
        ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/Joback.jl"
        self.group_traces = [
            GroupExtension(Group('-CH3', '[CX4H3]', [0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('-CH2-', '[!R;CX4H2]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('>CH-', '[!R;CX4H]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('>C<', '[!R;CX4H0]', [0,0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH2=CH-', '[CX3H2][CX3H1]', [1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('-CH=CH-', '[CX3H1][CX3H1]', [0,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('=C<', '[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('=C=', '[$([CX2H0](=*)=*)]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH', '[$([CX2H1]#[!#7])]', [0], 'SMARTS'), self.doi, 'C#', None),
            GroupExtension(Group('C', '[$([CX2H0]#[!#7])]', [0,0], 'SMARTS'), self.doi, 'C#', None),
            GroupExtension(Group('ring-CH2-', '[R;CX4H2]', [0,0], 'SMARTS'), self.doi, 'Cring', None),
            GroupExtension(Group('ring>CH-', '[R;CX4H]', [0,0,0], 'SMARTS'), self.doi, 'Cring', None),
            GroupExtension(Group('ring>C<', '[R;CX4H0]', [0,0,0,0], 'SMARTS'), self.doi, 'Cring', None),
            GroupExtension(Group('ring=CH-', '[R;CX3H1,cX3H1]', [0,0], 'SMARTS'), self.doi, 'Cring', None),
            GroupExtension(Group('ring=C<', '[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]', [0,0,0], 'SMARTS'), self.doi, 'Cring', None),
            GroupExtension(Group('-F', '[F;X1]', [0], 'SMARTS'), self.doi, 'F', None),
            GroupExtension(Group('-Cl', '[Cl;X1]', [0], 'SMARTS'), self.doi, 'Cl', None),
            GroupExtension(Group('-Br', '[Br;X1]', [0], 'SMARTS'), self.doi, 'Br', None),
            GroupExtension(Group('-I', '[I;X1]', [0], 'SMARTS'), self.doi, 'I', None),
            GroupExtension(Group('-OH (alcohol)', '[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]', [0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('-OH (phenol)', '[O;H1;$(O-!@c)]', [0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('-O- (non-ring)', '[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]', [0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('-O- (ring)', '[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]', [0,0], 'SMARTS'), self.doi, 'Oring', None),
            GroupExtension(Group('>C=O (non-ring)', '[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O', [0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('>C=O (ring)', '[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O', [0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('O=CH- (aldehyde)', '[CH;D2](=O)', [0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('-COOH (acid)', '[OX2H]-[C]=O', [1], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('-COO- (ester)', '[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]', [0,2], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('=O (other than above)', '[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]', [0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('-NH2', '[NX3H2]', [0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('>NH (non-ring)', '[NX3H1;!R]', [0,0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('>NH (ring)', '[#7X3H1;R]', [0,0], 'SMARTS'), self.doi, 'Nring', None),
            GroupExtension(Group('>N- (non-ring)', '[#7X3H0;!$([#7](~O)~O)]', [0,0,0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('-N= (non-ring)', '[#7X2H0;!R]', [0,0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('-N= (ring)', '[#7X2H0;R]', [0,0], 'SMARTS'), self.doi, 'Nring', None),
            GroupExtension(Group('=NH', '[#7X2H1]', [0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('-CN', '[#6X2]#[#7X1H0]', [0], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('-NO2', '[$([#7X3,#7X3+][!#8])](=[O])[O-]', [0], 'SMARTS'), self.doi, 'N(=O)O', None),
            GroupExtension(Group('-SH', '[SX2H]', [0], 'SMARTS'), self.doi, 'S', None),
            GroupExtension(Group('-S- (non-ring)', '[#16X2H0;!R]', [0,0], 'SMARTS'), self.doi, 'S', None),
            GroupExtension(Group('-S- (ring)', '[#16X2H0;R]', [0,0], 'SMARTS'), self.doi, 'Sring', None),
        ]


class Unifac(BasisSet):
    """UNIFAC group contribution method for estimating activity coefficients.

    This class predefines a library of chemical groups based on the UNIFAC method.
    """

    def __init__(self):
        """Initialize the UNIFAC library with predefined groups."""
        super().__init__()
        self.doi = "https://self.doi.org/10.1021/i200013a024"
        self.ref = "https://github.com/ClapeyronThermo/GCIdentifier.jl/tree/main/src/database/UNIFAC.jl"
        self.group_traces = [
            GroupExtension(Group('CH3', '[CX4;H3;!R]', [0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH2', '[CX4;H2;!R]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH', '[CX4;H1;!R]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('C', '[CX4;H0;!R]', [0,0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH=CH2', '[CX3;H2]=[CX3;H1]', [1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('CH=CH', '[CX3;H1]=[CX3;H1]', [0,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('C=CH2', '[CX3;H2]=[CX3;H0]', [1,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('C=CH', '[CX3;H1]=[CX3;H0]', [0,1,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('OH (P)', '[OH1;$([OH1][CX4H2])]', [0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('CH3OH', '[CX4;H3][OX2;H1]', [], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('H2O', '[OH2]', [], 'SMARTS'), self.doi, 'OH2', None),
            GroupExtension(Group('ACOH', '[cX3;H0;R][OX2;H1]', [0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('CH3CO', '[CX4;H3][CX3;!H1](=O)', [1], 'SMARTS'), self.doi, 'CC=O', None),
            GroupExtension(Group('CH2CO', '[CX4;H2][CX3;!H1](=O)', [0,1], 'SMARTS'), self.doi, 'CC=O', None),
            GroupExtension(Group('CHO', '[CX3H1](=O)', [0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('CH3COO', '[CH3][CX3;H0](=[O])[OH0]', [3], 'SMARTS'), self.doi, 'CC(=O)O', None),
            GroupExtension(Group('CH2COO', '[CX4;H2][CX3](=[OX1])[OX2]', [0,3], 'SMARTS'), self.doi, 'CC(=O)O', None),
            GroupExtension(Group('HCOO', '[CX3;H1](=[OX1])[OX2]', [2], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('CH3O', '[CH3;!R][OH0;!R]', [1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CH2O', '[CH2;!R][OH0;!R]', [0,1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CHO', '[C;H1;!R][OH0;!R]', [0,0,1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('ACH', '[cX3;H1]', [0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('AC', '[cX3;H0]', [0,0,0], 'SMARTS'), self.doi, 'C=O', None),
            GroupExtension(Group('ACCH3', '[cX3;H0][CX4;H3]', [0,0], 'SMARTS'), self.doi, 'CringCH3', None),
            GroupExtension(Group('ACCH2', '[cX3;H0][CX4;H2]', [0,0,1], 'SMARTS'), self.doi, 'CringCH2', None),
            GroupExtension(Group('ACCH', '[cX3;H0][CX4;H1]', [0,0,1,1], 'SMARTS'), self.doi, 'CringCH', None),
            GroupExtension(Group('THF', '[CX4;H2;R;$(C(C)OCC)][OX2;R][CX4;H2;R]', [], 'SMARTS'), self.doi, 'C1CCOC1', None),
            GroupExtension(Group('CH3NH2', '[CX4;H3][NX3;H2]', [], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH2NH2', '[CX4;H2][NX3;H2]', [0], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CHNH2', '[CX4;H1][NX3;H2]', [0,0], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH3NH', '[CX4;H3][NX3;H1]', [1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH2NH', '[CX4;H2][NX3;H1]', [0,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CHNH', '[CX4;H1][NX3;H1]', [0,0,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH3N', '[CX4;H3][NX3;H0]', [1,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('CH2N', '[CX4;H2][NX3;H0]', [0,1,1], 'SMARTS'), self.doi, 'CN', None),
            GroupExtension(Group('ACNH2', '[c][NX3;H2]', [0,0], 'SMARTS'), self.doi, 'CringNH2', None),
            GroupExtension(Group('AC2H2N', '[cX3H1][n][cX3H1]', [0,1,2], 'SMARTS'), self.doi, 'CringNringCring', None),
            GroupExtension(Group('AC2HN', '[cX3H0][n][cX3H1]', [0,0,1,2], 'SMARTS'), self.doi, 'CringNringCring', None),
            GroupExtension(Group('AC2N', '[cX3H0][n][cX3H0]', [0,0,1,2,2], 'SMARTS'), self.doi, 'CringNringCring', None),
            GroupExtension(Group('CH3CN', '[CX4;H3][CX2]#[NX1]', [], 'SMARTS'), self.doi, 'CC#N', None),
            GroupExtension(Group('CH2CN', '[CX4;H2][CX2]#[NX1]', [0], 'SMARTS'), self.doi, 'CC#N', None),
            GroupExtension(Group('C(=O)O', '[CX3,cX3](=[OX1])[OX2H0,oX2H0]', [0,2], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('COOH', '[CX3](=[OX1])[O;H1]', [0], 'SMARTS'), self.doi, 'C(=O)OH', None),
            GroupExtension(Group('HCOOH', '[CX3;H1](=[OX1])[OX2;H1]', [], 'SMARTS'), self.doi, 'HC(=O)OH', None),
            GroupExtension(Group('CH2CL', '[CX4;H2;!$(C(Cl)(Cl))](Cl)', [0], 'SMARTS'), self.doi, 'C(Cl)', None),
            GroupExtension(Group('CHCL', '[CX4;H1;!$(C(Cl)(Cl))](Cl)', [0,0], 'SMARTS'), self.doi, 'C(Cl)', None),
            GroupExtension(Group('CCL', '[CX4;H0](Cl)', [0,0,0], 'SMARTS'), self.doi, 'C(Cl)', None),
            GroupExtension(Group('CH2CL2', '[CX4;H2;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', [], 'SMARTS'), self.doi, 'C(Cl)(Cl)', None),
            GroupExtension(Group('CHCL2', '[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', [0], 'SMARTS'), self.doi, 'C(Cl)(Cl)', None),
            GroupExtension(Group('CCL2', '[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', [0,0], 'SMARTS'), self.doi, 'C(Cl)(Cl)', None),
            GroupExtension(Group('CHCL3', '[CX4;H1;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)', [], 'SMARTS'), self.doi, 'C(Cl)(Cl)(Cl)', None),
            GroupExtension(Group('CCL3', '[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)', [0], 'SMARTS'), self.doi, 'C(Cl)(Cl)(Cl)', None),
            GroupExtension(Group('CCL4', '[CX4;H0]([Cl])([Cl])([Cl])([Cl])', [], 'SMARTS'), self.doi, 'C(Cl)(Cl)(Cl)(Cl', None),
            GroupExtension(Group('ACCL', '[c][Cl]', [0,0], 'SMARTS'), self.doi, 'C(Cl)', None),
            GroupExtension(Group('CH3NO2', '[CX4;H3][N+X3](=[OX1])([O-X1])', [], 'SMARTS'), self.doi, 'CN(=O)O', None),
            GroupExtension(Group('CH2NO2', '[CX4;H2][N+X3](=[OX1])([O-X1])', [0], 'SMARTS'), self.doi, 'CN(=O)O', None),
            GroupExtension(Group('CHNO2', '[CX4;H1][N+X3](=[OX1])([O-X1])', [0,0], 'SMARTS'), self.doi, 'CN(=O)O', None),
            GroupExtension(Group('ACNO2', '[cX3][N+X3](=[OX1])([O-X1])', [0,0], 'SMARTS'), self.doi, 'C(=O)N(=O)O', None),
            GroupExtension(Group('CS2', 'C(=S)=S', [], 'SMARTS'), self.doi, 'S=C=S', None),
            GroupExtension(Group('CH3SH', '[SX2H][CX4;H3]', [], 'SMARTS'), self.doi, 'CH3SH', None),
            GroupExtension(Group('CH2SH', '[SX2H][CX4;H2]', [0], 'SMARTS'), self.doi, 'CH2SH', None),
            GroupExtension(Group('FURFURAL', 'c1cc(oc1)C=O', [], 'SMARTS'), self.doi, 'c1cc(oc1)C=O', None),
            GroupExtension(Group('DOH', '[OX2;H1][CX4;H2][CX4;H2][OX2;H1]', [], 'SMARTS'), self.doi, 'OHCH2CH2OH', None),
            GroupExtension(Group('I', '[I]', [0], 'SMARTS'), self.doi, 'I', None),
            GroupExtension(Group('BR', '[Br]', [0], 'SMARTS'), self.doi, '(Br)', None),
            GroupExtension(Group('CH=-C', '[CX2;H1]#[CX2;H0]', [1], 'SMARTS'), self.doi, 'C#C', None),
            GroupExtension(Group('C=-C', '[CX2;H0]#[CX2;H0]', [0,1], 'SMARTS'), self.doi, 'C#C', None),
            GroupExtension(Group('DMSO', '[SX3H0](=[OX1])([CX4;H3])[CX4;H3]', [], 'SMARTS'), self.doi, 'CS(=O)C', None),
            GroupExtension(Group('ACRY', '[CX3;H2]=[CX3;H1][CX2;H0]#[NX1;H0]', [], 'SMARTS'), self.doi, 'C=CC#N', None),
            GroupExtension(Group('CL-(C=C)', '[$([Cl;H0]([C]=[C]))]', [0], 'SMARTS'), self.doi, 'Cl', None),
            GroupExtension(Group('C=C', '[CX3;H0]=[CX3;H0]', [0,0,1,1], 'SMARTS'), self.doi, 'C=C', None),
            GroupExtension(Group('ACF', '[cX3][F]', [0,0], 'SMARTS'), self.doi, 'CringF', None),
            GroupExtension(Group('DMF', '[CX4;H3][N]([CX4;H3])[CX3;H1]=[O]', [], 'SMARTS'), self.doi, 'CN(C)C=O', None),
            GroupExtension(Group('HCON(CH2)2', '[NX3]([CX4;H2])([CX4;H2])[CX3;H1](=[OX1])', [1,2], 'SMARTS'), self.doi, 'C(=O)N(C)C', None),
            GroupExtension(Group('CF3', 'C(F)(F)F', [0], 'SMARTS'), self.doi, 'C(F)(F)F', None),
            GroupExtension(Group('CF2', 'C(F)F', [0,0], 'SMARTS'), self.doi, 'C(F)F', None),
            GroupExtension(Group('CF', 'C(F)', [0,0,0], 'SMARTS'), self.doi, 'CF', None),
            GroupExtension(Group('CY-CH2', '[CH2;R]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CY-CH', '[CH1;R]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CY-C', '[CH0;R]', [0,0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('OH (S)', '[OH1;$([OH1][CX4H1])]', [0], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('OH (T)', '[OH1;$([OH1][CX4H0])]', [0], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CY-CH2O', '[CX4H2;R][OX2;R;$(O(CC)C)][CX4H2;R][OX2;R][CX4H2;R]', [0,4], 'SMARTS'), self.doi, 'COCOC', None),
            GroupExtension(Group('TRIOXAN', '[CX4H2;R][OX2;R;$(O(C)C)]', [0,1], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CNH2', '[CX4H0][NH2]', [0,0,0], 'SMARTS'), self.doi, 'CNH2', None),
            GroupExtension(Group('NMP', '[OX1H0]=[C;R][NX3H0;R][CH3]', [1,2], 'SMARTS'), self.doi, 'O=CNCH3', None),
            GroupExtension(Group('NEP', '[OX1H0]=[CH0X3;R][NX3H0;R][CH2]', [1,2,3], 'SMARTS'), self.doi, 'NEP', None),
            GroupExtension(Group('NIPP', '[OX1H0;!R]=[CX3H0;R][NX3H0;R][C;!R]', [1,2,3,3,3], 'SMARTS'), self.doi, 'O=CNC', None),
            GroupExtension(Group('NTBP', '[OX1H0;!R]=[CH0X3;R][NX3H0;R][CH0;!R]', [0,1,2,2,2], 'SMARTS'), self.doi, 'O=CNC', None),
            GroupExtension(Group('CONH2', '[CX3](=[OX1H0])[NX3H2]', [0], 'SMARTS'), self.doi, 'C(=O)N', None),
            GroupExtension(Group('CONHCH3', '[OX1H0;!R]=[CX3;!R][NH1X3;!R][CH3;!R]', [1], 'SMARTS'), self.doi, 'C(=O)NC', None),
            GroupExtension(Group('CONHCH2', '[CH2X4;!R][NH1X3;!R][CX3;!R]=[OX1H0;!R]', [0,2], 'SMARTS'), self.doi, 'CNC(=O)', None),
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
        self.group_traces = [
            GroupExtension(Group('CH3', '[CX4H3]', [0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH2', '[!R;CX4H2]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH', '[!R;CX4H]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('C', '[!R;CX4H0]', [0,0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('aCH', '[cX3;H1]', [0], 'SMARTS'), self.doi, 'C(=O)', None),
            GroupExtension(Group('aCCH2', '[cX3;H0][CX4;H2]', [0,0,1], 'SMARTS'), self.doi, 'CC(=O)', None),
            GroupExtension(Group('aCCH', '[cX3;H0][CX4;H1]', [0,0,1,1], 'SMARTS'), self.doi, 'CC(=O)', None),
            GroupExtension(Group('CH2=', '[CX3H2]', [0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CH=', '[!R;CX3H1;!$([CX3H1](=O))]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('cCH2', '[CH2;R]', [0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('COOH', '[OX2H]-[C]=O', [1], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('COO', '[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]', [0,2], 'SMARTS'), self.doi, 'C(=O)O', None),
            GroupExtension(Group('OH', '[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]', [0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('CH2OH', '[CX4;H2;!R][OH1]', [0], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CHOH', '[CX4;H1;!R][OH1]', [0,0], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('NH2', '[NX3H2]', [0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('NH', '[NX3H1;!R]', [0,0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('N', '[#7X3H0;!$([#7](~O)~O)]', [0,0,0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('cNH', '[#7X3H1;R]', [0,0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('cN', '[#7X3H0;R]', [0,0,0], 'SMARTS'), self.doi, 'N', None),
            GroupExtension(Group('C=', '[!R;CX3H0;!$([CX3H0](=O))]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('aCCH3', '[cX3;H0][CX4;H3]', [0,0], 'SMARTS'), self.doi, 'CringC', None),
            GroupExtension(Group('aCOH', '[cX3;H0;R][OX2;H1]', [0,0], 'SMARTS'), self.doi, 'CringO', None),
            GroupExtension(Group('cCH', '[CH1;R]', [0,0,0], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('cCHNH', '[CH1;R][NH1;!R]', [0,0,1], 'SMARTS'), self.doi, 'CringN', None),
            GroupExtension(Group('cCHN', '[CH1;R][NH0;!R]', [0,0,1,1], 'SMARTS'), self.doi, 'CringN', None),
            GroupExtension(Group('aCCOaC', '[cH0][C;!R](=O)[cH0]', [0,0,0,3,3,3], 'SMARTS'), self.doi, 'CC(=O)C', None),
            GroupExtension(Group('aCCOOH', '[OX2H]-[C](=O)[cH0]', [3,3,3], 'SMARTS'), self.doi, 'OC(=O)C', None),
            GroupExtension(Group('aCNHaC', '[cH0][NH1;!R][cH0]', [0,0,0,2,2,2], 'SMARTS'), self.doi, 'CNC', None),
            GroupExtension(Group('CH3CO', '[CH3][CX3](=O)', [1], 'SMARTS'), self.doi, 'CC=O', None),
            GroupExtension(Group('eO', '[OH0;!R;$([OH0;!R][CH3;!R]);$([OH0;!R][CH2;!R])]', [0,0], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('cO', '[OH0;!R;$([OH0;!R][CH2;!R])]', [0,0], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('CH3OH', '[CX4;H3][OX2;H1]', [], 'SMARTS'), self.doi, 'CO', None),
            GroupExtension(Group('H2O', '[OH2]', [], 'SMARTS'), self.doi, 'O', None),
            GroupExtension(Group('CH4', '[CH4]', [], 'SMARTS'), self.doi, 'C', None),
            GroupExtension(Group('CO2', '[C;X2](=O)(=O)', [], 'SMARTS'), self.doi, 'C(=O)=O', None),
        ]

        print(
            "This library is not fully implemented yet... please raise a GitHub issue if you wish to use this."
        )
