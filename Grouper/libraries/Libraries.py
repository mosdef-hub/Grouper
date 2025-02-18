"""Module for pre-loading libraries of GroupGraphs."""

# from genGrouper.library import Library, Basis, Joback
from collections import namedtuple
from typing import Dict, Union

from Grouper import Group

# Declaring namedtuple for defining adding nodes and their SMARTS to a library.
NodeTrace = namedtuple("NodeTrace", ["node", "doi", "smarts", "priority"])


class BasisSet(object):
    def __init__(self):
        # trace = (Group(name, smiles, hubs), doi,  smarts, priorty) # could also name moiety
        self.node_traces_ = set()

    def __repr__(self):
        return f"{self.__class__.__name__} has {len(self.node_traces_)} unique nodes."

    def add_node(
        self, node: Group, doi: str = None, smarts: str = None, priority: int = None
    ):
        trace = NodeTrace(node, doi, smarts, priority)
        self.node_traces_.add(trace)

    def query_nodes(self, query: Dict[str, Union[str, int]]):
        """Query the library for nodes that match the query."""
        matched_nodes = []
        for node in self.node_traces:
            if any([getattr(node, var) == value for var, value in query.items()]):
                matched_nodes.append(node)
        return matched_nodes

    def get_nodes(self):
        for trace in self.node_traces:
            yield trace.node

    def visualize_library(self):
        "Do rdkit visualization"
        # rdkit.Chem.Draw.MolToBlockImage(list_of_mols, n_mols_per_row)
        pass

    @property
    def node_traces(self):
        return self.node_traces_

    @node_traces.setter
    def node_traces(self, value):
        self.node_traces_ = value

    @property
    def n_nodes(self):
        return len(self.node_traces_)


class Joback(BasisSet):
    """Joback group contribution method for estimating physical properties of organic compounds."""

    def __init__(self):
        super().__init__()
        self.node_traces = {
            NodeTrace(Group("-CH3", "CH3", [0]), "", "CX4H3", None),
            NodeTrace(Group("-CH2-", "CH2", [0, 0]), "", "[!R;CX4H2]", None),
            NodeTrace(Group(">CH-", "CH", [0, 0, 0]), "", "[!R;CX4H]", None),
            NodeTrace(Group(">C<", "C", [0, 0, 0, 0]), "", "[!R;CX4H0]", None),
            NodeTrace(Group("CH2=CH-", "C=C", [0]), "", "[CX3H2][CX3H1]", None),
            NodeTrace(Group("-CH=CH-", "C=C", [0, 1]), "", "[CX3H1][CX3H1]", None),
            NodeTrace(
                Group("=C<", "C", [0, 0, 0]),
                "",
                "[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]",
                None,
            ),
            NodeTrace(Group("=C=", "C", [0, 0]), "", "[$([CX2H0](=*)=*)]", None),
            NodeTrace(Group("CH", "C", [0]), "", "[$([CX2H1]#[!#7])]", None),
            NodeTrace(Group("C", "C", [0, 0]), "", "[$([CX2H0]#[!#7])]", None),
            NodeTrace(Group("ring-CH2-", "C", [0, 0]), "", "[R;CX4H2]", None),
            NodeTrace(Group("ring>CH-", "C", [0, 0, 0]), "", "[R;CX4H]", None),
            NodeTrace(Group("ring>C<", "C", [0, 0, 0, 0]), "", "[R;CX4H0]", None),
            NodeTrace(Group("ring=CH-", "C", [0, 0]), "", "[R;CX3H1,cX3H1]", None),
            NodeTrace(
                Group("ring=C<", "C", [0, 0, 0]),
                "",
                "[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]",
                None,
            ),
            NodeTrace(Group("-F", "F", [0]), "", "[F;X1]", None),
            NodeTrace(Group("-Cl", "Cl", [0]), "", "[Cl;X1]", None),
            NodeTrace(Group("-Br", "Br", [0]), "", "[Br;X1]", None),
            NodeTrace(Group("-I", "I", [0]), "", "[I;X1]", None),
            NodeTrace(
                Group("-OH (alcohol)", "O", [0]),
                "",
                "[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]",
                None,
            ),
            NodeTrace(Group("-OH (phenol)", "O", [0]), "", "[O;H1;$(O-!@c)]", None),
            NodeTrace(
                Group("-O- (non-ring)", "O", [0, 0]),
                "",
                "[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]",
                None,
            ),
            NodeTrace(
                Group("-O- (ring)", "O", [0, 0]),
                "",
                "[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]",
                None,
            ),
            NodeTrace(
                Group(">C=O (non-ring)", "C=O", [0, 0]),
                "",
                "[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O",
                None,
            ),
            NodeTrace(
                Group(">C=O (ring)", "C=O", [0, 0]),
                "",
                "[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O",
                None,
            ),
            NodeTrace(Group("O=CH- (aldehyde)", "O=C", [1]), "", "[CH;D2](=O)", None),
            NodeTrace(Group("-COOH (acid)", "C(=O)O", [0]), "", "[OX2H]-[C]=O", None),
            NodeTrace(
                Group("-COO- (ester)", "C(=O)O", [0, 2]),
                "",
                "[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]",
                None,
            ),
            NodeTrace(
                Group("=O (other than above)", "O", [0]),
                "",
                "[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]",
                None,
            ),
            NodeTrace(Group("-NH2", "N", [0]), "", "[NX3H2]", None),
            NodeTrace(Group(">NH (non-ring)", "N", [0, 0]), "", "[NX3H1;!R]", None),
            NodeTrace(Group(">NH (ring)", "N", [0, 0]), "", "[#7X3H1;R]", None),
            NodeTrace(
                Group(">N- (non-ring)", "N", [0, 0, 0]),
                "",
                "[#7X3H0;!$([#7](~O)~O)]",
                None,
            ),
            NodeTrace(Group("-N= (non-ring)", "N", [0, 0]), "", "[#7X2H0;!R]", None),
            NodeTrace(Group("-N= (ring)", "N", [0, 0]), "", "[#7X2H0;R]", None),
            NodeTrace(Group("=NH", "N", [0]), "", "[#7X2H1]", None),
            NodeTrace(Group("-CN", "CN", [0]), "", "[#6X2]#[#7X1H0]", None),
            NodeTrace(
                Group("-NO2", "NO2", [0]),
                "",
                "[$([#7X3,#7X3+][!#8])](=[O])~[O-]",
                None,
            ),
            NodeTrace(Group("-SH", "S", [0]), "", "[SX2H]", None),
            NodeTrace(Group("-S- (non-ring)", "S", [0, 0]), "", "[#16X2H0;!R]", None),
            NodeTrace(Group("-S- (ring)", "S", [0, 0]), "", "[#16X2H0;R]", None),
        }


class Unifac(BasisSet):
    def __init__(self):
        super().__init__()
        self.node_traces = {
            NodeTrace(Group('CH3', 'C', [0]), '', '[CX4;H3;!R]', None),
            NodeTrace(Group('CH2', 'C', [0,0]), '', '[CX4;H2;!R]', None),
            NodeTrace(Group('CH', 'C', [0,0,0]), '', '[CX4;H1;!R]', None),
            NodeTrace(Group('C', 'C', [0,0,0,0]), '', '[CX4;H0;!R]', None),
            NodeTrace(Group('CH2=CH', 'C=C', [0]), '', '[CX3;H2]=[CX3;H1]', None),
            NodeTrace(Group('CH=CH', 'C=C', [0,1]), '', '[CX3;H1]=[CX3;H1]', None),
            NodeTrace(Group('CH2=C', 'C=C', [0,0]), '', '[CX3;H2]=[CX3;H0]', None),
            NodeTrace(Group('CH=C', 'C=C', [0,0,1]), '', '[CX3;H1]=[CX3;H0]', None),
            NodeTrace(Group('OH (P)', 'O', [0]), '', '[OH1;$([OH1][CX4H2])]', None),
            NodeTrace(Group('CH3OH', 'CO', [0]), '', '[CX4;H3][OX2;H1]', None),
            NodeTrace(Group('H2O', 'O', []), '', '[OH2]', None),
            NodeTrace(Group('ACOH', 'C=O', [0,0]), '', '[cX3;H0;R][OX2;H1]', None),
            NodeTrace(Group('CH3CO', 'CC=O', []), '', '[CX4;H3][CX3;!H1](=O)', None),
            NodeTrace(Group('CH2CO', 'CC=O', [0]), '', '[CX4;H2][CX3;!H1](=O)', None),
            NodeTrace(Group('CHO', 'C=O', [0,0]), '', '[CX3H1](=O)', None),
            NodeTrace(Group('CH3COO', 'CC(=O)O', [3]), '', '[CH3][CX3;H0](=[O])[OH0]', None),
            NodeTrace(Group('CH2COO', 'CC(=O)O', [0,3]), '', '[CX4;H2][CX3](=[OX1])[OX2]', None),
            NodeTrace(Group('HCOO', 'C(=O)O', [0,3]), '', '[CX3;H1](=[OX1])[OX2]', None),
            NodeTrace(Group('CH3O', 'CO', [0,1]), '', '[CH3;!R][OH0;!R]', None),
            NodeTrace(Group('CH2O', 'CO', [0,1]), '', '[CH2;!R][OH0;!R]', None),
            NodeTrace(Group('CHO', 'CO', [0,1]), '', '[C;H1;!R][OH0;!R]', None),
            NodeTrace(Group('ACH', 'C=O', [0]), '', '[cX3;H1]', None),
            NodeTrace(Group('AC', 'C=O', [0]), '', '[cX3;H0]', None),
            NodeTrace(Group('ACCH3', 'C(=O)C', [0]), '', '[cX3;H0][CX4;H3]', None),
            NodeTrace(Group('ACCH2', 'C(=O)C', [0,2]), '', '[cX3;H0][CX4;H2]', None),
            NodeTrace(Group('ACCH', 'C(=O)C', [0,2,2]), '', '[cX3;H0][CX4;H1]', None),
            NodeTrace(Group('THF', 'C1CCOC1', []), '', '[CX4;H2;R;$(C(C)OCC)][OX2;R][CX4;H2;R]', None),
            NodeTrace(Group('CH3NH2', 'CN', []), '', '[CX4;H3][NX3;H2]', None),
            NodeTrace(Group('CH2NH2', 'CN', [0]), '', '[CX4;H2][NX3;H2]', None),
            NodeTrace(Group('CHNH2', 'CN', [0,0]), '', '[CX4;H1][NX3;H2]', None),
            NodeTrace(Group('CH3NH', 'CN', [1]), '', '[CX4;H3][NX3;H1]', None),
            NodeTrace(Group('CH2NH', 'CN', [0,1]), '', '[CX4;H2][NX3;H1]', None),
            NodeTrace(Group('CHNH', 'CN', [0,0,1]), '', '[CX4;H1][NX3;H1]', None),
            NodeTrace(Group('CH3N', 'CN', [1,1]), '', '[CX4;H3][NX3;H0]', None),
            NodeTrace(Group('CH2N', 'CN', [0,1,1]), '', '[CX4;H2][NX3;H0]', None),
            NodeTrace(Group('ACNH2', 'C(=O)N', [0]), '', '[c][NX3;H2]', None),
            NodeTrace(Group('AC2H2N', 'C(=O)NC=O', [2]), '', '[cX3H1][n][cX3H1]', None),
            NodeTrace(Group('AC2HN', 'C(=O)NC=O', [0,2]), '', '[cX3H0][n][cX3H1]', None),
            NodeTrace(Group('AC2N', 'C(=O)NC=O', [0,2,3]), '', '[cX3H0][n][cX3H0]', None),
            NodeTrace(Group('CH3CN', 'CC#N', []), '', '[CX4;H3][CX2]#[NX1]', None),
            NodeTrace(Group('CH2CN', 'CC#N', [0]), '', '[CX4;H2][CX2]#[NX1]', None),
            NodeTrace(Group('COO', 'C(=O)O', [0,1]), '', '[CX3,cX3](=[OX1])[OX2H0,oX2H0]', None),
            NodeTrace(Group('COOH', 'C(=O)O', [0]), '', '[CX3](=[OX1])[O;H1]', None),
            NodeTrace(Group('HCOOH', 'C(=O)O', []), '', '[CX3;H1](=[OX1])[OX2;H1]', None),
            NodeTrace(Group('CH2CL', 'C(Cl)', [0]), '', '[CX4;H2;!$(C(Cl)(Cl))](Cl)', None),
            NodeTrace(Group('CHCL', 'C(Cl)', [0,0]), '', '[CX4;H1;!$(C(Cl)(Cl))](Cl)', None),
            NodeTrace(Group('CCL', 'C(Cl)', [0,0,0]), '', '[CX4;H0](Cl)', None),
            NodeTrace(Group('CH2CL2', 'C(Cl)(Cl)', []), '', '[CX4;H2;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', None),
            NodeTrace(Group('CH2CL', 'C(Cl)', [0]), '', '[CX4;H2;!$(C(Cl)(Cl))](Cl)', None),
            NodeTrace(Group('CHCL', 'C(Cl)', [0,0]), '', '[CX4;H1;!$(C(Cl)(Cl))](Cl)', None),
            NodeTrace(Group('CCL', 'C(Cl)', [0,0,0]), '', '[CX4;H0](Cl)', None),
            NodeTrace(Group('CH2CL2', 'C(Cl)(Cl)', []), '', '[CX4;H2;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', None),
            NodeTrace(Group('CHCL2', 'C(Cl)(Cl)', [0]), '', '[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', None),
            NodeTrace(Group('CCL2', 'C(Cl)(Cl)', [0,0]), '', '[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)', None),
            NodeTrace(Group('CHCL3', 'C(Cl)(Cl)(Cl)', []), '', '[CX4;H1;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)', None),
            NodeTrace(Group('CCL3', 'C(Cl)(Cl)(Cl)', [0]), '', '[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)', None),
            NodeTrace(Group('CCL4', 'C(Cl)(Cl)(Cl)(Cl', []), '', '[CX4;H0]([Cl])([Cl])([Cl])([Cl])', None),
            NodeTrace(Group('ACCL', 'C(=O)(Cl)', [0]), '', '[c][Cl]', None),
            NodeTrace(Group('CH3NO2', 'CN(=O)O', []), '', '[CX4;H3][NX3](=[OX1])([OX1])', None),
            NodeTrace(Group('CH2NO2', 'CN(=O)O', [0]), '', '[CX4;H2][NX3](=[OX1])([OX1])', None),
            NodeTrace(Group('CHNO2', 'CN(=O)O', [0,0]), '', '[CX4;H1][NX3](=[OX1])([OX1])', None),
            NodeTrace(Group('ACNO2', 'C(=O)N(=O)O', [0]), '', '[cX3][NX3](=[OX1])([OX1])', None),
            NodeTrace(Group('CS2', 'S=C=S', []), '', 'C(=S)=S', None),
            NodeTrace(Group('CH3SH', 'CS', []), '', '[SX2H][CX4;H3]', None),
            NodeTrace(Group('CH2SH', 'CS', [0]), '', '[SX2H][CX4;H2]', None),
            NodeTrace(Group('FURFURAL', 'c1cc(oc1)C=O', []), '', 'c1cc(oc1)C=O', None),
            NodeTrace(Group('DOH', 'OCCO', []), '', '[OX2;H1][CX4;H2][CX4;H2][OX2;H1]', None),
            NodeTrace(Group('I', 'I', [0]), '', '[I]', None),
            NodeTrace(Group('BR', '(Br)', [0]), '', '[Br]', None),
            NodeTrace(Group('CH=-C', 'C#C', [0]), '', '[CX2;H1]#[CX2;H0]', None),
            NodeTrace(Group('C=-C', 'C#C', [0,1]), '', '[CX2;H0]#[CX2;H0]', None),
            NodeTrace(Group('DMSO', 'CS(=O)C', []), '', '[SX3H0](=[OX1])([CX4;H3])[CX4;H3]', None),
            NodeTrace(Group('ACRY', 'C=CC#N', []), '', '[CX3;H2]=[CX3;H1][CX2;H0]#[NX1;H0]', None),
            NodeTrace(Group('CL-(C=C)', '(Cl)C=C', [1,2,2]), '', '[$([Cl;H0]([C]=[C]))]', None),
            NodeTrace(Group('C=C', 'C=C', [0,0,1,1]), '', '[CX3;H0]=[CX3;H0]', None),
            NodeTrace(Group('ACF', 'C(=O)F', [0]), '', '[cX3][F]', None),
            NodeTrace(Group('DMF', 'CN(C)C=O', []), '', '[CX4;H3][N]([CX4;H3])[CX3;H1]=[O]', None),
            NodeTrace(Group('HCON(CH2)2', 'C(=O)N(C)C', []), '', '[NX3]([CX4;H2])([CX4;H2])[CX3;H1](=[OX1])', None),
            NodeTrace(Group('CF3', 'C(F)(F)F', [0]), '', 'C(F)(F)F', None),
            NodeTrace(Group('CF2', 'C(F)F', [0,0]), '', 'C(F)F', None),
            NodeTrace(Group('CF', 'CF', [0,0,0]), '', 'C(F)', None),
            NodeTrace(Group('CY-CH2', 'C', [0,0]), '', '[CH2;R]', None),
            NodeTrace(Group('CY-CH', 'C', [0,0,0]), '', '[CH1;R]', None),
            NodeTrace(Group('CY-C', 'C', [0,0,0,0]), '', '[CH0;R]', None),
            NodeTrace(Group('OH (S)', 'CO', [0,0]), '', '[OH1;$([OH1][CX4H1])]', None),
            NodeTrace(Group('OH (T)', 'CO', [0,0,0]), '', '[OH1;$([OH1][CX4H0])]', None),
            NodeTrace(Group('CY-CH2O', 'CO', [0,1]), '', '[CX4H2;R][OX2;R;$(O(CC)C)][CX4H2;R][OX2;R][CX4H2;R]', None),
            NodeTrace(Group('TRIOXAN', 'C1OCOCO1', []), '', '[CX4H2;R][OX2;R;$(O(C)C)]', None),
            NodeTrace(Group('CNH2', 'CN', [0,0,0]), '', '[CX4H0][NH2]', None),
            NodeTrace(Group('NMP', 'CN1CCCC1=O.CN1CCCC1=O', []), '', '[OX1H0]=[C;R][NX3H0;R][CH3]', None),
            NodeTrace(Group('CONH2', 'C(=O)N', [0]), '', '[CX3H0](=[OX1H0])[NX3H2]', None),
            NodeTrace(Group('CONHCH3', 'C(=O)NC', [0]), '', '[OX1H0;!R]=[CX3H0;!R][NH1X3;!R][CH3;!R]', None),
            NodeTrace(Group('CONHCH2', 'C(=O)NC', [0,3]), '', '[CH2X4;!R][NH1X3;!R][CX3H0;!R]=[OX1H0;!R]', None),
        }


class SaftGammaMie(BasisSet):
    def __init__(self):
        super().__init__()
        self.node_traces = {
            NodeTrace(Group("-CH3", "CH3", [0]), "", "[CX4H3]", None),
        }


Libraries = {"saftgm": SaftGammaMie, "joback": Joback, "UNIFAC": Unifac}
