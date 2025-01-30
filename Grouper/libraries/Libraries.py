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
            NodeTrace(Group("[CX4;H3;!R]", "C", [0]), "", "CH3", None),
            NodeTrace(Group("[CX4;H2;!R]", "C", [0, 0]), "", "CH2", None),
            NodeTrace(Group("[CX4;H1;!R]", "C", [0, 0, 0]), "", "CH", None),
            NodeTrace(Group("[CX4;H0;!R]", "C", [0, 0, 0, 0]), "", "C", None),
            NodeTrace(Group("[CX3;H2]=[CX3;H1]", "C=C", [0]), "", "CH2=CH", None),
            NodeTrace(Group("[CX3;H1]=[CX3;H1]", "C=C", [0, 1]), "", "CH=CH", None),
            NodeTrace(Group("[CX3;H2]=[CX3;H0]", "C=C", [0, 0]), "", "CH2=C", None),
            NodeTrace(Group("[CX3;H1]=[CX3;H0]", "C=C", [0, 0, 1]), "", "CH=C", None),
            NodeTrace(Group("[OH1;$([OH1][CX4H2])]", "O", [0]), "", "OH (P)", None),
            NodeTrace(Group("[CX4;H3][OX2;H1]", "CO", [0]), "", "CH3OH", None),
            NodeTrace(Group("[OH2]", "O", []), "", "H2O", None),
            NodeTrace(Group("[cX3;H0;R][OX2;H1]", "C=O", [0, 0]), "", "ACOH", None),
            NodeTrace(Group("[CX4;H3][CX3;!H1](=O)", "CC=O", []), "", "CH3CO", None),
            NodeTrace(Group("[CX4;H2][CX3;!H1](=O)", "CC=O", [0]), "", "CH2CO", None),
            NodeTrace(Group("[CX3H1](=O)", "C=O", [0, 0]), "", "CHO", None),
            NodeTrace(
                Group("[CH3][CX3;H0](=[O])[OH0]", "CC(=O)O", [3]), "", "CH3COO", None
            ),
            NodeTrace(
                Group("[CX4;H2][CX3](=[OX1])[OX2]", "CC(=O)O", [0, 3]),
                "",
                "CH2COO",
                None,
            ),
            NodeTrace(
                Group("[CX3;H1](=[OX1])[OX2]", "C(=O)O", [0, 3]), "", "HCOO", None
            ),
            NodeTrace(Group("[CH3;!R][OH0;!R]", "CO", [0, 1]), "", "CH3O", None),
            NodeTrace(Group("[CH2;!R][OH0;!R]", "CO", [0, 1]), "", "CH2O", None),
            NodeTrace(Group("[C;H1;!R][OH0;!R]", "CO", [0, 1]), "", "CHO", None),
            NodeTrace(Group("[cX3;H1]", "C=O", [0]), "", "ACH", None),
            NodeTrace(Group("[cX3;H0]", "C=O", [0]), "", "AC", None),
            NodeTrace(Group("[cX3;H0][CX4;H3]", "C(=O)C", [0]), "", "ACCH3", None),
            NodeTrace(Group("[cX3;H0][CX4;H2]", "C(=O)C", [0, 2]), "", "ACCH2", None),
            NodeTrace(Group("[cX3;H0][CX4;H1]", "C(=O)C", [0, 2, 2]), "", "ACCH", None),
            NodeTrace(
                Group("[CX4;H2;R;$(C(C)OCC)][OX2;R][CX4;H2;R]", "C1CCOC1", []),
                "",
                "THF",
                None,
            ),
            NodeTrace(Group("[CX4;H3][NX3;H2]", "CN", []), "", "CH3NH2", None),
            NodeTrace(Group("[CX4;H2][NX3;H2]", "CN", [0]), "", "CH2NH2", None),
            NodeTrace(Group("[CX4;H1][NX3;H2]", "CN", [0, 0]), "", "CHNH2", None),
            NodeTrace(Group("[CX4;H3][NX3;H1]", "CN", [1]), "", "CH3NH", None),
            NodeTrace(Group("[CX4;H2][NX3;H1]", "CN", [0, 1]), "", "CH2NH", None),
            NodeTrace(Group("[CX4;H1][NX3;H1]", "CN", [0, 0, 1]), "", "CHNH", None),
            NodeTrace(Group("[CX4;H3][NX3;H0]", "CN", [1, 1]), "", "CH3N", None),
            NodeTrace(Group("[CX4;H2][NX3;H0]", "CN", [0, 1, 1]), "", "CH2N", None),
            NodeTrace(Group("[c][NX3;H2]", "C(=O)N", [0]), "", "ACNH2", None),
            NodeTrace(Group("[cX3H1][n][cX3H1]", "C(=O)NC=O", [2]), "", "AC2H2N", None),
            NodeTrace(
                Group("[cX3H0][n][cX3H1]", "C(=O)NC=O", [0, 2]), "", "AC2HN", None
            ),
            NodeTrace(
                Group("[cX3H0][n][cX3H0]", "C(=O)NC=O", [0, 2, 3]), "", "AC2N", None
            ),
            NodeTrace(Group("[CX4;H3][CX2]#[NX1]", "CC#N", []), "", "CH3CN", None),
            NodeTrace(Group("[CX4;H2][CX2]#[NX1]", "CC#N", [0]), "", "CH2CN", None),
            NodeTrace(
                Group("[CX3,cX3](=[OX1])[OX2H0,oX2H0]", "C(=O)O", [0, 1]),
                "",
                "COO",
                None,
            ),
            NodeTrace(Group("[CX3](=[OX1])[O;H1]", "C(=O)O", [0]), "", "COOH", None),
            NodeTrace(
                Group("[CX3;H1](=[OX1])[OX2;H1]", "C(=O)O", []), "", "HCOOH", None
            ),
            NodeTrace(
                Group("[CX4;H2;!$(C(Cl)(Cl))](Cl)", "C(Cl)", [0]), "", "CH2CL", None
            ),
            NodeTrace(
                Group("[CX4;H1;!$(C(Cl)(Cl))](Cl)", "C(Cl)", [0, 0]),
                "",
                "CHCL",
                None,
            ),
            NodeTrace(Group("[CX4;H0](Cl)", "C(Cl)", [0, 0, 0]), "", "CCL", None),
            NodeTrace(
                Group("[CX4;H2;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)", "C(Cl)(Cl)", []),
                "",
                "CH2CL2",
                None,
            ),
            NodeTrace(
                Group("[CX4;H1;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)", "C(Cl)(Cl)", [0]),
                "",
                "CHCL2",
                None,
            ),
            NodeTrace(
                Group("[CX4;H0;!$(C(Cl)(Cl)(Cl))](Cl)(Cl)", "C(Cl)(Cl)", [0, 0]),
                "",
                "CCL2",
                None,
            ),
            NodeTrace(
                Group(
                    "[CX4;H1;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)",
                    "C(Cl)(Cl)(Cl)",
                    [],
                ),
                "",
                "CHCL3",
                None,
            ),
            NodeTrace(
                Group(
                    "[CX4;H0;!$([CX4;H0](Cl)(Cl)(Cl)(Cl))](Cl)(Cl)(Cl)",
                    "C(Cl)(Cl)(Cl)",
                    [0],
                ),
                "",
                "CCL3",
                None,
            ),
            NodeTrace(

                Group("[CX4;H0]([Cl])([Cl])([Cl])([Cl])", "C(Cl)(Cl)(Cl)(Cl", []),
                "",
                "CCL4",
                None,
            ),
            NodeTrace(Group("[c][Cl]", "C(=O)(Cl)", [0]), "", "ACCL", None),
            NodeTrace(

                Group("[CX4;H3][NX3](=[OX1])([OX1])", "CN(=O)O", []),
                "",
                "CH3NO2",
                None,
            ),
            NodeTrace(
                Group("[CX4;H2][NX3](=[OX1])([OX1])", "CN(=O)O", [0]),
                "",
                "CH2NO2",
                None,
            ),
            NodeTrace(
                Group("[CX4;H1][NX3](=[OX1])([OX1])", "CN(=O)O", [0, 0]),
                "",
                "CHNO2",
                None,
            ),
            NodeTrace(
                Group("[cX3][NX3](=[OX1])([OX1])", "C(=O)N(=O)O", [0]),
                "",
                "ACNO2",
                None,
            ),
            NodeTrace(Group("C(=S)=S", "S=C=S", []), "", "CS2", None),
            NodeTrace(Group("[SX2H][CX4;H3]", "CS", []), "", "CH3SH", None),
            NodeTrace(Group("[SX2H][CX4;H2]", "CS", [0]), "", "CH2SH", None),
            NodeTrace(Group("c1cc(oc1)C=O", "c1cc(oc1)C=O", []), "", "FURFURAL", None),
            NodeTrace(
                Group("[OX2;H1][CX4;H2][CX4;H2][OX2;H1]", "OCCO", []),
                "",
                "DOH",
                None,
            ),
            NodeTrace(Group("[I]", "I", [0]), "", "I", None),
            NodeTrace(Group("[Br]", "(Br)", [0]), "", "BR", None),
            NodeTrace(Group("[CX2;H1]#[CX2;H0]", "C#C", [0]), "", "CH=-C", None),
            NodeTrace(Group("[CX2;H0]#[CX2;H0]", "C#C", [0, 1]), "", "C=-C", None),
            NodeTrace(
                Group("[SX3H0](=[OX1])([CX4;H3])[CX4;H3]", "CS(=O)C", []),
                "",
                "DMSO",
                None,
            ),
            NodeTrace(
                Group("[CX3;H2]=[CX3;H1][CX2;H0]#[NX1;H0]", "C=CC#N", []),
                "",
                "ACRY",
                None,
            ),
            NodeTrace(
                Group("[$([Cl;H0]([C]=[C]))]", "(Cl)C=C", [1, 2, 2]),
                "",
                "CL-(C=C)",
                None,
            ),
            NodeTrace(Group("[CX3;H0]=[CX3;H0]", "C=C", [0, 0, 1, 1]), "", "C=C", None),
            NodeTrace(Group("[cX3][F]", "C(=O)F", [0]), "", "ACF", None),
            NodeTrace(
                Group("[CX4;H3][N]([CX4;H3])[CX3;H1]=[O]", "CN(C)C=O", []),
                "",
                "DMF",
                None,
            ),
            NodeTrace(
                Group("[NX3]([CX4;H2])([CX4;H2])[CX3;H1](=[OX1])", "C(=O)N(C)C", []),
                "",
                "HCON(CH2)2",
                None,
            ),
            NodeTrace(Group("C(F)(F)F", "C(F)(F)F", [0]), "", "CF3", None),
            NodeTrace(Group("C(F)F", "C(F)F", [0, 0]), "", "CF2", None),
            NodeTrace(Group("C(F)", "CF", [0, 0, 0]), "", "CF", None),
            NodeTrace(Group("[CH2;R]", "C", [0, 0]), "", "CY-CH2", None),
            NodeTrace(Group("[CH1;R]", "C", [0, 0, 0]), "", "CY-CH", None),
            NodeTrace(Group("[CH0;R]", "C", [0, 0, 0, 0]), "", "CY-C", None),
            NodeTrace(Group("[OH1;$([OH1][CX4H1])]", "CO", [0, 0]), "", "OH (S)", None),
            NodeTrace(
                Group("[OH1;$([OH1][CX4H0])]", "CO", [0, 0, 0]), "", "OH (T)", None
            ),
            NodeTrace(
                Group(
                    "[CX4H2;R][OX2;R;$(O(CC)C)][CX4H2;R][OX2;R][CX4H2;R]",
                    "CO",
                    [0, 1],
                ),
                "",
                "CY-CH2O",
                None,
            ),
            NodeTrace(
                Group("[CX4H2;R][OX2;R;$(O(C)C)]", "C1OCOCO1", []),
                "",
                "TRIOXAN",
                None,
            ),
            NodeTrace(Group("[CX4H0][NH2]", "CN", [0, 0, 0]), "", "CNH2", None),
            NodeTrace(
                Group("[OX1H0]=[C;R][NX3H0;R][CH3]", "CN1CCCC1=O.CN1CCCC1=O", []),
                "",
                "NMP",
                None,
            ),
            NodeTrace(
                Group("[CX3H0](=[OX1H0])[NX3H2]", "C(=O)N", [0]), "", "CONH2", None
            ),
            NodeTrace(
                Group("[OX1H0;!R]=[CX3H0;!R][NH1X3;!R][CH3;!R]", "C(=O)NC", [0]),
                "",
                "CONHCH3",
                None,
            ),
            NodeTrace(
                Group("[CH2X4;!R][NH1X3;!R][CX3H0;!R]=[OX1H0;!R]", "C(=O)NC", [0, 3]),
                "",
                "CONHCH2",
                None,
            ),
        }


class SaftGammaMie(BasisSet):
    def __init__(self):
        super().__init__()
        self.node_traces = {
            NodeTrace(Group("-CH3", "CH3", [0]), "", "[CX4H3]", None),
        }


Libraries = {"saftgm": SaftGammaMie, "joback": Joback, "UNIFAC": Unifac}
