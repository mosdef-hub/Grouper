"""Module for pre-loading libraries of GroupGraphs."""

# from genGrouper.library import Library, Basis, Joback
from collections import namedtuple
from typing import Dict, Union

from Grouper import Node

# Declaring namedtuple for defining adding nodes and their SMARTS to a library.
NodeTrace = namedtuple("NodeTrace", ["node", "doi", "smarts", "priority"])


class BasisSet(object):
    def __init__(self):
        # trace = (Node(name, smiles, hubs), doi,  smarts, priorty) # could also name moiety
        self.node_traces_ = set()

    def __repr__(self):
        return f"{self.__class__.__name__} has {len(self.node_traces_)} unique nodes."

    def add_node(
        self, node: Node, doi: str = None, smarts: str = None, priority: int = None
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
            NodeTrace(Node(0, "-CH3", "CH3", [0]), "", "CX4H3", None),
            NodeTrace(Node(1, "-CH2-", "CH2", [0, 0]), "", "[!R;CX4H2]", None),
            NodeTrace(Node(2, ">CH-", "CH", [0, 0, 0]), "", "[!R;CX4H]", None),
            NodeTrace(Node(3, ">C<", "C", [0, 0, 0, 0]), "", "[!R;CX4H0]", None),
            NodeTrace(Node(4, "CH2=CH-", "C=C", [0]), "", "[CX3H2][CX3H1]", None),
            NodeTrace(Node(5, "-CH=CH-", "C=C", [0, 1]), "", "[CX3H1][CX3H1]", None),
            NodeTrace(
                Node(6, "=C<", "C", [0, 0, 0]),
                "",
                "[$([!R;#6X3H0]);!$([!R;#6X3H0]=[#8])]",
                None,
            ),
            NodeTrace(Node(7, "=C=", "C", [0, 0]), "", "[$([CX2H0](=*)=*)]", None),
            NodeTrace(Node(8, "CH", "C", [0]), "", "[$([CX2H1]#[!#7])]", None),
            NodeTrace(Node(9, "C", "C", [0, 0]), "", "[$([CX2H0]#[!#7])]", None),
            NodeTrace(Node(10, "ring-CH2-", "C", [0, 0]), "", "[R;CX4H2]", None),
            NodeTrace(Node(11, "ring>CH-", "C", [0, 0, 0]), "", "[R;CX4H]", None),
            NodeTrace(Node(12, "ring>C<", "C", [0, 0, 0, 0]), "", "[R;CX4H0]", None),
            NodeTrace(Node(13, "ring=CH-", "C", [0, 0]), "", "[R;CX3H1,cX3H1]", None),
            NodeTrace(
                Node(14, "ring=C<", "C", [0, 0, 0]),
                "",
                "[$([R;#6X3H0]);!$([R;#6X3H0]=[#8])]",
                None,
            ),
            NodeTrace(Node(15, "-F", "F", [0]), "", "[F;X1]", None),
            NodeTrace(Node(16, "-Cl", "Cl", [0]), "", "[Cl;X1]", None),
            NodeTrace(Node(17, "-Br", "Br", [0]), "", "[Br;X1]", None),
            NodeTrace(Node(18, "-I", "I", [0]), "", "[I;X1]", None),
            NodeTrace(
                Node(19, "-OH (alcohol)", "O", [0]),
                "",
                "[OX2H;!$([OX2H]-[#6]=[O]);!$([OX2H]-a)]",
                None,
            ),
            NodeTrace(Node(20, "-OH (phenol)", "O", [0]), "", "[O;H1;$(O-!@c)]", None),
            NodeTrace(
                Node(21, "-O- (non-ring)", "O", [0, 0]),
                "",
                "[OX2H0;!R;!$([OX2H0]-[#6]=[#8])]",
                None,
            ),
            NodeTrace(
                Node(22, "-O- (ring)", "O", [0, 0]),
                "",
                "[#8X2H0;R;!$([#8X2H0]~[#6]=[#8])]",
                None,
            ),
            NodeTrace(
                Node(23, ">C=O (non-ring)", "C=O", [0, 0]),
                "",
                "[$([CX3H0](=[OX1]));!$([CX3](=[OX1])-[OX2]);!R]=O",
                None,
            ),
            NodeTrace(
                Node(24, ">C=O (ring)", "C=O", [0, 0]),
                "",
                "[$([#6X3H0](=[OX1]));!$([#6X3](=[#8X1])~[#8X2]);R]=O",
                None,
            ),
            NodeTrace(
                Node(25, "O=CH- (aldehyde)", "O=C", [1]), "", "[CH;D2](=O)", None
            ),
            NodeTrace(
                Node(26, "-COOH (acid)", "C(=O)O", [0]), "", "[OX2H]-[C]=O", None
            ),
            NodeTrace(
                Node(27, "-COO- (ester)", "C(=O)O", [0, 2]),
                "",
                "[#6X3H0;!$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]",
                None,
            ),
            NodeTrace(
                Node(28, "=O (other than above)", "O", [0]),
                "",
                "[OX1H0;!$([OX1H0]~[#6X3]);!$([OX1H0]~[#7X3]~[#8])]",
                None,
            ),
            NodeTrace(Node(29, "-NH2", "N", [0]), "", "[NX3H2]", None),
            NodeTrace(Node(30, ">NH (non-ring)", "N", [0, 0]), "", "[NX3H1;!R]", None),
            NodeTrace(Node(31, ">NH (ring)", "N", [0, 0]), "", "[#7X3H1;R]", None),
            NodeTrace(
                Node(32, ">N- (non-ring)", "N", [0, 0, 0]),
                "",
                "[#7X3H0;!$([#7](~O)~O)]",
                None,
            ),
            NodeTrace(Node(33, "-N= (non-ring)", "N", [0, 0]), "", "[#7X2H0;!R]", None),
            NodeTrace(Node(34, "-N= (ring)", "N", [0, 0]), "", "[#7X2H0;R]", None),
            NodeTrace(Node(35, "=NH", "N", [0]), "", "[#7X2H1]", None),
            NodeTrace(Node(36, "-CN", "CN", [0]), "", "[#6X2]#[#7X1H0]", None),
            NodeTrace(
                Node(37, "-NO2", "NO2", [0]),
                "",
                "[$([#7X3,#7X3+][!#8])](=[O])~[O-]",
                None,
            ),
            NodeTrace(Node(38, "-SH", "S", [0]), "", "[SX2H]", None),
            NodeTrace(
                Node(39, "-S- (non-ring)", "S", [0, 0]), "", "[#16X2H0;!R]", None
            ),
            NodeTrace(Node(40, "-S- (ring)", "S", [0, 0]), "", "[#16X2H0;R]", None),
        }


class SaftGammaMie(BasisSet):
    def __init__(self):
        super().__init__()
        self.node_traces = {
            NodeTrace(Node(0, "-CH3", "CH3", [0]), "", "[CX4H3]", None),
        }


class Unifac(BasisSet):
    def __init__(self):
        super().__init__()
        self.node_traces = {
            NodeTrace(Node(0, "-CH3", "CH3", [0]), "", "[CX4H3]", None),
        }


Libraries = {"saftgm": SaftGammaMie, "joback": Joback, "UNIFAC": Unifac}
