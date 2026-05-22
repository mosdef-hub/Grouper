from Grouper import GroupGraph
from Grouper.utils import (
    convert_to_nx,
    data_to_gg_edge,
    parse_vcolg_line,
)
from Grouper.tests.base_test import BaseTest
import networkx as nx
import pytest
import torch
import rdkit.Chem
from rdkit.Chem import Descriptors
from math import isnan
import torch_geometric
import torch_geometric.data


def node_descriptor_generator(node_smiles):
    mol = rdkit.Chem.MolFromSmiles(node_smiles)

    desc = Descriptors.CalcMolDescriptors(mol)
    desc = {k: desc[k] for k in desc.keys() if (not isnan(desc[k]) or desc[k] is not None)}
    choosen_keys = [
        "MaxPartialCharge", "MinPartialCharge"," MaxAbsPartialCharge", "MinAbsPartialCharge", # Partial charge descriptors
        'BCUT2D_CHGHI, BCUT2D_CHGLO', # Burden eigenvalues from 2D matrix weighted by atomic charges
        'BCUT2D_LOGPHI', # Burden eigenvalues from 2D matrix weighted by logP
        "BCUT2D_MRHI", "BCUT2D_MRLO", # Burden eigenvalues from 2D matrix weighted by molar refractivity
        "BCUT2D_MWHI", "BCUT2D_MWLO", # Burden eigenvalues from 2D matrix weighted by molecular weight
        "Chi0", "Chi0n", "Chi0v", "Chi1", "Chi1n", "Chi1v", "Chi2n", "Chi2v", "Chi3n", "Chi3v", "Chi4n", "Chi4v", # Kier and Hall molecular connectivity indices
        "MolLogP", # octanol-water partition coefficient
        "LabuteASA", # Labute's Approximate Surface Area
        "TPSA", # topological polar surface area
        "PEOE_VSA1", "PEOE_VSA10", "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13", "PEOE_VSA14", "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5", "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9", # PEOE_VSA properties
        "MaxEStateIndex", "MinEStateIndex", # E-state indices
        "NumHAcceptors", "NumHDonors", # Number of hydrogen bond acceptors and donors
    ]
    # flatten descriptors into single vector
    desc = [v for k,v in desc.items() if k in choosen_keys]
    desc = torch.tensor(desc, dtype=torch.float64)
    return desc


class TestUtils(BaseTest):
    @pytest.mark.parametrize(
        "graph_fixture", ["basic_graph", "single_node_graph", "single_edge_graph", "five_member_ring_graph"]
    )
    def test_convert_to_nx(self, request, graph_fixture):
        graph = request.getfixturevalue(graph_fixture)
        nxgraph = convert_to_nx(graph)
        assert isinstance(nxgraph, nx.Graph)

    def test_to_data(self):
        graph = GroupGraph()
        graph.add_node("oxygen", "O", hubs = [0,0])
        graph.add_node("oxygen", "O", hubs = [0,0])
        graph.add_node("carbon", "C", hubs = [0,0,0,0])
        graph.add_node("carbon", "C", hubs = [0,0,0,0])
        graph.add_node("nitrogen", "N", hubs = [0,0,0])
        graph.add_edge((2, 0), (4, 0), 1)
        graph.add_edge((1, 0), (4, 2), 1)
        graph.add_edge((0, 0), (4, 1), 1)
        graph.add_edge((0, 1), (3, 0), 1)
        max_ports = max([len(group.hubs) for group in graph.nodes.values()])
        nxG = convert_to_nx(graph)
        data = nxG.to_PyG_Data(
            node_descriptor_generator = node_descriptor_generator,
            max_ports = max_ports
        )
        assert isinstance(data, torch_geometric.data.Data)
        edges = data_to_gg_edge(data, max_ports)
        # Add bond order
        edges = [(e[0], e[1], e[2], e[3], 1) for e in edges]

        # Check that the edges are the same
        assert set(edges) == set(graph.edges)

        graph = GroupGraph()
        graph.add_node("carbon", "C", hubs = [0,0,0,0])
        graph.add_node("oxygen", "O", hubs = [0,0])
        graph.add_node("carbon", "C", hubs = [0,0,0,0])
        graph.add_node("nitrogen", "N", hubs = [0,0,0])
        graph.add_node("nitrogen", "N", hubs = [0,0,0])
        graph.add_edge((2,1),(4,2),1)
        graph.add_edge((1,0),(4,0),1)
        graph.add_edge((0,1),(4,1),1)
        graph.add_edge((0,0),(3,0),1)

        max_ports = max([len(group.hubs) for group in graph.nodes.values()])
        nxG = convert_to_nx(graph)
        data = nxG.to_PyG_Data(
            node_descriptor_generator = node_descriptor_generator,
            max_ports = max_ports + 1
        )
        assert isinstance(data, torch_geometric.data.Data)
        edges = data_to_gg_edge(data, max_ports + 1)
        # Add bond order
        edges = [(e[0], e[1], e[2], e[3], 1) for e in edges]
        assert set(edges) == set(graph.edges)


class TestParseVcolgLine:
    """Public parser for vcolg's terse-format (``-T``) lines.

    The pin-down values come from running ``geng -c <n> | vcolg -T -m<c>``
    on small inputs and copying the resulting lines verbatim.
    """

    def test_simple_line(self):
        """Three-node path, single color — n_vertices, colors, edges
        all parse cleanly."""
        # `geng -c 3 | vcolg -T -m1` on the path graph
        n, colors, edges = parse_vcolg_line("3 2 0 0 0   0 1 1 2")
        assert n == 3
        assert colors == [0, 0, 0]
        assert edges == [(0, 1), (1, 2)]

    def test_multi_color_line(self):
        """Multiple colors land in the right slots, edges still parse
        correctly. From a typical n=4, c=2 vcolg sample."""
        n, colors, edges = parse_vcolg_line("4 3 0 1 0 1   0 1 1 2 2 3")
        assert n == 4
        assert colors == [0, 1, 0, 1]
        assert edges == [(0, 1), (1, 2), (2, 3)]

    def test_no_edges_line(self):
        """An isolated single node has no edges — make sure the empty
        edge list parses cleanly (the line still has the double-space
        separator, just nothing after it)."""
        n, colors, edges = parse_vcolg_line("1 0 0  ")
        assert n == 1
        assert colors == [0]
        assert edges == []

    def test_missing_separator_raises(self):
        """A line without the double-space header/edge separator is
        malformed and the parser must say so clearly rather than
        silently returning garbage."""
        with pytest.raises(ValueError, match="double-space"):
            parse_vcolg_line("3 2 0 0 0 0 1 1 2")  # only single spaces

    def test_odd_edge_count_raises(self):
        """vcolg always emits edges as (src, dst) pairs — an odd number
        of edge tokens means the input is corrupted, not that the last
        edge is half-formed."""
        with pytest.raises(ValueError, match="odd token count"):
            parse_vcolg_line("3 2 0 0 0   0 1 1")  # 3 edge tokens

    def test_round_trip_with_exhaustive_generate(self):
        """A regression check that hardwires us against the actual
        vcolg format: shell out to vcolg on a tiny input and confirm
        the parser handles every line.

        Belt-and-suspenders since the simple-line tests above already
        pin the format, but if a future vcolg release changes the
        spacing or column order this catches it.
        """
        import os
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            geng_path = os.path.join(td, "g.txt")
            vcolg_path = os.path.join(td, "v.txt")
            os.system(f"geng -c 3 > {geng_path} 2>/dev/null")
            os.system(f"vcolg {geng_path} -T -m2 > {vcolg_path} 2>/dev/null")
            with open(vcolg_path) as f:
                lines = [ln.rstrip("\n") for ln in f if ln.strip()]
        assert lines, "vcolg produced no output — is nauty on PATH?"
        for line in lines:
            n, colors, edges = parse_vcolg_line(line)
            assert n == 3
            assert len(colors) == n
            # Every edge endpoint must reference a valid node index.
            for src, dst in edges:
                assert 0 <= src < n and 0 <= dst < n
