from Grouper import GroupGraph
from Grouper.utils import (
    convert_to_nx, 
    data_to_gg_edge
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
