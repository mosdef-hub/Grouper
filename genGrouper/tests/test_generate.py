from mbuild.tests.base_test import BaseTest
from genGrouper import GroupGraph, exhaustive_generate 
# from pysmiles import write_smiles
# import pytest

class TestGroupGraph(BaseTest):
    def test_chiral_group_graph_generation(self):
        node_defs = set()
        node_defs.add(('OH', 'O', ['C11', 'C12'], ['O1']))
        node_defs.add(('CC', 'C', ['C11', 'C12'], ['C21', 'C22']))
        out = exhaustive_generate(3, node_defs)

        def check_if_graph_has_ports(graph, ports):
            all_ports = set()
            for e in graph.edges(data=True):
                for p in e[2]['ports']:
                    all_ports.add(p[0].split('.')[1])
                    all_ports.add(p[1].split('.')[1])
            return all_ports == ports

        cc_2_oh = []
        for g in out:
            n_OH, n_CC = 0, 0
            for name, info in g.nodes(data=True):
                if info['type'] == 'OH':
                    n_OH += 1
                if info['type'] == 'CC':
                    n_CC += 1
            if n_OH == 2 and n_CC == 1:
                cc_2_oh.append(g)

        # cis TODO verifiy actually cis 
        assert any([True if check_if_graph_has_ports(g, ports={'C11', 'C22', 'O1'}) else False for g in cc_2_oh])

        # trans
        assert any([True if check_if_graph_has_ports(g, ports={'C11', 'C21', 'O1'}) else False for g in cc_2_oh])

        # terminal alkene
        assert any([True if check_if_graph_has_ports(g, ports={'C11', 'C12', 'O1'}) else False for g in cc_2_oh])

