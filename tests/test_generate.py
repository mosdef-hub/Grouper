from mbuild.tests.base_test import BaseTest
from molGrouper import GroupGraph
from molGrouper import generate_group_graph_space
import pytest

class TestGroupGraph(BaseTest):

    @pytest.fixture(autouse=True)
    def test_creation(self):
        # Define node types with ports
        self.node_types = {
            'NH2': ['N1'], # amine
            'CO': ['C1', 'C2'], # carbonyl
            'CC': ['C11', 'C12', 'C21', 'C22'], # alkene
        }
        self.node_type_to_smiles = {
            'NH2': 'N',
            'CO': 'C=O',
            'CC': 'C=C',
        }
        self.node_port_to_atom_index = {
            'NH2': {'N1': 0},
            'CO': {'C1': 0, 'C2': 0},
            'CC': {'C11': 0, 'C12': 0, 'C21': 1, 'C22': 1},
        }

        # Create an instance of GroupGraph for testing
        self.graph = GroupGraph(self.node_types)

    def test_generate_group_graph_space(self):
        out = generate_group_graph_space(3, self.node_types)
        # assert len(out) == 539
        pass

    def test_single_node_type(self):
        # Define node types with ports
        node_types = {
            'CC': ['C11', 'C12', 'C13', 'C21', 'C22', 'C23'] #alkane
        }
        node_types_to_smiles = {
            'CC': 'CC'
        }
        node_port_to_atom_index = {
            'CC': {'C11': 0, 'C12': 0, 'C13': 0, 'C21': 1, 'C22': 1, 'C23': 1}
        }

        graph = GroupGraph(node_types)

        out = generate_group_graph_space(3, node_types)


