import networkx as nx
import rdkit
import rdkit.Chem
import group_selfies
from pysmiles import read_smiles, write_smiles
#   taken from a good answer by Gambit1614 on StackOverflow https://stackoverflow.com/questions/57095809/networkx-connecting-nodes-using-ports

class GroupGraph(nx.Graph):
    """A graph with ports as parts of nodes that can be connected to other ports."""

    def __init__(self, node_types=None):
        super(GroupGraph, self).__init__()
        if node_types is None:
            node_types = {}

        if not isinstance(node_types, dict):
            raise TypeError("node_types must be a dictionary of node types and their ports")
        if len(node_types) != 0:
            holder = next(iter(node_types.keys())) # raises error if empty
            holder_type = type(holder)
            if not all(isinstance(k, holder_type) for k in node_types.keys()):
                raise ValueError("All keys in node_types must be of the same type")
            if not all(isinstance(v, list) for v in node_types.values()):
                raise ValueError("All values in node_types must be lists")
            if not all(isinstance(p, holder_type) for v in node_types.values() for p in v):
                raise ValueError("All values in node_types must be lists of the same type as the keys")
        
        self.node_types = node_types

    def __str__(self):
        return f"Nodes: {['{} ({}) {}'.format(d[0], d[1]['type'], d[1]['ports']) for d in self.nodes.data()]}\n\nEdges: {','.join(str(tuple(*e[-1])) for e in self.edges.data('ports'))}\n\n"

    def __repr__(self):
        return f"GroupGraph({','.join(str(tuple(*e[-1])) for e in self.edges.data('ports'))})"
    
    def __eq__(self, other):
        return self.nodes == other.nodes and self.edges == other.edges
    
    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __bool__(self):
        return len(self.nodes) > 0 and len(self.edges) > 0

    def add_node(self, nodeID, node_type):
        super(GroupGraph, self).add_node(nodeID)
        self.nodes[nodeID]['type'] = node_type
        self.nodes[nodeID]['ports'] = self.node_types[self.nodes[nodeID]['type']]
    
    def add_edge(self, node1, port1, node2, port2):

        edge_ports = []

        for n, p in [(node1, port1), (node2, port2)]:
            # Sanity check to see if the nodes and ports are present in Graph
            if n not in self.nodes:
                raise Exception(f"Node: {p} is not present in Graph")
            if p not in self.nodes(data=True)[n]['ports']:
                raise Exception(f"Port: {p} is incorrect for Node: {n}!")

            edge_ports.append(n + '.' + str(p))

        # Add the port points as edge attributes
        if self.has_edge(node1, node2):
            self.edges[node1, node2]['ports'].append(edge_ports)
        else:
            super(GroupGraph, self).add_edge(node1, node2, ports=[edge_ports])

    def make_undirected(self):
        for edge in self.edges:
            node_port = tuple(self.edges[edge]['ports'][0])
            node_s, node_t = edge
            port_s, port_t = edge['ports'][0], edge['ports'][1]
            self.add_bond(node_t, port_t, node_s, port_s)
            self.edges[edge]['ports'].append(f'{node_port[1]}.{node_port[0]}')
        
    def n_free_ports(self, node):
        # num ports - num edges - num edges with node as target
        n_edges_with_node_target = 0
        if len(self.edges()) == 0:
            return len(self.nodes[node]['ports'])
        for e in self.edges():
            for port_edge in self.edges[e]['ports']:
                if port_edge[-1].split('.')[0] == node:
                    n_edges_with_node_target += 1
        return len(self.nodes[node]['ports']) - len(list(list(self.edges(data=True))[0][-1].values())[0]) - n_edges_with_node_target
    
    def group_graph_from_smiles(self, smiles, groups):
        """
        Generates the group graph from a smiles string by identifying the groups in the molecule and making the nessisary connections between these groups
        """
        vocab_fragment = dict([(f'frag{idx}',group_selfies.Group(f'frag{idx}', g.canonsmiles)) for idx, g in enumerate(groups)])
        grammar_fragment = group_selfies.GroupGrammar(vocab=vocab_fragment)
        mol = rdkit.Chem.MolFromSmiles(smiles)
        extracted_groups = grammar_fragment.extract_groups(mol)
        group_graph = GroupGraph(node_types= { g[0].name : g[0].attachment_points for g in extracted_groups } )
        for i, group in enumerate(extracted_groups):
            group_graph.add_node(f'node{i}', group[0].name)
        for i, group in enumerate(extracted_groups):
            for j, group2 in enumerate(extracted_groups):
                if i != j:
                    for port in group[0].attachment_points:
                        for port2 in group2[0].attachment_points:
                            if port == port2:
                                group_graph.add_edge(f'node{i}', port, f'node{j}', port2)

        return group_graph
    
    def from_mbuild(self, compound, groups):
        """Create a GroupGraph from a mb.Compound and a list of Groups."""
        smiles = compound.get_smiles()
        return self.group_graph_from_smiles(smiles, groups)
    
    def to_vector(self):
        """Convert the GroupGraph to a vector representation."""
        type_to_idx = { t : i for i, t in enumerate(self.node_types.keys()) }
        histogram = [0] * len(self.node_types)
        for n in self.nodes(data=True):
            histogram[type_to_idx[n[1]['type']]] += 1
        return histogram
    
    def to_molecular_graph(self, node_type_to_smiles, node_type_port_to_index):
        """Converts the GroupGraph to a molecular graph using SMILES notation."""
        molecular_graph = nx.Graph()

        atom_id = -1
        node_port_to_atom_index = {}
        for node, data in self.nodes(data=True):
            node_type = data['type']
            ports = data['ports']
            smiles = node_type_to_smiles[node_type]
            mole_graph = read_smiles(smiles)
            node_port_to_atom_index[node] = {}
            group_atom_id = -1
            for atom in mole_graph.nodes(data=True):
                atom_data = atom[1]['element']
                atom_id += 1
                group_atom_id += 1
                if group_atom_id in node_type_port_to_index[node_type].values():
                    for i, port in enumerate(ports):
                        node_port_to_atom_index[node][port] = atom_id
        
        # Need conversion from node and subgraph indices to molecular graph indices
        atom_id = -1
        node_sub_graph_indices_to_molecular_graph_indices = {}
        for node, data in self.nodes(data=True):
            node_sub_graph_indices_to_molecular_graph_indices[node] = {}
            mole_graph = read_smiles(node_type_to_smiles[data['type']])
            for atom in mole_graph.nodes(data=True):
                atom_id += 1
                node_sub_graph_indices_to_molecular_graph_indices[node][atom[0]] = atom_id

                    

        # Iterate over nodes in the GroupGraph
        atom_id = -1
        for node, data in self.nodes(data=True):
            node_type = data['type']
            ports = data['ports']

            # Get SMILES notation for the node type
            smiles = node_type_to_smiles[node_type]

            # Read the SMILES notation to create a molecular graph
            mole_graph = read_smiles(smiles)

            # Add atoms to the molecular graph
            for atom in mole_graph.nodes(data=True):
                atom_id += 1
                molecular_graph.add_node(
                    atom_id, 
                    element=atom[1]['element'], 
                    charge=atom[1]['charge'], 
                    aromatic=atom[1]['aromatic'], 
                    # hcount=atom[1]['hcount']
                    hcount=0
                )
            
            # Add bonds from the subgraph to the molecular graph
            for bond in mole_graph.edges(data=True):
                sub_atom1, sub_atom2, data = bond
                atom1 = node_sub_graph_indices_to_molecular_graph_indices[node][sub_atom1]
                atom2 = node_sub_graph_indices_to_molecular_graph_indices[node][sub_atom2]
                molecular_graph.add_edge(atom1, atom2, order=data['order'])
                

        # Iterate over edges in the GroupGraph
        for edge in self.edges(data=True):
            node1, node2, data = edge
            edge_ports_list = data['ports']

            # Add bonds to the molecular graph
            for edge_ports in edge_ports_list:
                node1, node2 = edge_ports
                node1, port1 = node1.split('.')[0], node1.split('.')[1]
                node2, port2 = node2.split('.')[0], node2.split('.')[1]
                # need to convert node_port_index to index in the molecular graph
                atom1 = node_port_to_atom_index[node1][port1]
                atom2 = node_port_to_atom_index[node2][port2]
                molecular_graph.add_edge(atom1, atom2, order=1)

        return molecular_graph
    
    def to_smiles(self, node_type_to_smiles, node_type_port_to_index):
        """Converts the GroupGraph to a SMILES notation."""
        molecular_graph = self.to_molecular_graph(node_type_to_smiles, node_type_port_to_index)
        return write_smiles(molecular_graph)
    
    def to_PyG_Data(self, node_descriptor_generater, max_n_attachments):
        """Convert the GroupGraph to a data representation."""
        from torch_geometric.data import Data
        import torch
        one_hot_vector = lambda index, num_classes : torch.eye(num_classes)[index]

        # Create the node features
        dummy_feature = node_descriptor_generater(list(self.node_types.keys())[0])
        node_features = torch.zeros(len(self.nodes), len(dummy_feature))
        for i, n in enumerate(self.nodes(data=True)):
            node_features[i] = node_descriptor_generater(n[1]['type'])
        
        # Create the edge index
        edge_index = torch.zeros(2, len(self.edges))
        for i, e in enumerate(self.edges):
            edge_index[0, i] = list(self.nodes).index(e[0])
            edge_index[1, i] = list(self.nodes).index(e[1])
        
        # Create the edge features
        edge_features = torch.zeros(len(self.edges), max_n_attachments*2)

        for i, e in enumerate(self.edges(data=True)):
            
            # get nodes and ports
            node_ports = [node_port.split('.') for node_port in e[2]['ports'][0]]
            # convert the ports to one-hot vectors
            port_index_s = self.nodes(data=True)[node_ports[0][0]]['ports'].index(node_ports[0][1])
            port_index_t = self.nodes(data=True)[node_ports[1][0]]['ports'].index(node_ports[1][1])
            edge_features[i] = torch.cat([
                one_hot_vector(port_index_s, max_n_attachments), 
                one_hot_vector(port_index_t, max_n_attachments)
            ])

        return Data(x=node_features, edge_index=edge_index, edge_attr=edge_features)
            
