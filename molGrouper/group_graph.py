import networkx as nx
import rdkit
import rdkit.Chem
import group_selfies
from pysmiles import read_smiles, write_smiles
import mbuild
import networkx 
import torch_geometric
from typing import List, Dict, Tuple, Union, Any, Callable, Sequence
#   taken from a good answer by Gambit1614 on StackOverflow https://stackoverflow.com/questions/57095809/networkx-connecting-nodes-using-ports

class GroupGraph(nx.Graph):
    """A graph with ports as parts of nodes that can be connected to other ports."""

    def __init__(self, node_types: Dict[str, List] = None):
        """
        Initialize a GroupGraph.

        Parameters:
        - node_types (dict): Dictionary of node types and their corresponding ports.

        Raises:
        - TypeError: If node_types is not a dictionary.
        - ValueError: If keys in node_types are not of the same type or values are not lists.
        """
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

        for k, v in node_types.items(): # check if ports are unique
            if len(v) != len(set(v)):
                raise ValueError("All ports in node_types must be unique")
            
        self.node_types = node_types

    def __str__(self) -> str:
        """
        Return a string representation of the GroupGraph.

        Returns:
        - str: String representation of nodes and edges.
        """
        return f"Nodes: {['{} ({}) {}'.format(d[0], d[1]['type'], d[1]['ports']) for d in self.nodes.data()]}\nEdges: {','.join(str(tuple(*e[-1])) for e in self.edges.data('ports'))}\n"

    def __repr__(self) -> str:
        """
        Return a representation of the GroupGraph.

        Returns:
        - str: Representation of the GroupGraph.
        """
        return f"GroupGraph({','.join(str(tuple(*e[-1])) for e in self.edges.data('ports'))})"
    
    def __eq__(self, other) -> bool:
        """
        Check if two GroupGraph objects are equal.

        Parameters:
        - other: Another GroupGraph object.

        Returns:
        - bool: True if equal, False otherwise.
        """
        # return self.nodes(data=True) == other.nodes(data=True) and self.edges(data=True) == other.edges(data=True)
        return nx.utils.misc.graphs_equal(self, other)
    
    def __ne__(self, other) -> bool:
        """
        Check if two GroupGraph objects are not equal.

        Parameters:
        - other: Another GroupGraph object.

        Returns:
        - bool: True if not equal, False otherwise.
        """
        return not self.__eq__(other)
    
    def __bool__(self) -> bool:
        """
        Check if the GroupGraph is non-empty.

        Returns:
        - bool: True if non-empty, False otherwise.
        """
        return len(self.nodes) > 0 and len(self.edges) > 0

    def add_node(self, nodeID: Any, node_type: str) -> None:
        """
        Add a node to the GroupGraph.

        Parameters:
        - nodeID: Identifier for the node.
        - node_type: Type of the node.

        Raises:
        - Exception: If the node is already present in the Graph.
        """
        # Sanity check to see if the node is already present in Graph
        if nodeID in self.nodes:
            raise AttributeError(f"Node: {nodeID} is already present in Graph")
        super(GroupGraph, self).add_node(nodeID)
        self.nodes[nodeID]['type'] = node_type
        self.nodes[nodeID]['ports'] = self.node_types[self.nodes[nodeID]['type']]
    
    def add_edge(self, node1: Any, port1: Any, node2: Any, port2: Any) -> None:
        """
        Add an edge between two nodes in the GroupGraph.

        Parameters:
        - node1: First node.
        - port1: Port on the first node.
        - node2: Second node.
        - port2: Port on the second node.

        Raises:
        - Exception: If a node has no free ports or if nodes or ports are not present in the Graph.
        """
        if self.n_free_ports(node1) <= 0:
            raise AttributeError(f"Node: {node1} has no free ports!")
        if self.n_free_ports(node2) <= 0:
            raise AttributeError(f"Node: {node2} has no free ports!")
        # check if port is already occupied
        for edge in self.edges(data=True):
            for node_port in edge[-1]['ports']:
                src_node, src_port = node_port[0].split('.')
                dst_node, dst_port = node_port[1].split('.')
                if src_node == node1:
                    if src_port == port1:
                        raise AttributeError(f"Node: {node1}.{port1} is already occupied!")
                if src_node == node2:
                    if src_port == port2:
                        raise AttributeError(f"Node: {node2}.{port2} is already occupied!")
                if dst_node == node1:
                    if dst_port == port1:
                        raise AttributeError(f"Node: {node1}.{port1} is already occupied!")
                if dst_node == node2:
                    if dst_port == port2:
                        raise AttributeError(f"Node: {node2}.{port2} is already occupied!")
                

        edge_ports = []

        for n, p in [(node1, port1), (node2, port2)]:
            # Sanity check to see if the nodes and ports are present in Graph
            if n not in self.nodes:
                raise AttributeError(f"Node: {p} is not present in Graph")
            if p not in self.nodes(data=True)[n]['ports']:
                raise AttributeError(f"Port: {p} is incorrect for Node: {n}!")

            edge_ports.append(str(n) + '.' + str(p))

        # Add the port points as edge attributes
        if self.has_edge(node1, node2):
            self.edges[node1, node2]['ports'].append(edge_ports)
        else:
            super(GroupGraph, self).add_edge(node1, node2, ports=[edge_ports])

    def make_undirected(self):
        """
        Convert the GroupGraph to an undirected graph.

        Adds bonds and updates edge ports accordingly.
        """
        for edge in self.edges:
            node_port = tuple(self.edges[edge]['ports'][0])
            node_s, node_t = edge
            port_s, port_t = edge['ports'][0], edge['ports'][1]
            self.add_bond(node_t, port_t, node_s, port_s)
            self.edges[edge]['ports'].append(f'{node_port[1]}.{node_port[0]}')
        
    def n_free_ports(self, nodeID: Any) -> int:
        """
        Get the number of free ports on a node.

        Parameters:
        - nodeID: Identifier for the node.

        Returns:
        - int: Number of free ports.
        """
        # num ports - num edges - num edges with node as target
        occupied_ports = 0
        if len(self.edges()) == 0:
            return len(self.nodes[nodeID]['ports'])
        for e in self.edges(data=True):
            for node_port in e[-1]['ports']:
                if node_port[0].split('.')[0] == nodeID or node_port[1].split('.')[0] == nodeID:
                    occupied_ports += 1
        return len(self.nodes[nodeID]['ports']) - occupied_ports # total number of ports - number of ports with edges
    
    def group_graph_from_smiles(self, smiles: str, groups: List[group_selfies.Group]) -> 'GroupGraph':
        """
        Generate a GroupGraph from a SMILES string and a list of group_selfies.Groups.

        Parameters:
        - smiles (str): SMILES string.
        - groups (List[group_selfies.Group]): List of group_selfies.Groups.

        Returns:
        - GroupGraph: Generated GroupGraph.
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
    
    def from_mbuild(self, compound: mbuild.Compound, groups: List[group_selfies.Group]) -> 'GroupGraph':
        """
        Create a GroupGraph from an mbuild.Compound and a list of group_selfies.Groups.

        Parameters:
        - compound (mbuild.Compound): Input mbuild.Compound.
        - groups (List[group_selfies.Group]): List of group_selfies.Groups.

        Returns:
        - GroupGraph: Generated GroupGraph.
        """
        smiles = compound.get_smiles()
        return self.group_graph_from_smiles(smiles, groups)
    
    def to_vector(self) -> List[int]:
        """
        Convert the GroupGraph to a vector representation.

        Returns:
        - List[int]: Vector representation of the GroupGraph.
        """
        type_to_idx = { t : i for i, t in enumerate(self.node_types.keys()) }
        histogram = [0] * len(self.node_types)
        for n in self.nodes(data=True):
            histogram[type_to_idx[n[1]['type']]] += 1
        return histogram
    
    def to_molecular_graph(
            self, 
            node_type_to_smiles: Dict[str, str], 
            node_type_port_to_index: Dict[str , Dict[str, int]]
    ) -> networkx.Graph:
        """
        Convert the GroupGraph to a molecular graph using SMILES notation.

        Parameters:
        - node_type_to_smiles (Dict[str, str]): Mapping from node type to SMILES notation.
        - node_type_port_to_index (Dict[str, Dict[any, int]]): Mapping from node type and port to index.

        Returns:
        - networkx.Graph: Molecular graph.
        """
        molecular_graph = nx.Graph()

        # Need conversion from node and port in the group graph to atom index in the molecular graph
        node_port_to_atom_index = {}
        atom_count = 0
        for node, data in self.nodes(data=True):
            node_type = data['type']
            ports = data['ports']

            smiles = node_type_to_smiles[node_type]
            mole_graph = read_smiles(smiles)

            node_port_to_atom_index[str(node)] = {}
            for i, port in enumerate(ports):
                node_port_to_atom_index[str(node)][str(port)] = atom_count + node_type_port_to_index[node_type][port]
            
            atom_count += len(mole_graph.nodes)
        
        # Need conversion from node and subgraph indices to molecular graph indices
        atom_id = -1
        node_sub_graph_indices_to_molecular_graph_indices = {}
        for node, data in self.nodes(data=True):
            node_sub_graph_indices_to_molecular_graph_indices[node] = {}
            mole_graph = read_smiles(node_type_to_smiles[data['type']])
            for atom in mole_graph.nodes(data=True):
                atom_id += 1
                node_sub_graph_indices_to_molecular_graph_indices[node][atom[0]] = atom_id

        # Add atoms and bonds from the subgraphs to the molecular graph
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
                

        # Add bonds from the group graph to the molecular graph
        for edge in self.edges(data=True):
            node1, node2, data = edge
            edge_ports_list = data['ports']

            # Add bonds to the molecular graph
            #TODO check if port1 and port2 are of same marriage group
            for edge_ports in edge_ports_list:
                node1, node2 = edge_ports
                node1, port1 = node1.split('.')[0], node1.split('.')[1]
                node2, port2 = node2.split('.')[0], node2.split('.')[1]
                # need to convert node_port_index to index in the molecular graph
                atom1 = node_port_to_atom_index[node1][port1]
                atom2 = node_port_to_atom_index[node2][port2]
                molecular_graph.add_edge(atom1, atom2, order=1)

        return molecular_graph
    
    def to_smiles(self, node_type_to_smiles: Dict[str, str], node_type_port_to_index: Dict[str, Dict[any, int]]) -> str:
        """
        Convert the GroupGraph to a SMILES notation.

        Parameters:
        - node_type_to_smiles (Dict[str, str]): Mapping from node type to SMILES notation.
        - node_type_port_to_index (Dict[str, Dict[any, int]]): Mapping from node type and port to index.

        Returns:
        - str: SMILES notation.
        """
        molecular_graph = self.to_molecular_graph(node_type_to_smiles, node_type_port_to_index)
        return write_smiles(molecular_graph)
    
    def to_PyG_Data(self, node_descriptor_generater: Callable[[str], Sequence[float]]) -> torch_geometric.data.Data:
        """
        Convert the GroupGraph to a PyG Data representation.

        Parameters:
        - node_descriptor_generater (Callable[[str], Sequence[float]]): Callable to generate node descriptors.

        Returns:
        - torch_geometric.data.Data: PyG Data representation.
        """
        from torch_geometric.data import Data
        import torch
        max_n_attachments = max(len(v) for v in self.node_types.values())
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
    
            
