import networkx as nx
import rdkit
import rdkit.Chem
import group_selfies
#   taken from a good answer by Gambit1614 on StackOverflow https://stackoverflow.com/questions/57095809/networkx-connecting-nodes-using-ports

class GroupGraph(nx.Graph):
    """A graph with ports as parts of nodes that can be connected to other ports."""

    def __init__(self, node_types=None):
        super(GroupGraph, self).__init__()
        self.node_types = node_types

    def __str__(self):
        return f"Nodes: {self.nodes.data('ports')} \nEdges: {self.edges.data('ports')}"

    def __repr__(self):
        return f"Nodes: {self.nodes.data('type')}, {self.nodes.data('ports')} \nEdges: {self.edges.data('ports')}"
    
    def __eq__(self, other):
        return self.nodes == other.nodes and self.edges == other.edges
    
    def __ne__(self, other):
        return not self.__eq__(other)

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
    
    def from_compound(self, compound, groups):
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
            
