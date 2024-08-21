from genGrouper.utils import run_performance_eval
from genGrouper import Node


node_defs = set()
node_defs.add(Node(0, 'benzene', 'C1=CC=CC=C1', [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]))
node_defs.add(Node(0, 'hydroxyl', 'O', [0], [0]))
node_defs.add(Node(0, 'secondary_amine', 'N', [0,1], [0,0]))
node_defs.add(Node(0, 'tertiary_amine', 'N', [0,1,2], [0,0,0]))
node_defs.add(Node(0, 'alkene', 'C=C', [0], [0]))
node_defs.add(Node(0, 'ester', 'C(=O)O', [0,1], [0,2]))

performance = run_performance_eval(
    nauty_path="/raid6/homes/kierannp/projects/molGrouper/packages/nauty2_8_8",
    node_defs=node_defs,
    n_runs=3,
    max_nodes=5,
    n_cpus=30,
    verbose=False
)
print(performance)
