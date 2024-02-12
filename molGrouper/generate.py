import typing as t
import ctypes
import pathlib
import os
import subprocess
from molGrouper.group_graph import GroupGraph

_libgeng = str(pathlib.Path(__file__).parent) + "/../packages/nauty2_8_8/libgeng.so"
libgeng = ctypes.CDLL(_libgeng) # Load the shared library into c types.
libgeng.main.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)] # for now just use the main function of geng
libgeng.main.restype = None

_libvcolg = str(pathlib.Path(__file__).parent) + "/../packages/nauty2_8_8/libvcolg.so"
libvcolg = ctypes.CDLL(_libvcolg) # Load the shared library into c types.
libvcolg.main.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)] # for now just use the main function of geng
libvcolg.main.restype = None


_libmultig = str(pathlib.Path(__file__).parent) + "/../packages/nauty2_8_8/libmultig.so"
libmultig = ctypes.CDLL(_libmultig) # Load the shared library into c types.
libmultig.main.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)] # for now just use the main function of geng


"""
To create libgeng.so:
    1. Download nauty.tar.gz from http://pallini.di.uniroma1.it/
    2. Extract the tar file and cd into the directory
    3. ./configure
    4. make
    5. gcc -o geng.o -c -O3  -march=native -DMAXN=WORDSIZE -DWORDSIZE=32 geng.c gtoolsW.o nautyW1.o \
                nautilW1.o naugraphW1.o schreier.o naurng.o
    6. gcc -shared -O3  -march=native -DMAXN=WORDSIZE -DWORDSIZE=32 geng.o gtoolsW.o nautyW1.o \
                nautilW1.o naugraphW1.o schreier.o naurng.o -o libgeng.so

"""


def generate_group_graph_space(n_nodes, max_n_attachments, node_types):
    unfiltered_space = _all_possible_group_graphs(n_nodes, max_n_attachments, len(node_types))
    return _process_multig_output(node_types, max_n_attachments)
    # filtered_space = _apply_filters(unfiltered_space)
    # return _convert_space_to_group_graphs(filtered_space)

def _all_possible_group_graphs(n_nodes, max_n_attachments, n_node_types):
    _call_geng(n_nodes, max_n_attachments)
    _call_vcolg(n_node_types)
    _call_multig(max_n_attachments**2)

def subset_possible():
    pass

def _apply_filters():
    pass

def _call_geng(n_nodes, max_edges):
    geng_path = str(pathlib.Path(__file__).parent) + "/../packages/nauty2_8_8/geng"
    args = [geng_path, str(n_nodes), f"1:{max_edges}", "geng_out.txt", "-ctf"]
    subprocess.run(args)
    if not os.path.exists("geng_out.txt"):
        raise FileNotFoundError("geng failed to create all possible graphs. Check the input parameters.")

def _call_vcolg(n_colors):
    vcolg_path = str(pathlib.Path(__file__).parent) + "/../packages/nauty2_8_8/vcolg"
    args = [vcolg_path, "geng_out.txt", "-T", f"-m{n_colors}", "vcolg_out.txt"]
    subprocess.run(args)
    if not os.path.exists("vcolg_out.txt"):
        raise FileNotFoundError("vcolg failed to create all graphs with colors. Check the input parameters.")

def _call_multig(edge_multiplicity):
    multig_path = str(pathlib.Path(__file__).parent) + "/../packages/nauty2_8_8/multig"
    args = [multig_path, "vcolg_out.txt", "-T", "-V", f"-m{edge_multiplicity}", "multig_out.txt"]
    subprocess.run(args)
    if not os.path.exists("multig_out.txt"):
        raise FileNotFoundError("multig failed to create all graphs with different edges. Check the input parameters.")

def _multi_to_pair(multi, max_multi):
    if 1 <= multi <= max_multi:
        # Calculate x and y values based on the input multi
        x = (multi - 1) % int(max_multi**.5) + 1
        y = (multi - 1) // int(max_multi**.5) + 1
        return x, y
    else:
        raise ValueError("Input multi must be in the range 1-16")

def _multig_output_to_graphs(line, node_types, max_n_attachments):
    node_description, edge_description = line.split("  ")
    edge_description = edge_description.split(" ")

    # edges are in the form of (node1, node2, edge_type)
    edge_list = []
    for i in range(0, len(edge_description), 3):
        edge_list.append((int(edge_description[i]), int(edge_description[i+1]), int(edge_description[i+2]))) # source, destination, edge_type

    # graph is in the form of (node_type, edges {color})
    node_description = node_description.split(" ")
    n_vertices = int(node_description[0])
    n_edges = int(node_description[1])
    colors = node_description[2:]
    
    gG = GroupGraph(node_types)
    # Add nodes
    for i in range(n_vertices):
        gG.add_node(f"node{i}", colors[i])
    # Add edges
    try:
        for e in edge_list:
            port1, port2 = _multi_to_pair(e[2], max_n_attachments**2)
            gG.add_edge(f'node{e[0]}', f'port{port1}', f'node{e[1]}', f'port{port2}')
        return gG
    except:
        print("Couldn't produce graph from multig output")

def _process_multig_output(node_types, max_n_attachments):
    for line in open("multig_out.txt"):
        yield _multig_output_to_graphs(line, node_types, max_n_attachments)
    os.remove("multig_out.txt")
    os.remove("vcolg_out.txt")
    os.remove("geng_out.txt")


if __name__ == "__main__":
    node_types = {
        '0': ['port1', 'port2'],
        '1': ['port1', 'port2'],
        '2': ['port1', 'port2', 'port3'],
    }
    print(list(generate_group_graph_space(3, 3, node_types)))