import typing as t
import ctypes
import pathlib

_libname = str(pathlib.Path(__file__).parent) + "/../packages/nauty2_8_8/libgeng.so"
libgeng = ctypes.CDLL(_libname) # Load the shared library into c types.
libgeng.main.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_char_p)] # for now just use the main function of geng


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

def all_possible(n_nodes, max_edges):
    _call_geng(n_nodes, max_edges)

def subset_possible():
    pass

def apply_filters():
    pass

def _node_sub():
    pass

def _edge_sub():
    pass

def _call_geng(n_nodes, max_edges):
    args = ["geng", str(n_nodes), f"1:{max_edges}" , "-ctf"]
    argv = [ctypes.c_char_p(arg.encode('utf-8')) for arg in args]
    libgeng.main(len(args), (ctypes.c_char_p * len(argv))(*argv))