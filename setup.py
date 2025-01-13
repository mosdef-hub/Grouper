from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import pybind11

condabase = os.environ['CONDA_PREFIX']

# Define the extension module
molgrouper_module = Extension(
    'Grouper._Grouper',
    sources=[
        'Grouper/dataStructures.cpp',
        'Grouper/binding.cpp',
        'Grouper/processColoredGraphs.cpp',
        'Grouper/generate.cpp',
        'Grouper/fragmentation.cpp',
        'Grouper/autUtils.cpp',
    ],
    include_dirs = [
        os.path.join(condabase, 'include'),
        os.path.join(condabase, "include/cairo"),
        os.path.join(condabase, "include/boost"),
        os.path.join(condabase, "include/rdkit"),
        os.path.join(condabase, "include/omp"),
        os.path.join(condabase, "include/nauty"),
        os.path.join(condabase, "include/libpq"),
        pybind11.get_include(),
        pybind11.get_include(user=True),
        'Grouper',
    ],
    library_dirs = [
        os.path.join(condabase, 'lib'),
        os.path.join(condabase, "lib/cairo"),
    ],
    libraries = [
        'RDKitFileParsers', 
        'RDKitSmilesParse', 
        'RDKitGraphMol', 
        'RDKitRDGeneral', 
        'omp', 
        'nauty', 
        'pq'
    ],
    extra_compile_args = [
        '-Xpreprocessor', '-fopenmp', '-std=c++17', '-g',
        '-fPIC',  # Ensure proper position-independent code for shared libraries
    ],
    language='c++',
    # Optionally add separate flags for C compilation (for nauty or other C files)
    undef_macros=['_Thread_local'],  # This will undefine any conflicting macro definitions if needed
    extra_link_args=[],
)

setup(
    name='Grouper',
    version='0.1',
    packages=find_packages(),
    author='Kieran Nehil-Puleo',
    description='A package for working with group graphs',
    ext_modules=[molgrouper_module],
    cmdclass={
        build_ext: build_ext
    }
)
