from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import pybind11

condabase = os.environ['CONDA_PREFIX']

# Define the extension module
molgrouper_module = Extension(
    'genGrouper._genGrouper',
    sources=[
        'genGrouper/dataStructures.cpp',
        'genGrouper/binding.cpp',
        'genGrouper/processColoredGraphs.cpp',
        'genGrouper/generate.cpp',
        'genGrouper/fragmentation.cpp',
        'genGrouper/autUtils.cpp',
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
            'genGrouper',
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
    extra_compile_args = ['-Xpreprocessor', '-fopenmp', '-std=c++17', '-g'], # '-mmacosx-version-min=10.13'
    language='c++')

setup(
    name='genGrouper',
    version='0.1',
    packages=find_packages(),
    author='Kieran Nehil-Puleo',
    description='A package for working with group graphs',
    ext_modules=[molgrouper_module],
    cmdclass={
        build_ext: build_ext
    }
)