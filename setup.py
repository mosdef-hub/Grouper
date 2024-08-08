from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import pybind11

condabase = os.environ['CONDA_PREFIX']

# Define the extension module
molgrouper_module = Extension(
    'molGrouper._molGrouper',
    sources=[
        'molGrouper/dataStructures.cpp',
        'molGrouper/binding.cpp',
        'molGrouper/processColoredGraphs.cpp',
        'molGrouper/generate.cpp',
    ],
    include_dirs = [
            os.path.join(condabase, 'include'),
            os.path.join(condabase, "include/cairo"),
            os.path.join(condabase, "include/boost"),
            os.path.join(condabase, "include/rdkit"),
            os.path.join(condabase, "include/omp"),
            pybind11.get_include(),
            pybind11.get_include(user=True),
            'molGrouper',
    ],
    library_dirs = [
            os.path.join(condabase, 'lib'),
            os.path.join(condabase, "lib/cairo")
    ],
    libraries = ['RDKitFileParsers', 'RDKitSmilesParse', 'RDKitGraphMol', 'RDKitRDGeneral', 'omp'],
    extra_compile_args = ['-Xpreprocessor', '-fopenmp', '-std=c++17', '-mmacosx-version-min=10.13'],
    language='c++',
    extra_link_args=['-Wl'])

setup(
    name='molGrouper',
    version='0.1',
    packages=find_packages(),
    author='Kieran Nehil-Puleo',
    description='A Python package with C++ extensions',
    ext_modules=[molgrouper_module],
    cmdclass={
        build_ext: build_ext
    }
)
