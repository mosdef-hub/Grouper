from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import pybind11

rdbase = os.environ['RDBASE']

# Define the extension module
molgrouper_module = Extension(
    'molGrouper',
    sources=[
        'molGrouper/dataStructures.cpp',
        'molGrouper/binding.cpp',
    ],
    include_dirs = [
            os.path.join(rdbase, 'Code'),
            os.path.join(rdbase, 'lib'),
            '/opt/homebrew/Cellar/cairo/1.18.0/include/cairo',
            '/opt/homebrew/Cellar/boost/1.85.0/include',
            '/opt/homebrew/Cellar/libomp/18.1.8/include',
            pybind11.get_include(),
            pybind11.get_include(user=True),
            'molGrouper',
    ],
    library_dirs = [
            os.path.join(rdbase, 'lib'),
            '/opt/homebrew/Cellar/libomp/18.1.8/lib',
            '/opt/homebrew/Cellar/boost/1.85.0/lib',
    ],
    libraries = ['RDKitFileParsers', 'RDKitSmilesParse', 'RDKitGraphMol', 'RDKitRDGeneral', 'omp'],
    extra_compile_args = ['-Xpreprocessor', '-fopenmp', '-std=c++17'],
    language='c++',
    extra_link_args=['-Wl', '-rpath', os.path.join(rdbase, 'lib')],
)

setup(
    name='molGrouper',
    version='0.1',
    packages=find_packages('molGrouper'),
    author='Kieran Nehil-Puleo',
    description='A Python package with C++ extensions',
    ext_modules=[molgrouper_module],
    cmdclass={
        build_ext: build_ext
    }
)