from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import pybind11

condabase = os.environ['CONDA_PREFIX']

# Platform-specific compiler and linker flags
platform_compile_args = ['-std=c++17', '-fPIC', '-Xpreprocessor', '-fopenmp', '-g', '-O3']
platform_link_args = []
if sys.platform == 'darwin':  # macOS
    platform_compile_args += ['-arch', 'arm64']  # Ensure architecture compatibility
elif sys.platform.startswith('linux'):  # Linux
    platform_compile_args.append('-D_Thread_local=thread_local')

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
    include_dirs=[
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
    library_dirs=[
        os.path.join(condabase, 'lib'),
        os.path.join(condabase, "lib/cairo"),
    ],
    libraries=[
        'RDKitFileParsers', 
        'RDKitSmilesParse', 
        'RDKitGraphMol', 
        'RDKitRDGeneral', 
        'omp', 
        'nauty', 
        'pq',
    ],
    extra_compile_args=platform_compile_args,
    extra_link_args=platform_link_args,
    language='c++',
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
