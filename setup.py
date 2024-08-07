from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import sys
import pybind11

# class BuildExt(build_ext):
#     def build_extensions(self):
#         # Ensure RDBASE is set
#         if 'RDBASE' not in os.environ:
#             raise RuntimeError("RDBASE environment variable is not set")

#         rdbase = os.environ['RDBASE']
#         # Set the include and library directories
#         include_dirs = [
#             os.path.join(rdbase, 'Code'),
#             '/opt/homebrew/Cellar/cairo/1.18.0/include/cairo',
#             '/opt/homebrew/Cellar/boost/1.85.0/include/boost',
#             '/opt/homebrew/Cellar/libomp/18.1.8/include',
#             'molGrouper',
#         ]
#         library_dirs = [
#             os.path.join(rdbase, 'lib'),
#             '/opt/homebrew/Cellar/libomp/18.1.8/lib',
#             '/opt/homebrew/Cellar/boost/1.85.0/lib',
#         ]
#         libraries = [
#             'RDKitFileParsers',
#             'RDKitSmilesParse',
#             'RDKitGraphMol',
#             'RDKitRDGeneral',
#             'omp',
#         ]

#         # Set extra compile args
#         extra_compile_args = ['-Xpreprocessor', '-fopenmp', '-std=c++17'] # remove '-Xpreprocessor',

#         # # Set the library paths
#         # extra_link_args = [
#         #     os.path.join(rdbase, 'lib', 'libRDKitFileParsers.dylib'),
#         #     os.path.join(rdbase, 'lib', 'libRDKitSmilesParse.dylib'),
#         #     os.path.join(rdbase, 'lib', 'libRDKitGraphMol.dylib'),
#         #     os.path.join(rdbase, 'lib', 'libRDKitRDGeneral.dylib'),
#         #     '/opt/homebrew/Cellar/libomp/18.1.8/lib/libomp.dylib',
#         # ]

#         for ext in self.extensions:
#             ext.include_dirs = include_dirs
#             ext.library_dirs = library_dirs
#             ext.libraries = libraries
#             ext.extra_compile_args = extra_compile_args
#             # ext.extra_link_args = extra_link_args

#         build_ext.build_extensions(self)

#rdbase = os.environ['RDBASE']
condabase = os.environ['CONDA_PREFIX']

# Define the extension module
molgrouper_module = Extension(
    'molGrouper',
    sources=[
        'molGrouper/dataStructures.cpp',
        'molGrouper/binding.cpp',
    ],
    include_dirs = [
            os.path.join(condabase, 'include'),
            os.path.join(condabase, "include/cairo"),
            os.path.join(condabase, "include/boost"),
            os.path.join(condabase, "include/rdkit"),
            #os.path.join(condabase, "include/omp")
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
    packages=find_packages('molGrouper'),
    author='Kieran Nehil-Puleo',
    description='A Python package with C++ extensions',
    ext_modules=[molgrouper_module],
    cmdclass={
        build_ext: build_ext
    }
)
