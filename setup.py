from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import sys
import os
import subprocess
import pybind11

class CUDAExtension(Extension):
    """Extension class for CUDA sources."""
    def __init__(self, name, sources, **kwargs):
        super().__init__(name, sources, **kwargs)

class CustomBuildExt(build_ext):
    """Custom build_ext class to handle CUDA files and macro definitions."""
    def build_extensions(self):
        # Retrieve Conda prefix
        condabase = os.environ.get('CONDA_PREFIX')
        if not condabase:
            raise EnvironmentError("CONDA_PREFIX environment variable is not set. Activate your Conda environment.")

        # Define CUDA home
        cuda_home = os.environ.get('CUDA_HOME', '/usr/local/cuda')
        if not os.path.exists(cuda_home):
            raise EnvironmentError(f"CUDA_HOME is not set or CUDA not found at {cuda_home}")

        # Path to nvcc
        nvcc = os.path.join(cuda_home, 'bin', 'nvcc')
        if not os.path.isfile(nvcc):
            raise EnvironmentError(f"nvcc not found at {nvcc}")

        for ext in self.extensions:
            print(f"Processing extension: {ext.name}")
            print(f"Original sources: {ext.sources}")

            # Separate CUDA and C++ sources
            cuda_sources = [s for s in ext.sources if s.endswith('.cu')]
            cpp_sources = [s for s in ext.sources if not s.endswith('.cu')]

            print(f"CUDA sources: {cuda_sources}")
            print(f"C++ sources: {cpp_sources}")

            # Compile CUDA sources with nvcc
            if cuda_sources:
                for source in cuda_sources:
                    # Determine the object file path
                    obj = os.path.splitext(source)[0] + '.o'
                    obj_path = os.path.join(self.build_temp, obj)

                    # Ensure the directory exists
                    obj_dir = os.path.dirname(obj_path)
                    if not os.path.exists(obj_dir):
                        os.makedirs(obj_dir)

                    # Compile with nvcc
                    cmd = [
                        nvcc,
                        '-O2',
                        '-std=c++17',
                        '-c', source,
                        '-o', obj_path,
                        '-Xcompiler', '-fPIC',
                        '-D_FILE_OFFSET_BITS=64',  # Override the macro definition
                        '-I', os.path.join(cuda_home, 'include'),
                    ]

                    # Add include directories
                    for inc_dir in ext.include_dirs:
                        cmd += ['-I', inc_dir]

                    # Debug: Print the compilation command
                    print(f"Compiling CUDA source: {source}")
                    print(f"Compilation command: {' '.join(cmd)}")

                    # Execute the compilation
                    subprocess.check_call(cmd)

                # Remove .cu sources from ext.sources
                ext.sources = [s for s in ext.sources if not s.endswith('.cu')]
                print(f"Updated sources after removing .cu files: {ext.sources}")

                # Add object files to extra_objects
                object_files = [os.path.splitext(s)[0] + '.o' for s in cuda_sources]
                object_files_full = [os.path.join(self.build_temp, obj) for obj in object_files]
                ext.extra_objects = getattr(ext, 'extra_objects', []) + object_files_full
                print(f"Added extra objects: {object_files_full}")

                # Add CUDA library directories and libraries
                ext.library_dirs.append(os.path.join(cuda_home, 'lib64'))
                ext.libraries.append('cudart')
                print(f"Updated library_dirs: {ext.library_dirs}")
                print(f"Updated libraries: {ext.libraries}")

                # Embed RPATH to find CUDA libraries at runtime
                if not hasattr(ext, 'extra_link_args') or ext.extra_link_args is None:
                    ext.extra_link_args = []
                ext.extra_link_args += [f'-Wl,-rpath,{os.path.join(cuda_home, "lib64")}']
                print(f"Updated extra_link_args: {ext.extra_link_args}")

            # Ensure all C++ sources are compiled with -fPIC and _FILE_OFFSET_BITS=64
            if not hasattr(ext, 'extra_compile_args') or ext.extra_compile_args is None:
                ext.extra_compile_args = []
            ext.extra_compile_args += ['-fPIC', '-D_FILE_OFFSET_BITS=64', '-std=c++17']
            print(f"Updated extra_compile_args: {ext.extra_compile_args}")

        # Proceed with the standard build
        build_ext.build_extensions(self)

# Retrieve Conda prefix
condabase = os.environ.get('CONDA_PREFIX')
if not condabase:
    raise EnvironmentError("CONDA_PREFIX environment variable is not set. Activate your Conda environment.")

# Define CUDA home
cuda_home = os.environ.get('CUDA_HOME', '/usr/local/cuda')
if not os.path.exists(cuda_home):
    raise EnvironmentError(f"CUDA_HOME is not set or CUDA not found at {cuda_home}")

# Path to CUB
cub_include = os.path.join(cuda_home, 'include', 'cub')
if not os.path.exists(cub_include):
    # If CUB is not part of CUDA installation, you might need to install it manually
    # For example, download from https://github.com/NVIDIA/cub and place it accordingly
    raise EnvironmentError(f"CUB include directory not found at {cub_include}")

# Define the extension module with CUDA sources
molgrouper_module = CUDAExtension(
    'Grouper._Grouper',
    sources=[
        'Grouper/dataStructures.cpp',
        'Grouper/binding.cpp',
        'Grouper/processColoredGraphs.cpp',
        'Grouper/generate.cpp',
        'Grouper/fragmentation.cpp',
        'Grouper/autUtils.cpp',
        'Grouper/colorPermutations.cu',  # Add your CUDA source files here
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
        cub_include,
        os.path.join(cuda_home, 'include'),
    ],
    library_dirs=[
        os.path.join(condabase, 'lib'),
        os.path.join(condabase, "lib/cairo"),
        os.path.join(cuda_home, 'lib64'),
    ],
    libraries=[
        'RDKitFileParsers', 
        'RDKitSmilesParse', 
        'RDKitGraphMol', 
        'RDKitRDGeneral', 
        'omp', 
        'nauty', 
        'pq',
        'cudart',  # CUDA runtime library
    ],
    language='c++',
    extra_link_args=[
        f'-Wl,-rpath,{os.path.join(condabase, "lib")}',
        f'-Wl,-rpath,{os.path.join(cuda_home, "lib64")}',
    ],
)

setup(
    name='Grouper',
    version='0.1',
    packages=find_packages(),
    author='Kieran Nehil-Puleo',
    description='A package for working with group graphs',
    ext_modules=[molgrouper_module],
    cmdclass={
        'build_ext': CustomBuildExt
    },
    zip_safe=False,
)
