#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=error-%J.err
#SBATCH --output=output-%J.out
#SBATCH --job-name=cpu_gen
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=day-long-std
#SBATCH --mem=31G



source /raid6/homes/kierannp/.bashrc
module load anaconda/3.9
conda activate molGPU

cd /raid6/homes/kierannp/foo/gpufoo/molGrouper/molGrouper/cpp_code/rdkit_cpu

rm -rf build
mkdir build
cd build
cmake ..
make
./lineProcessor

