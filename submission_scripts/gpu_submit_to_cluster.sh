#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=error-%J.err
#SBATCH --output=output-%J.out
#SBATCH --job-name=gpu_gen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=day-long-tesla
#SBATCH --mem=9G  # Example to request 16GB of memory



source /raid6/homes/kierannp/.bashrc
module load anaconda/3.9
conda activate molGPU

cd /raid6/homes/kierannp/foo/gpufoo/genGrouper/genGrouper/cpp_code/rdkit_gpu

rm -rf build
mkdir build
cd build
cmake ..
make
./lineProcessor
