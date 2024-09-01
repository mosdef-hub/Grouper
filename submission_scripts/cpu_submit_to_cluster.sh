#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=error-%J.err
#SBATCH --output=output-%J.out
#SBATCH --job-name=genGrouper
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --partition=week-long-std
#SBATCH --mem=31G


RUN_DIR=/raid6/homes/kierannp/projects/genGrouper
cd RUN_DIR

source /raid6/homes/kierannp/.bashrc
module load anaconda/3.9
conda activate genGrouper

python setup.py build_ext --inplace
python setup.py install

python run_exhaustive_generate.py --n 7 --n_cpus 30
