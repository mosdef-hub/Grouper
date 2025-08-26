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



source /raid6/homes/kierannp/.bashrc
load
conda activate genGrouper
# micromamba activate /raid6/homes/kierannp/y/envs/pureGrouper

echo "--> Activated conda environment"

cd /raid6/homes/kierannp/projects/genGrouper
python setup.py clean --all
python setup.py build_ext --inplace
python setup.py install

echo "--> Compiled genGrouper"

python run_exhaustive_generate.py --n 7 --n_cpus 30 --config_path /raid6/homes/kierannp/projects/genGrouper/dfconfig.cfg
