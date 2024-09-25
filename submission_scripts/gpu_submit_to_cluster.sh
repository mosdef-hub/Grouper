#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=error-%J.err
#SBATCH --output=output-%J.out
#SBATCH --job-name=gpu_gen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:A100:1
#SBATCH --partition=day-long-tesla
#SBATCH --mem=9G  # Example to request 16GB of memory



source /raid6/homes/kierannp/.bashrc
module load anaconda/3.9
conda activate genGrouper

cd /raid6/homes/kierannp/projects/genGrouper 

python setup.py clean --all
python setup.py build_ext --inplace
python setup.py install

pytest genGrouper/tests --benchmark-save=gpu_run_submit
