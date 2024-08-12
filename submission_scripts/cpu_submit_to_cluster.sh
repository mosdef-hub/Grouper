#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=error-%J.err
#SBATCH --output=output-%J.out
#SBATCH --job-name=genGrouper
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=day-long-std
#SBATCH --mem=31G



source /raid6/homes/kierannp/.bashrc
module load anaconda/3.9
conda activate pureGrouper

RUN_DIR=

python run_exhaustive_generation.py
