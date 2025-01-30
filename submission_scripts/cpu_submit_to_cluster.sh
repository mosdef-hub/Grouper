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




module load anaconda/3.9
conda activate genGrouper
echo $CONDA_PREFIX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

echo "--> Activated conda environment"

cd /raid6/homes/kierannp/projects/genGrouper
python setup.py clean --all
python setup.py build_ext --inplace
python setup.py install --user

echo "--> Compiled genGrouper"

python run_exhaustive_generate.py --n 7 --n_cpus 30

