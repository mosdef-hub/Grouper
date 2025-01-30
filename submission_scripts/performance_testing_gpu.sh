#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=error-%J.err
#SBATCH --output=output-%J.out
#SBATCH --job-name=gpu_gen
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:A100:1
#SBATCH --partition=day-long-tesla
#SBATCH --mem=9G  # Example to request 16GB of memory



module load anaconda/3.9
conda activate genGrouper
echo $CONDA_PREFIX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

echo "--> Activated conda environment"
# -------------------------------------------------------------------------------   
cd /raid6/homes/kierannp/foo/grouper_performance/with_gpu/genGrouper
python setup.py clean --all
python setup.py build_ext --inplace
python setup.py install --user
echo "--> Compiled genGrouper"

echo "--> Running performance testing of exhaustive generate gpu and auto"
pytest genGrouper/tests/test_generate.py --benchmark-save=gpu_auto
# --------------------------------------------------------------------------------

echo "--> Done"