#!/bin/bash
#SBATCH --mail-user=kieran.d.nehil-puleo@vanderbilt.edu
#SBATCH --mail-type=end
#SBATCH --error=error-%J.err
#SBATCH --output=output-%J.out
#SBATCH --job-name=genGrouper
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --partition=week-long-std
#SBATCH --mem=31G
#SBATCH --cpus-per-task=32


module load anaconda/3.9
conda activate genGrouper
echo $CONDA_PREFIX
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH


echo "--> Activated conda environment"
# -------------------------------------------------------------------------------   
cd /raid6/homes/kierannp/foo/grouper_performance/with_orbit/genGrouper
python setup.py clean --all
python setup.py build_ext --inplace
python setup.py install --user
echo "--> Compiled genGrouper"

echo "--> Running performance testing of exhaustive generate with orbit"
pytest genGrouper/tests/test_generate.py --benchmark-save=with_orbit
# --------------------------------------------------------------------------------
# cd /raid6/homes/kierannp/foo/grouper_performance/with_auto/genGrouper
# python setup.py clean --all
# python setup.py build_ext --inplace
# python setup.py install --user
# echo "--> Compiled genGrouper"

# echo "--> Running performance testing of exhaustive generate with just auto"
# pytest genGrouper/tests/test_generate.py --benchmark-save=with_auto
# -------------------------------------------------------------------------------
# cd /raid6/homes/kierannp/foo/grouper_performance/without_either/genGrouper
# python setup.py clean --all
# python setup.py build_ext --inplace
# python setup.py install --user
# echo "--> Compiled genGrouper"

# echo "--> Running performance testing of exhaustive generate without either"
# pytest genGrouper/tests/test_generate.py --benchmark-save=without_either
# -------------------------------------------------------


echo "--> Done"
