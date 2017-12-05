#!/bin/bash
#SBATCH --exclusive
####SBATCH -p bnl
#SBATCH --nodes=1
#SBATCH -o sppipe.%N.out -e sppipe.%N.err
#SBATCH -t 36:00:00 
###SBATCH --reservation=andras_test

module load gcc
module load openmpi
module load lib/gsl/2.3
module load python3/3.6.2

export PYTHONPATH=$PYTHONPATH:/cm/shared/sw/pkg-old/devel/python3/3.6.2/lib/python3.6/site-packages

#srun python3 PowerSpectraPipeline.py 200 50 0 4 12 14
srun python3 PowerSpectraPipeline.py 20 100 >> out.txt
srun python3 PowerSpectraPipeline.py 100 20 >> out.txt
srun python3 PowerSpectraPipeline.py 200 10 >> out.txt
srun python3 PowerSpectraPipeline.py 400 5 >> out.txt
srun python3 PowerSpectraPipeline.py 500 4 >> out.txt



