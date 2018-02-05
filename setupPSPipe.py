directory = '/mnt/home/landerson/lyalpha/'
#directory = './'
for j in range(46):
    k = j + 5
    sbatchStringStart = """#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -o psV{0}.out -e psV{1}.err
#SBATCH --reservation andras_test
#SBATCH -p cca
#SBATCH --qos cca

module load gcc
module load openmpi
module load lib/gsl/2.3
module load python3

export OMP_NUM_THREADS=1
date
""".format(k, k)

    with open(directory + 'psVar{0}.sbatch'.format(k), 'w') as f:
        f.write(sbatchStringStart)
        #for i in range(5):
        num = k 
        pspipeLine = "srun python3 /mnt/ceph/users/landerson/LyA-InvertPhase/PowerSpectraPipeline_single.py 200 10 lyalphaVaried{0} || exit 1 \n".format(num)
        f.write(pspipeLine)
        f.write('date \n')


with open(directory + 'submit_psjobs.sbatch', 'w') as f:
    for j in range(46):
        f.write('sbatch psVar{0}.sbatch \n'.format(j+5))

