import os 

#directory = '/mnt/home/landerson/lyalpha/'
#directory = './'
directory = '/home/landerson/lyalpha/40Mpc_512/'
for j in range(50):
    try: os.mkdir(directory + str(j))
    except FileExistsError: print ('file exists')
    try: os.mkdir(directory + 'NCV_0_{0}'.format(j))
    except FileExistsError: print ('file exists')
    try: os.mkdir(directory + 'NCV_1_{0}'.format(j))
    except FileExistsError: print ('file exists')

    k = j 
    sbatchStringStart = """#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -o psV{0}.out -e psV{1}.err
#SBATCH -J ps{2}
#SBATCH --mail-user=landerson@flatironinstitute.org     
#SBATCH --mail-type=ALL  
#SBATCH -t 7-00:00
#SBATCH --partition=general
#SBATCH --export=ALL

module load intel
module load openmpi2_ib

export OMP_NUM_THREADS=1
date
""".format(k, k, k)

    with open(directory + 'psVar{0}.sbatch'.format(k), 'w') as f:
        f.write(sbatchStringStart)
        #for i in range(5):
        num = k 
        #pspipeLine = "srun python3 /mnt/ceph/users/landerson/LyA-InvertPhase/PowerSpectraPipeline_single.py 200 10 lyalphaVaried{0} || exit 1 \n".format(num)
        pspipeLineT = "srun python3 /home/landerson/src/LyA-InvertPhase/PowerSpectraPipeline_single.py 400 10 {0} || exit 1 \n".format(num)
        pspipeLineP1 = "srun python3 /home/landerson/src/LyA-InvertPhase/PowerSpectraPipeline_single.py 400 10 NCV_0_{0} || exit 1 \n".format(num)
        pspipeLineP2 = "srun python3 /home/landerson/src/LyA-InvertPhase/PowerSpectraPipeline_single.py 400 10 NCV_1_{0} || exit 1 \n".format(num)

        f.write(pspipeLineT)
        f.write(pspipeLineP1)
        f.write(pspipeLineP2)
        f.write('date \n')


with open(directory + 'submit_psjobs.sbatch', 'w') as f:
    for j in range(50):
        f.write('sbatch psVar{0}.sbatch \n'.format(j))



with open(directory + 'makePK.sbatch', 'w') as f:
    startscript = """#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1 
#SBATCH -o pk.out -e pk.err 
#SBATCH -J ps{2} 
#SBATCH --mail-user=landerson@flatironinstitute.org
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00 
#SBATCH --partition=general
#SBATCH --export=ALL

module load intel
module load hdf5

"""
    f.write(startscript)
    for j in range(50):
        snapshot_directory_pre = '/home/fvillaescusa/data/Lya_ncv/40Mpc_512/'
        snapshot_dirs = ['{0}/'.format(j), 'NCV_0_{0}/'.format(j), 'NCV_1_{0}/'.format(j)]
        snapshot_save_directory_pre = '/home/landerson/lyalpha/40Mpc_512/'
        snapshots = [0, 1, 2]
        for sdir in snapshot_dirs:
            for snap in snapshots:
                cmd = 'srun /home/landerson/src/GenPK/gen-pk -i {0}snap_{1:03d} -o {2} \n'.format(snapshot_directory_pre+sdir, snap, snapshot_save_directory_pre+sdir)
                f.write(cmd)
