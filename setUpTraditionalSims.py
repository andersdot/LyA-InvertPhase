import numpy as np

directory = '/mnt/home/landerson/lyalpha/'

for m in range(46):
    n = m + 5
    seed = np.random.randint(100000)

    #make paramfile for genic
    genicParamString = """OutputDir = output # Directory for output
FileBase = IC_varied_{0}  # Base-filename of output files

InvertPhase = 0
UnitaryAmplitude = 0

Ngrid = 256 # Size of cubic grid on which to create particles.

BoxSize = 20000   # Periodic box size of simulation



Omega0 = 0.2814      # Total matter density  (at z=0)
OmegaLambda = 0.7186      # Cosmological constant (at z=0)
OmegaBaryon = 0.0464     # Baryon density        (at z=0)
ProduceGas = 1         # 1 = Produce gas  0 = no gas, just DM.
HubbleParam = 0.697      # Hubble paramater (may be used for power spec parameterization)

Redshift = 99        # Starting redshift

Sigma8 = 0.810      # power spectrum normalization

WhichSpectrum = 2         # "1" selects Eisenstein & Hu spectrum,
		                   # "2" selects a tabulated power spectrum in
                           # the file 'FileWithInputSpectrum'
                           # otherwise, Efstathiou parametrization is used


FileWithInputSpectrum = /mnt/xfs1/home/landerson/src/MP-Gadget/examples/powerspectrum-wmap9.txt  # filename of tabulated input
                                                                   # spectrum (if used)
InputSpectrum_UnitLength_in_cm  = 3.085678e24 # defines length unit of tabulated
                                                                   # input spectrum in cm/h.
                                                # Note: This can be chosen different from UnitLength_in_cm


PrimordialIndex  0.971 = # may be used to tilt the primordial index

Seed = {1}    #  seed for IC-generator


UnitLength_in_cm = 3.085678e21   # defines length unit of output (in cm/h)
UnitMass_in_g = 1.989e43      # defines mass unit of output (in g/cm)
UnitVelocity_in_cm_per_s = 1e5 # defines velocity unit of output (in cm/sec)
""".format(n, seed)

    with open(directory + 'paramfileVaried{0}.genic'.format(n), 'w') as f:
        f.write(genicParamString)

    #make paramfile for gadget

    gadgetParamString = """#  Relevant files

InitCondFile = output/IC_varied_{0}
OutputDir =  /mnt/cephtest/landerson/lyalphaVaried{1}
#OutputDir = /mnt/ceph/users/landerson/lyalphaVaried1
TreeCoolFile = /mnt/xfs1/home/landerson/src/MP-Gadget/examples/TREECOOL_fg_june11
OutputList = 0.010101010101010102,0.02,0.1,0.192307692308,0.2,0.208333333333,0.217391304348,0.227272727273,0.238095238095,0.25,0.263157894737,0.277777777778,0.294117647059,0.3125


Nmesh = 512

# CPU time -limit

TimeLimitCPU = 43000 #= 8 hours

# Code options

#  Characteristics of run

TimeMax = 0.33333

Omega0 = 0.2814      # Total matter density  (at z=0)
OmegaLambda = 0.7186      # Cosmological constant (at z=0)
OmegaBaryon = 0.0464     # Baryon density        (at z=0)
HubbleParam = 0.697      # Hubble paramater (may be used for power spec parameterization)

CoolingOn = 1
StarformationOn = 1
RadiationOn = 1
HydroOn = 1
BlackHoleOn = 0
WindOn = 0
StarformationCriterion = density
MassiveNuLinRespOn = 0

#  Further parameters of SPH
#  #Only kernel supported by fake_spectra
DensityKernelType = cubic
InitGasTemp = 270.
MinGasTemp = 100

# Memory allocation

PartAllocFactor = 2.0

#----------------------SFR Stuff-------------------------

CritPhysDensity = 0       #  critical physical density for star formation in
#  hydrogen number density in cm^(-3)
CritOverDensity = 1000   #  overdensity threshold value
QuickLymanAlphaProbability = 1 # Set to 1.0 to turn dense gas directly into stars.

SnapshotWithFOF = 1
WindModel = nowind
""".format(n, n)

    with open(directory + 'paramfileVaried{0}.gadget'.format(n), 'w') as f:
        f.write(gadgetParamString)

#make sbatch file, with n/10 sims per sbatch file
#I think we should shoot for 50 sims so that's 5 per sbatch file


dateLine = "date \n"

for j in range(10):
    k = j + 5
    sbatchStringStart = """
#!/bin/bash
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -o lyaV{0}.out -e lyaV{1}.err
#SBATCH --reservation andras_test
#SBATCH -p cca
#SBATCH --qos cca

module load gcc
module load openmpi
module load lib/gsl/2.3


ROOT=/mnt/xfs1/home/landerson/src/MP-Gadget/build

export OMP_NUM_THREADS=1
date
""".format(k, k)

    with open(directory + 'Var{0}.sbatch'.format(k), 'w') as f:
        f.write(sbatchStringStart)
        for i in range(5):
            num = k + 10*i
            genicLine = "srun -n 28 -c 1 --mpi=pmi2 $ROOT/MP-GenIC  paramfileVaried{0}.genic || exit 1 \n".format(num)
            gadgetLine = "srun -n 28 -c 1 --mpi=pmi2 $ROOT/MP-Gadget paramfileVaried{0}.gadget || exit 1 \n".format(num)
            f.write(genicLine)
            f.write(gadgetLine)
            f.write(dateLine)
