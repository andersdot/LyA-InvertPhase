import numpy as np
import astropy.units as u
import math as mh
import sys
import power_spectra as spe
import fourier_estimators as fou
import boxes
import bigfile
import subprocess
import os.path


def get1dps(snapshot_dir = '.', snapshot_num=14, grid_width=20, spectral_res=50*u.km/u.s, reload_snapshot=True, label=None, boxsize=20.):

    if reload_snapshot == False:
        try:
            print('trying to load 1D ps')
            reload_snapshot=False
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res, reload_snapshot=reload_snapshot)
        except OSError:

            reload_snapshot = True
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res, reload_snapshot=reload_snapshot)
    else:
        spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res, reload_snapshot=reload_snapshot)



def get3dps(snapshot_directory, snapshot):
    filename = snapshot_directory + '/PK-DM-PART_{0:03d}'.format(snapshot)
    if os.path.exists(filename):
        data = np.genfromtxt(filename, names= ['k', 'p'])
    else:
        command = '/mnt/home/landerson/src/GenPK/gen-pk -i {0}/PART_{1:03d} -o {2}'.format(snapshot_dir_pre + s, sn, snapshot_dir_pre + s)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        exit_code = process.wait()
        data = np.genfromtxt(filename, names=['k', 'p'])
    return data['p'], data['k']

def lnMeanFlux(z):
    return np.log(0.8)*((1. + z)/3.25)**3.2


if __name__ == "__main__":
    import matplotlib as mpl
    mpl.use('pdf')
    import matplotlib.pyplot as plt

    #python3 PowerSpectrPipeline.py 200 50 lyalphaVaried1
    grid_width = int(sys.argv[1]) #200
    spectral_res = int(sys.argv[2])*u.km/u.s
    snapshot_dir = sys.argv[3] #50
    snap_nums = [0, 4, 7, 9, 12, 14] #[int(s) for s in sys.argv[3:]]
    boxsize = 20. #Mpc/h

    xlim1d = (0.3, 10)
    xlim3d = (0.3, 100)
    ylim_avg = (0.1, 10)

    #snapshot_dir_pre = '/mnt/ceph/users/landerson/'
    snapshot_dir_pre = '/mnt/cephtest/landerson/'

    #loop over redshift and the grid with associated with it
    #I currently set the grid width to be the same at each redshift, though the resolution is different at different redshifts
    #something to improve in the future
    for sn in snap_nums:
        get1dps(snapshot_num=sn, snapshot_dir = snapshot_dir_pre + snapshot_dir, reload_snapshot=False, label=l, grid_width=gw, spectral_res=sr)
        p3d, k3dnow = get3dps(snapshot_dir_pre + s, sn)
