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


def fluxSkewers(snapshot_dir = '.', snapshot_num=14, grid_width=20, spectral_res=50*u.km/u.s, reload_snapshot=True, spectra_savedir=None):
    """ generate spectra skewers, either from file or on the fly and save if they haven't been generated before""" 
    if reload_snapshot == False:
        try:
            print('trying to load 1D ps')
            reload_snapshot=False
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res, reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
        except OSError:
            print('tried to load spectra but doesnt exist, so generating new spectra')
            reload_snapshot = True
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res, reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
            spectra.spectra_instance.save_file()
    else:
        print('generating new spectra')
        spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res, reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
        spectra.spectra_instance.save_file()

    return spectra

def fluxps1d(snapshot_dir = '.', snapshot_num=14, grid_width=20, spectral_res=50*u.km/u.s, reload_snapshot=True, label=None, boxsize=20., spectra_savedir=None):
    spectra = fluxSkewers(snapshot_dir = snapshot_dir, snapshot_num=snapshot_num, 
                          grid_width=grid_width, spectral_res=spectral_res, 
                          reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
    spectra.convert_fourier_units_to_distance = True
    tau = spectra.get_optical_depth()
    tau_scaling = 1.
    mean_flux = np.mean(np.exp(-1.*tau*tau_scaling))
    print('The mean flux is: ', mean_flux)                                                                                                                                        
    spectra_box = spectra.skewers_realisation()
    fourier_estimator_instance = fou.FourierEstimator1D(spectra_box)
    result = fourier_estimator_instance.get_flux_power_1D()
    x = np.arange(len(result))*2*np.pi/boxsize
    return result, x, mean_flux, spectra._redshift
    
def fluxps3d(snapshot_num=14, snapshot_dir='.', grid_width_in_samps=20, spectrum_resolution=50, reload_snapshot=False, spectra_savefile_root='gridded_spectra', spectra_savedir=None, power_spectra_savefile='/power_spectra.npz'):

    simulation_box_instance = box.SimulationBox(snapshot_num, snapshot_dir, grid_width_in_samps, spectrum_resolution, reload_snapshot=reload_snapshot, spectra_savefile_root=spectra_savefile_root, spectra_savedir=specra_savedir)

    simulation_box_instance.convert_fourier_units_to_distance = True
    delta_flux_box = simulation_box_instance.skewers_realisation()
    k_box = simulation_box_instance.k_box()
    mu_box = simulation_box_instance.mu_box()

    #Binning to match GenPK
    n_k_bins = 15
    n_mu_bins = 1
    k_max = np.max(k_box) #0.704 / u.Mpc

    k_min = np.min(k_box[k_box > 0. / u.Mpc])
    k_bin_max = mh.exp(mh.log(k_max.value) + ((mh.log(k_max.value) - mh.log(k_min.value)) / (n_k_bins - 1))) / u.Mpc
    k_bin_edges = np.exp(np.linspace(mh.log(k_min.value), mh.log(k_bin_max.value), n_k_bins + 1)) / u.Mpc

    mu_bin_edges = np.linspace(0., 1., n_mu_bins + 1)

    fourier_estimator_instance = fou.FourierEstimator3D(delta_flux_box)
    power_binned, k_binned, mu_binned, bin_counts = fourier_estimator_instance.get_power_3D_two_coords_binned(k_box,np.absolute(mu_box),k_bin_edges,mu_bin_edges,count=True)
    np.savez(spectra_savedir + power_spectra_savefile, power_binned, k_binned, mu_binned, bin_counts)
    return power_binned, k_binned, mu_binned, bin_counts

def matterps3d(snapshot_directory, snapshot, snapshot_save_directory):
    filename = snapshot_save_directory + '/PK-DM-snap_{0:03d}'.format(snapshot)
    try:
        data = np.genfromtxt(filename, names= ['k', 'p'])
    except OSError:
        command = '/home/landerson/src/GenPK/gen-pk -i {0}/PART_{1:03d} -o {2}'.format(snapshot_directory, snapshot, snapshot_save_directory)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        exit_code = process.wait()
        try: data = np.genfromtxt(filename, names=['k', 'p'])
        except OSError: return 0, 0
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
    #snap_nums = [4, 7, 9, 12, 14] #[int(s) for s in sys.argv[3:]]
    snap_nums = [0, 1, 2]
    boxsize = 20. #Mpc/h

    xlim1d = (0.3, 10)
    xlim3d = (0.3, 100)
    ylim_avg = (0.1, 10)

    p_corr = boxsize**3.
    k_corr = 2*np.pi/boxsize

    fig, ax = plt.subplots(3, figsize=(5, 15))

    #snapshot_dir_pre = '/mnt/ceph/users/landerson/'
    snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/40Mpc_512/'
    #snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/'
    spectra_savedir_pre = '/home/landerson/lyalpha/40Mpc_512/'
    #loop over redshift and the grid with associated with it
    #I currently set the grid width to be the same at each redshift, though the resolution is different at different redshifts
    #something to improve in the future
    for sn in snap_nums:
        print('now doing ', sn)
        
        p1d, k1d, mean_flux, redshift = get1dps(snapshot_num=sn, snapshot_dir = snapshot_dir_pre + snapshot_dir, reload_snapshot=False, 
                                                grid_width=grid_width, spectral_res=spectral_res, boxsize=boxsize, spectra_savedir=spectra_savedir_pre+snapshot_dir+'/SPECTRA_{0:03d}/'.format(sn))
        p3d, k3d = get3dps(snapshot_dir_pre + snapshot_dir, sn, spectra_savedir_pre+snapshot_dir)
        ax[0].loglog(k1d, p1d, label=redshift)
        ax[1].loglog(k3d*k_corr, p3d*p_corr, label=redshift)
        ax[2].plot(redshift, np.log(mean_flux))
    zz = np.linspace(0, 4, 100)
    ax[2].plot(zz, lnMeanFlux(zz), lw=2, color='black')
    ax[0].set_ylabel('LyA 1D PS')
    ax[1].set_ylabel('Matter 3D PS')
    ax[1].set_ylabel('k [h/Mpc]')
    ax[2].set_ylabel('ln Mean Flux')
    ax[2].set_xlabel('redshift')
    plt.tight_layout()
    plt.legend()
    fig.savefig('PS' + snapshot_dir + '.pdf')
