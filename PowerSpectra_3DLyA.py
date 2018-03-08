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


def get3dps(snapshot_dir = '.', snapshot_num=14, grid_width=20, spectral_res=50*u.km/u.s, reload_snapshot=True, label=None, boxsize=20.,
            spectra_savedir=None):

    if reload_snapshot == False:
        try:
            print('trying to load 1D ps from ', spectra_savedir)
            reload_snapshot=False
            simulation_box_instance = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                          reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
        except OSError:
            print('tried but failed to load 1d ps doesnt exist, now calculating')
            reload_snapshot = True
            simulation_box_instance = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                          reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
        #except KeyError:
        #    print('tried but failed to load 1d ps due to key error, now calculating')
        #    reload_snapshot = True
        #    spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
        #                                  reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)

    else:
        print('calculating 1d power')
        simulation_box_instance = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                      reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)

    simulation_box_instance.convert_fourier_units_to_distance = True
    delta_flux_box = simulation_box_instance.skewers_realisation()
    k_box = simulation_box_instance.k_box()
    mu_box = simulation_box_instance.mu_box()

    #Binning to match GenPK
    n_k_bins = 6
    n_mu_bins = 4
    k_max = 20. / u.Mpc

    k_min = np.min(k_box[k_box > 0. / u.Mpc])
    k_bin_max = mh.exp(mh.log(k_max.value) + ((mh.log(k_max.value) - mh.log(k_min.value)) / (n_k_bins - 1))) / u.Mpc
    k_bin_edges = np.exp(np.linspace(mh.log(k_min.value), mh.log(k_bin_max.value), n_k_bins + 1)) / u.Mpc
    #k_bin_edges[-2] = k_max #HACK TO FIX BINNING OF NYQUIST FREQUENCY


    mu_bin_edges = np.linspace(0., 1., n_mu_bins + 1)

    print(len(k_box), len(mu_box), len(k_bin_edges), len(mu_bin_edges))
    fourier_estimator_instance = fou.FourierEstimator3D(delta_flux_box)
    #power_binned_k_mu, k_binned_2D = fourier_instance.get_flux_power_3D_two_coords_hist_binned(k_box, np.absolute(mu_box), k_bin_edges, mu_bin_edges, bin_coord2=False, count=False, std_err=False, norm=True)
    #import pdb; pdb.set_trace()
    power_binned, k_binned, mu_binned = fourier_estimator_instance.get_power_3D_two_coords_binned(k_box,np.absolute(mu_box),k_bin_edges,mu_bin_edges,bin_coord2=True)

    return power_binned, k_binned, mu_binned, simulation_box_instance._redshift

if __name__ == "__main__":
    import matplotlib as mpl
    mpl.use('pdf')
    import matplotlib.pyplot as plt

    snap_nums = [1, 2]
    boxsize = 40. #Mpc/h
    grid_width = 400
    spectral_res = 10*u.km/u.s

    xlim1d = (0.3, 10)
    xlim3d = (0.3, 100)
    ylim_avg = (0.1, 10)


    #spectra_dir_pre = '/home/landerson/lyalpha/'
    #snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/'

    snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/40Mpc_512/'
    spectra_dir_pre = '/home/landerson/lyalpha/40Mpc_512/'

    simulation_numbers = np.arange(0, 49.001, 1).astype(int)


    savenpz_pre = ['T', 'NCV_0', 'NCV_1']
    simulation_pre = ['', 'NCV_0_', 'NCV_1_']
    #loop over redshift and the grid with associated with it
    #I currently set the grid width to be the same at each redshift, though the resolution is different at different redshifts
    #something to improve in the future
    for snapshot_number in snap_nums:

        fig, ax = plt.subplots(4, figsize=(6, 8)) #len(snap_nums))

        for sim_pre, snpz in zip(simulation_pre, savenpz_pre):

            specSaveFile = 'spec3d{0}_{1}_{2}.npz'.format(grid_width, snpz, snapshot_number)
            try:
                data = np.load(specSaveFile)
                spectra = data['p3d']
                mu = data['mu']
                k = data['k3d']
                z = data['z']

                dataSaved = True
                print('data loaded from npz for ', specSaveFile)
            except IOError:
                dataSaved = False
                print('loading all data for ', specSaveFile)
                spectra = []
                mu = []
                k = []
                legend = []
                z = []

            for simulation_number in simulation_numbers:
                if (sim_pre == 'NCV_1_') and (simulation_number > 24): continue 
                if (sim_pre == 'NCV_0_') and (simulation_number > 24): continue 

                snapshot_directory = snapshot_dir_pre + sim_pre + str(simulation_number)
                save_directory = spectra_dir_pre + sim_pre + str(simulation_number)
                spectra_directory = save_directory + '/SPECTRA_{0:03d}'.format(snapshot_number)
                if not dataSaved:
                    ps, k_array, mu_array, redshift = get3dps(snapshot_num=snapshot_number, snapshot_dir = snapshot_directory, 
                                                                          reload_snapshot=False, grid_width=grid_width,
                                                                          spectral_res=spectral_res, spectra_savedir=spectra_directory)
                   
                    spectra.append(ps)
                    k.append(k_array)
                    z.append(redshift)
                    mu.append(mu_array)

                if dataSaved:
                    k_array = k[simulation_number]
                    ps = spectra[simulation_number]
                    redshift = z[simulation_number]
                    mu_array = mu[simulation_number]


            if not dataSaved: np.savez(specSaveFile, p3d=spectra, k3d=k, z=z, mu=mu)

