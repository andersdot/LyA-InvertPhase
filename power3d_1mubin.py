import numpy as np
import astropy.units as u
import math as mh
import sys

import power_spectra as spe
import boxes as box
import fourier_estimators as fou
import utils as uti

def get_k_bin_edges_logspace(n_k_bins, k_box):
    k_max = np.max(k_box) #0.704 / u.Mpc

    k_min = np.min(k_box[k_box > 0. / u.Mpc])
    k_bin_max = mh.exp(mh.log(k_max.value) + ((mh.log(k_max.value) - mh.log(k_min.value)) / (n_k_bins - 1))) / u.Mpc
    return np.exp(np.linspace(mh.log(k_min.value), mh.log(k_bin_max.value), n_k_bins + 1)) / u.Mpc

def get_mu_bin_edges_linspace(n_mu_bins):
    return np.linspace(0., 1., n_mu_bins + 1)


if __name__ == "__main__":
    """Input arguments: Snapshot number; Snapshot directory path; Width of skewer grid in samples;
    Resolution of spectra in km s^{-1}; Spectra directory path (with '/snapdir_XXX' if necessary)"""


import boxes
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
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt

boxsize = 40. #Mpc/h
grid_width = 400
spectral_res = 10*u.km/u.s
n = 50
snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/40Mpc_512/'
spectra_dir_pre = '/home/landerson/lyalpha/40Mpc_512/'
sim_pres = ['', 'NCV_0_', 'NCV_1_']

alphas = [0.4, 0.6, 0.8, 1.0]
simulation_numbers = np.arange(0, n)
colors = mpl.cm.Blues(np.linspace(0.5, 1.0, n))
snapshot_numbers = [0, 1,2]

for sim_pre in sim_pres:

    for snapshot_number in snapshot_numbers:

        snapshot_num = snapshot_number
        p1 = []

        k1 = []

        mu1 = []

        for simulation_number, c in zip(simulation_numbers, colors):
            if(simulation_number > 24) & (sim_pre == 'NCV_0_'):continue
            if(simulation_number > 24) & (sim_pre == 'NCV_1_'):continue

            snapshot_directory = snapshot_dir_pre + sim_pre + str(simulation_number)
            save_directory = spectra_dir_pre + sim_pre + str(simulation_number)
            spectra_directory = save_directory + '/SPECTRA_{0:03d}'.format(snapshot_number)
            simulation_box_instance = boxes.SimulationBox(snapshot_num, snapshot_directory, grid_width, spectral_res,reload_snapshot=False, spectra_savedir=spectra_directory)

            simulation_box_instance.convert_fourier_units_to_distance = True
            delta_flux_box = simulation_box_instance.skewers_realisation()
            k_box = simulation_box_instance.k_box()
            mu_box = simulation_box_instance.mu_box()

            #Binning to match GenPK
            n_k_bins = 15
            n_mu_bins = 1

            k_bin_edges = get_k_bin_edges_logspace(n_k_bins, k_box)
            mu_bin_edges = get_mu_bin_edges_linspace(n_mu_bins)

            fourier_estimator_instance = fou.FourierEstimator3D(delta_flux_box)
            power_binned, k_binned, mu_binned, bin_counts = fourier_estimator_instance.get_power_3D_two_coords_binned(k_box,np.absolute(mu_box),k_bin_edges,mu_bin_edges,count=True)


            np.savez(SPECTRA_SAVEDIR + POWER_SPECTRA_SAVEFILE, power_binned, k_binned, mu_binned, bin_counts)

            p1.append(power_binned)
            k1.append(k_binned)
            m1.append(mu_binned)
            plt.plot(np.log10(k_binned.value), np.log10(power_binned), color=c, label=None)
            import pdb; pdb.set_trace()
        np.savez('spec3d{0}_{1}{2}_1mubin'.format(grid_width, sim_pre, snapshot_number), power=p1, k=k1, mu=mu1)

    plt.legend(loc='best')
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('P(k)')
    plt.tight_layout()
    plt.savefig('test{0}.pdf'.format(sim_pre))
