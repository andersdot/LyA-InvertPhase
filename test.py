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

snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/40Mpc_512/'
spectra_dir_pre = '/home/landerson/lyalpha/40Mpc_512/'
sim_pres = ['', 'NCV_0_', 'NCV_1_']
n = 10
alphas = [0.4, 0.6, 0.8, 1.0]
simulation_numbers = np.arange(0, n)
colors = mpl.cm.Blues(np.linspace(0.5, 1.0, n))
snapshot_numbers = [1,2]

for sim_pre in sim_pres:

    for snapshot_number in [snapshot_numbers[-1]]:
        
        snapshot_num = snapshot_number
        p1 = []
        p2 = [] 
        p3 = []
        p4 = []
        
        k1 = []
        k2 = []
        k3 = []
        k4 = []
        
        mu1 = []
        mu2 = []
        mu3 = []
        mu4 = []
        

        for simulation_number, c in zip(simulation_numbers, colors):
            snapshot_directory = snapshot_dir_pre + sim_pre + str(simulation_number)
            save_directory = spectra_dir_pre + sim_pre + str(simulation_number)
            spectra_directory = save_directory + '/SPECTRA_{0:03d}'.format(snapshot_number)
            
            
            simulation_box_instance = boxes.SimulationBox(snapshot_num, snapshot_directory, grid_width, spectral_res,reload_snapshot=False, spectra_savedir=spectra_directory)
            simulation_box_instance.convert_fourier_units_to_distance = True
            delta_flux_box = simulation_box_instance.skewers_realisation()
            k_box = simulation_box_instance.k_box()
            mu_box = simulation_box_instance.mu_box()
            
            n_k_bins = 6
            n_mu_bins = 4
            k_max = 20. / u.Mpc
            k_min = np.min(k_box[k_box > 0. / u.Mpc])
            k_bin_max = mh.exp(mh.log(k_max.value) + ((mh.log(k_max.value) - mh.log(k_min.value)) / (n_k_bins - 1))) / u.Mpc
            k_bin_edges = np.exp(np.linspace(mh.log(k_min.value), mh.log(k_bin_max.value), n_k_bins + 1)) / u.Mpc
            #k_bin_edges[-2] = k_max #HACK TO FIX BINNING OF NYQUIST FREQUENCY                                                                             
            mu_bin_edges = np.linspace(0., 1., n_mu_bins + 1)
            fourier_estimator_instance = fou.FourierEstimator3D(delta_flux_box)
            result = fourier_estimator_instance.get_power_3D_two_coords_binned(k_box,np.absolute(mu_box),k_bin_edges,mu_bin_edges,bin_coord2=True)
            
            power = result[0]
            k = result[1]
            mu = result[2]
            
            for i, (plist, klist, mulist) in enumerate(zip([p1, p2, p3, p4], [k1, k2, k3, k4], [mu1, mu2, mu3, mu4])):
                plist.append(power[:,i])
                klist.append(k[:,i])
                mulist.append(mu[:,i])
                
            for i in range(n_mu_bins):
                if simulation_number == 0:
                    minmu = np.min(mu[:,i])
                    maxmu = np.max(mu[:,i])
                    label = 'mu{0:0.2f}-{1:0.2f}'.format(minmu, maxmu)
                else: label = None
                plt.plot(np.log10(k[:,i]), np.log10(power[:,i]), color=c, alpha=alphas[i], label=label)

    np.savez('spec3d{0}_{1}_{2}'.format(grid_width, sim_pre, snapshot_number), p1=p1, p2=p2, p3=p3, p4=p4, k1=k1, k2=k2, k3=k3, k4=k4, mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4)

    plt.legend(loc='best')
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('P(k)')
    plt.tight_layout()
    plt.savefig('test{0}.pdf'.format(sim_pre))
