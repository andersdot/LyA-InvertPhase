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


def get1dps(snapshot_dir = '.', snapshot_num=14, grid_width=20, spectral_res=50*u.km/u.s, reload_snapshot=True, label=None, boxsize=20.,
            spectra_savedir=None):

    if reload_snapshot == False:
        try:
            print('trying to load 1D ps from ', spectra_savedir)
            reload_snapshot=False
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                          reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
        except OSError:
            print('tried but failed to load 1d ps doesnt exist, now calculating')
            reload_snapshot = True
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                          reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
        #except KeyError:
        #    print('tried but failed to load 1d ps due to key error, now calculating')
        #    reload_snapshot = True
        #    spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
        #                                  reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)

    else:
        print('calculating 1d power')
        spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                      reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)

    spectra.convert_fourier_units_to_distance = True
    tau = spectra.get_optical_depth()
    tau_scaling = 1.
    mean_flux = np.mean(np.exp(-1.*tau*tau_scaling))
    #print('The mean flux is: ', mean_flux)
    #spectra_box = spectra.skewers_realisation_hydrogen_overdensity()
    spectra_box = spectra.skewers_realisation()
    fourier_estimator_instance = fou.FourierEstimator1D(spectra_box)
    result = fourier_estimator_instance.get_flux_power_1D()
    x = np.arange(len(result))*2*np.pi/boxsize
    return result, x, mean_flux, spectra._redshift


def get3dps(snapshot_directory, snapshot, save_directory):
    filename = save_directory + '/PK-DM-snap_{0:03d}'.format(snapshot)
    if os.path.exists(filename):
        print('loading 3d power from ', filename)
        data = np.genfromtxt(filename, names= ['k', 'p'])
    else:
        print('generating 3d power')
        command = '/mnt/home/landerson/src/GenPK/gen-pk -i {0}/snap_{1:03d} -o {2}'.format(snapshot_directory, snapshot, save_directory)
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

    #python3 PowerSpectrPipeline.py 200 50 0 4 12 14
    #gw = int(sys.argv[1]) #200
    #sr = int(sys.argv[2]) #50
    #snap_nums = [4, 7, 9, 12, 14] #[int(s) for s in sys.argv[3:]]
    snap_nums = [0, 1, 2]
    boxsize = 20. #Mpc/h
    grid_width = 200
    spectral_res = 10*u.km/u.s

    xlim1d = (0.3, 10)
    xlim3d = (0.3, 100)
    ylim_avg = (0.1, 10)


    snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/'
    spectra_savedir_pre = '/home/landerson/lyalpha/'

    spectra_dir_pre = '/home/landerson/lyalpha/'
    snapshot_dir_pre = '/home/fvillaescusa/data/Lya_ncv/'
    simNum = np.arange(0, 49.001, 1).astype(int)

    #snapshot_dir = ['lyalphaVaried1', 'lyalphaVaried2', 'lyalphaVaried3', 'lyalphaVaried4', 'lyalphaFixedA', 'lyalphaFixedB']
    #labels = ['V1', 'V2', 'V3', 'V4', 'FA', 'FB']
    #colors = ['#fdcc8a', '#fc8d59', '#e34a33', '#b30000', '#08519c', '#252525']
    #snap_nums = [0, 2]

    p_corr = boxsize**3.
    k_corr = 2*np.pi/boxsize

    figmf, axmf = plt.subplots()
    savenpz_pre = ['T', 'NCV_0', 'NCV_1']
    snap_pre = ['', 'NCV_0_', 'NCV_1_']
    #loop over redshift and the grid with associated with it
    #I currently set the grid width to be the same at each redshift, though the resolution is different at different redshifts
    #something to improve in the future
    for sn in snap_nums:

        fig, ax = plt.subplots(4, figsize=(6, 8)) #len(snap_nums))

        for sp, snpz in zip(snap_pre, savenpz_pre):
            try:
                data = np.load('spec200_{0}_{1}.npz'.format(snpz, sn) )
                spectra1d = data['p1d']
                spectra3d = data['p3d']
                k1d = data['k1d']
                k3d = data['k3d']
                lnmeanflux = data['lnmeanflux']
                z = data['z']
                dataSaved = True
                print('data loaded from npz for ','spec200_{0}_{1}'.format(snpz, sn))
            except IOError:
                dataSaved = False
                print('loading all data for ', 'spec200_{0}_{1}'.format(snpz, sn))
                spectra1d= []
                k1d = []
                spectra3d = []
                k3d = []
                legend = []
                lnmeanflux = []
                z = []
                
            for s in simNum:
                snap = snapshot_dir_pre + sp + str(s) 
                spec = spectra_dir_pre + sp + str(s) + '/SPECTRA_{0:03d}'.format(sn)
                if not dataSaved:
                    p1d, k1dnow, mean_flux, redshift = get1dps(snapshot_num=sn, snapshot_dir = snap, reload_snapshot=False, grid_width=grid_width,
                                                           spectral_res=spectral_res, spectra_savedir=spec)
                    lnmf = np.log(mean_flux)
                    spectra1d.append(p1d)
                    k1d.append(k1dnow)
                    z.append(redshift)
                    lnmeanflux.append(lnmf)
                if dataSaved:
                    k1dnow = k1d[s-1]
                    p1d = spectra1d[s-1]
                    lnmf = lnmeanflux[s-1]
                    redshift = z[s-1]
                lnow, = ax[0].loglog(k1dnow, p1d, color='black', lw=0.5, alpha=0.1)
                ax[0].set_xlim(xlim1d)
                #legend.append(lnow)
                if not dataSaved:
                    p3d, k3dnow = get3dps(snap, sn, spectra_dir_pre + sp + str(s))
                    spectra3d.append(p3d*p_corr)
                    k3d.append(k3dnow*k_corr)
                if dataSaved:
                    p3d = spectra3d[s-1]/p_corr
                    k3dnow = k3d[s-1]/k_corr
                ax[2].loglog(k3dnow*k_corr, p3d*p_corr, color='black', lw=0.5, alpha=0.1)
                ax[2].set_xlim(xlim3d)

                print(len(p1d), len(p3d))
                axmf.scatter(redshift, lnmf, edgecolors='black', facecolors='none', alpha=0.5)

            if not dataSaved: np.savez('spec200_{0}_{1}'.format(snpz, sn), p1d=spectra1d, p3d=spectra3d, k1d=k1d, k3d=k3d, lnmeanflux=lnmeanflux, z=z)    

"""
            for kmode, spec, axis, xlim in zip([k1d, k3d], [spectra1d, spectra3d], [ax[1], ax[3]], [xlim1d, xlim3d]):
                spec = np.vstack(spec)
                meanspec = np.mean(spec, axis=0)
                for kk, ss in zip(kmode, spec):
                    lt, = axis.loglog(kk, ss/meanspec, color='black', lw=0.5, alpha=0.5)
                    axis.set_xlim(xlim)



        for i in [1,3]:
            xlim = ax[i].get_xlim()
            ax[i].plot(xlim, [1.0, 1.0], linestyle='--', alpha=0.5, color='black')
            if i == 3: ax[i].set_ylim(ylim_avg)
            if i == 1: ax[i].set_ylim(0.7, 1.5)
            ax[i].set_xlim(xlim)
        ax[3].set_xlabel('k [h/Mpc]')
        ax[0].set_ylabel('1DP')
        ax[1].set_ylabel('1DP/<P>')
        ax[2].set_ylabel('3DP')
        ax[3].set_ylabel('3DP/<P>')
        legend = [la, lb, lm, lt]
        labels = ['FA', 'FB', 'FixMean', 'Trad']
        fig.legend(legend, labels,
                   ncol=len(labels), frameon=False, mode="expand", borderaxespad=0.2, bbox_to_anchor=(0., 0.95, 0.95, 0.))#, loc=3
        ax[0].set_title('z={0:0.2f}  spectral resolution {1}km/s\n'.format(redshift, spectral_res), y=1.25)
        fig.tight_layout()
        fig.savefig('ps_{0:03d}_ngrid{1:03d}_specres{2:03d}.pdf'.format(sn, grid_width, int(spectral_res.value)))
    #colors = [l.get_c() for l in legend]

    zz = np.linspace(1, 5, 100)
    axmf.plot(zz, lnMeanFlux(zz), lw=2, color='black')
    axmf.set_ylim(-2.5, 0)
    axmf.set_xlim(4.5, 1.5)
    axmf.set_xlabel('redshift')
    axmf.set_ylabel('ln mean flux')
    figmf.savefig('meanFlux_ngrid{0:03d}_specres{1:03d}.pdf'.format(grid_width, int(spectral_res.value)))
    print('saved meanflux plot for gridwidth {0} and spectral resolution {1}'.format(grid_width, int(spectral_res.value)))
"""
