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
            print('trying to load 1D ps')
            reload_snapshot=False
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                          reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
        except OSError:

            reload_snapshot = True
            spectra = boxes.SimulationBox(snapshot_num, snapshot_dir, grid_width, spectral_res,
                                          reload_snapshot=reload_snapshot, spectra_savedir=spectra_savedir)
    else:
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

    #python3 PowerSpectrPipeline.py 200 50 0 4 12 14
    #gw = int(sys.argv[1]) #200
    #sr = int(sys.argv[2]) #50
    snap_nums = [4, 7, 9, 12, 14] #[int(s) for s in sys.argv[3:]]
    boxsize = 20. #Mpc/h
    grid_width = 200
    spectral_res = 10*u.km/u.s

    xlim1d = (0.3, 10)
    xlim3d = (0.3, 100)
    ylim_avg = (0.1, 10)

    spectra_dir_pre = '/mnt/ceph/users/landerson/lyalphaVaried'
    snapshot_dir_pre = '/mnt/cephtest/landerson/lyalphaVaried'
    snapshots = np.arange(1, 50.001, 1)
    #snapshot_dir = ['lyalphaVaried1', 'lyalphaVaried2', 'lyalphaVaried3', 'lyalphaVaried4', 'lyalphaFixedA', 'lyalphaFixedB']
    #labels = ['V1', 'V2', 'V3', 'V4', 'FA', 'FB']
    #colors = ['#fdcc8a', '#fc8d59', '#e34a33', '#b30000', '#08519c', '#252525']
    #snap_nums = [0, 2]

    p_corr = boxsize**3.
    k_corr = 2*np.pi/boxsize

    figmf, axmf = plt.subplots()

    #loop over redshift and the grid with associated with it
    #I currently set the grid width to be the same at each redshift, though the resolution is different at different redshifts
    #something to improve in the future
    for sn in snap_nums:
        snapA = '/mnt/cephtest/landerson/lyalphaFixedA/'
        snapB = '/mnt/cephtest/landerson/lyalphaFixedB/'
        p1dA, k1dA, mean_fluxA, redshiftA = get1dps(snapshot_num=sn, snapshot_dir = snapA, reload_snapshot=True, grid_width=grid_width,
                                                    spectral_res=spectral_res, spectra_savedir=snapA)
        p1dB, k1dB, mean_fluxB, redshiftB = get1dps(snapshot_num=sn, snapshot_dir = snapB, reload_snapshot=True, grid_width=grid_width,
                                                    spectral_res=spectral_res, spectra_savedir=snapB)
        p3dA, k3dA = get3dps(snapA, sn)
        p3dB, k3dB = get3dps(snapB, sn)

        spectra1dF = [p1dA, p1dB]
        spectra3dF = [p3dA, p3dB]

        spectra1d= []
        k1d = []
        spectra3d = []
        k3d = []
        legend = []
        #meanflux = []
        #z = []
        fig, ax = plt.subplots(4, figsize=(6, 8)) #len(snap_nums))
        for s in snapshots:
            snap = snapshot_dir_pre + str(s)
            spec = spectra_dir_pre + str(s)
            p1d, k1dnow, mean_flux, redshift = get1dps(snapshot_num=sn, snapshot_dir = snap, reload_snapshot=True, grid_width=grid_width,
                                                       spectral_res=spectral_res, spectra_savedir=spec)
            lnow, = ax[0].loglog(k1dnow, p1d, label=l, color='black', lw=0.5, alpha=0.5)
            ax[0].set_xlim(xlim1d)
            legend.append(lnow)
            p3d, k3dnow = get3dps(snap, sn)
            ax[2].loglog(k3dnow*k_corr, p3d*p_corr, color='black', lw=0.5, alpha=0.5)
            ax[2].set_xlim(xlim3d)
            spectra1d.append(p1d)
            k1d.append(k1dnow)

            spectra3d.append(p3d*p_corr)
            k3d.append(k3dnow*k_corr)

            axmf.scatter(redshift, np.log(mean_flux), edgecolors=c, facecolors='none')

        for kmode, spec, specF, axis, xlim in zip([k1d, k3d], [spectra1d, spectra3d], [spectra1dF, spectra3dF], [ax[1], ax[3]], [xlim1d, xlim3d]):
            spec = np.vstack(spec)
            meanspec = np.mean(spec, axis=0)
            for kk, ss, c in zip(kmode, spec, colors):
                axis.loglog(kk, ss/meanspec, color='black', lw=0.5, alpha=0.5)
                axis.set_xlim(xlim)
            axis.loglog(kk, specF[0]/meanspec, color='red', lw=0.5, alpha=0.75)
            axis.loglog(kk, specF[1]/meanspec, color='blue', lw=0.5, alpha=0.75)
            axis.loglog(kk, 0.5*(specF[0] + specF[1])/meanspec, color='orange', lw=0.5, alpha=0.75)
        #ax[1].plot(k1d[4], 0.5*(spectra1d[4] + spectra1d[5])/np.mean(np.vstack(spectra1d), axis=0), linestyle=':', color='k', label='FixPair Mean')
        #ax[3].plot(k3d[4], 0.5*(spectra3d[4] + spectra3d[5])/np.mean(np.vstack(spectra3d), axis=0), linestyle=':', color='k', label='FixPair Mean')


        for i in [1,3]:
            xlim = ax[i].get_xlim()
            ax[i].plot(xlim, [1.0, 1.0], linestyle='--', alpha=0.5, color='black')
            if i == 3: ax[i].set_ylim(ylim_avg)
            ax[i].set_xlim(xlim)
        ax[3].set_xlabel('k [h/Mpc]')
        ax[0].set_ylabel('1DP')
        ax[1].set_ylabel('1DP/<P>')
        ax[2].set_ylabel('3DP')
        ax[3].set_ylabel('3DP/<P>')
        fig.legend(legend, labels,
                   ncol=len(labels), frameon=False, mode="expand", borderaxespad=0.2, bbox_to_anchor=(0., 0.95, 0.95, 0.))#, loc=3
        ax[0].set_title('z={0:0.2f}  spectral resolution {1}km/s\n'.format(redshift, sr), y=1.25)
        fig.tight_layout()
        fig.savefig('ps_{0:03d}_ngrid{1:03d}_specres{2:03d}.pdf'.format(sn, gw, sr))
    #colors = [l.get_c() for l in legend]

    zz = np.linspace(0, 4, 100)
    axmf.plot(zz, lnMeanFlux(zz), lw=2, color='black')
    axmf.set_ylim(-4, 1)
    axmf.set_xlabel('redshift')
    axmf.set_ylabel('ln mean flux')
    figmf.savefig('meanFlux_ngrid{0:03d}_specres{1:03d}.pdf'.format(gw, sr))
    print('saved meanflux plot for gridwidth {0} and spectral resolution {1}'.format(gw, sr))
