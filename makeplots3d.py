import numpy as np
import matplotlib as mpl
#mpl.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import ascii
from scipy.interpolate import interp1d

def lnMeanFlux(z):
    return np.log(0.8)*((1. + z)/3.25)**3.2

def rms(x):
    np.sqrt(np.sum(x**2.)/len(x))

snap_nums = [2]

#snap_nums = [0, 1, 2]
zz = [4]#, 3, 2]
boxsize = 20. #Mpc/h
spectral_res = 10
grid_width = 200
xlim1d = (0.3, 30)
xlim3d = (0.3, 30)
ylim_avg = (0.3, 2)

savenpz_pre = ['T', 'NCV_0', 'NCV_1']
snap_pre = ['', 'NCV_0_', 'NCV_1_']

grid_width = 400
boxsize = 40
fixK = 2.

alpha_line = 0.1
alpha_fill = 0.75

for sn, z in zip(snap_nums, zz):
    fig3d, ax3d = plt.subplots(3, figsize=(6, 6))
    figvar, axvar = plt.subplots(2)
    dataT = np.load('/Users/landerson/LyA-InvertPhase/spec3d{0}_T_{1}_small.npz'.format(grid_width, sn))
    data0 = np.load('/Users/landerson/LyA-InvertPhase/spec3d{0}_NCV_0_{1}_small.npz'.format(grid_width, sn))
    data1 = np.load('/Users/landerson/LyA-InvertPhase/spec3d{0}_NCV_1_{1}_small.npz'.format(grid_width, sn))

    #PkFranciscoFile = 'Pk_m_mean_z={0}.txt'.format(z)
    #PkFrancisco = ascii.read(PkFranciscoFile, names=['k', 'mean', 'var'])

    #PkFranciscoFile = 'Pk_m_mean_NCV_z={0}.txt'.format(z)
    #PkFranciscoNCV = ascii.read(PkFranciscoFile, names=['k', 'mean', 'var'])

    for pkey, kkey, mukey, c in zip(['p1', 'p2', 'p3', 'p4'], ['k1', 'k2', 'k3', 'k4'], ['mu1', 'mu2', 'mu3', 'mu4'], ['C0', 'C1', 'C2', 'C3']):
        minmu = np.min(data0[mukey])
        maxmu = np.max(data0[mukey])
        label = '${0:0.2f}<\mu<{1:0.2f}$'.format(minmu, maxmu)
        paired_p = 0.5*(data0[pkey] + data1[pkey])
        meanT = np.mean(np.vstack(dataT[pkey]), axis=0)
        meanP = np.mean(np.vstack(paired_p), axis=0)
        stdP = np.sqrt(np.sum((paired_p - meanP)**2., axis=0)/len(data0[pkey]))
        stdT = np.sqrt(np.sum((dataT[pkey] - meanT)**2., axis=0)/len(dataT[pkey]))

        axvar[0].semilogx(dataT[kkey][0], stdT/(stdP*np.sqrt(2)), label=label, color=c, lw=2)

        lT3 = ax3d[2].fill_between(dataT[kkey][0], (meanT + stdT)/meanT, (meanT - stdT)/meanT, color='black', alpha=alpha_fill-0.2)
        lP3 = ax3d[2].fill_between(data0[kkey][0], (meanP + stdP)/meanT, (meanP - stdP)/meanT, color='red', alpha=alpha_fill-0.3)
        for p1d, k1d in zip(paired_p, data0[kkey]):
            ax3d[0].loglog(k1d, p1d, color='red', lw=0.5, alpha=alpha_line)
            ax3d[1].semilogx(k1d, (p1d - meanT)/stdT, color='red', lw=0.5, alpha=alpha_line+0.2)


        for p1d, k1d in zip(dataT[pkey], dataT[kkey]):
            ax3d[0].loglog(k1d, p1d, color='black', lw=0.5, alpha=alpha_line)
            ax3d[1].semilogx(k1d, (p1d - meanT)/stdT, color='black', lw=0.5, alpha=alpha_line+0.2)
        axvar[1].semilogx(data0[kkey][0], (meanT - meanP)/np.sqrt(stdP**2. + stdT**2.), color=c, lw=2)
    ax3d[1].set_ylim(-2, 2)
    ax3d[2].set_ylim(0.5, 1.5)
    ax3d[2].set_xscale('log')

    ax3d[0].set_xlim(0.1, 40)
    ax3d[1].set_xlim(0.1, 40)
    ax3d[2].set_xlim(0.1, 40)

    for ax in [ax3d]:
        ax[2].axhline(1.0, linestyle='--', alpha=0.5, color='black')
        ax[1].axhline(0.0, linestyle='--', alpha=0.5, color='black')
        ax[0].set_title('z={0:0.2f}  spectral resolution {1}km/s\n'.format(z, spectral_res), y=1.25)


    ax3d[2].set_xlabel('k [h/Mpc]')
    ax3d[0].set_ylabel('3D Matter P')
    ax3d[1].set_ylabel('$\Delta/\sigma_T$')
    ax3d[2].set_ylabel('($<P>\pm\sigma)/<P_T>$')

    labels = ['traditional', 'paired']

    lT3 = mpatches.Patch(color='red', alpha=0.5, label='Traditional')
    lP3 = mpatches.Patch(color='black', alpha=0.5, label='Paired')
    legend3 = [lT3, lP3]

    fig3d.tight_layout()
    fig3d.savefig('ps3d_{0:03d}_{1}Mpc.pdf'.format(sn, boxsize))
    #colors = [l.get_c() for l in legend]

    axvar[0].legend(loc='upper right', ncol=2)
    axvar[0].set_xlim(0.1, 40)
    axvar[1].set_xlim(0.1, 40)
    axvar[0].set_ylim(0.9, 6)
    axvar[1].set_ylim(-1.0, 1.0)
    axvar[0].set_ylabel('$\sigma_T/\sigma_P$')
    axvar[0].axhline(1, linestyle='--', color='black')
    axvar[1].set_xlabel('k [h/Mpc]')
    axvar[1].set_ylabel(r'$(\overline{P_T} - \overline{P_P})/\sqrt{(\sigma_T^2 + \sigma_P^2)}$')
    axvar[0].set_title('z={0:0.2f}  1D spectral resolution {1}km/s\n'.format(z, spectral_res))
    figvar.savefig('varRatio_{0:03d}_{1}Mpc.pdf'.format(sn, boxsize))
