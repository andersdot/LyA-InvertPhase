import numpy as np
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt

def lnMeanFlux(z):
    return np.log(0.8)*((1. + z)/3.25)**3.2

def rms(x):
    np.sqrt(np.sum(x**2.)/len(x))

snap_nums = [0, 1, 2]
boxsize = 20. #Mpc/h
spectral_res = 10
grid_width = 200
xlim1d = (0.3, 10)
xlim3d = (0.3, 100)
ylim_avg = (0.1, 10)

spectra_savedir_pre = '/home/landerson/lyalpha/'
p_corr = boxsize**3.
k_corr = 2*np.pi/boxsize

figmf, axmf = plt.subplots()
savenpz_pre = ['T', 'NCV_0', 'NCV_1']
snap_pre = ['', 'NCV_0_', 'NCV_1_']

for sn in snap_nums:
    fig, ax = plt.subplots(4, figsize=(6, 8)) #len(snap_nums))
    
    dataT = np.load('spec200_T_{0}.npz'.format(sn))
    data0 = np.load('spec200_NCV_0_{0}.npz'.format(sn))
    data1 = np.load('spec200_NCV_1_{0}.npz'.format(sn))
    
    meanP1T = np.mean(np.vstack(dataT['p1d']), axis=0)
    meanP3T = np.mean(np.vstack(dataT['p3d']), axis=0)

    rmsP1T = np.sqrt(np.sum((dataT['p1d']/meanP1T - 1.)**2., axis=0)/len(dataT['p1d']))
    rmsP1P = np.sqrt(np.sum((0.5*(data0['p1d'] + data1['p1d'])/meanP1T - 1.)**2., axis=0)/len(dataT['p1d']))
    
    ax[1].fill_between(dataT['k1d'][0], 1. + rmsP1T, 1. - rmsP1T, color='black', alpha=0.5)
    ax[1].fill_between(data0['k1d'][0], 1. + rmsP1P, 1. - rmsP1P, color='red', alpha=0.5)

    for p1d, p3d, k1d, k3d in zip(dataT['p1d'], dataT['p3d'], dataT['k1d'], dataT['k3d']):
        lT, = ax[0].loglog(k1d, p1d, color='black', lw=0.5, alpha=0.1)
        ax[1].loglog(k1d, p1d/meanP1T, color='black', lw=0.5, alpha=0.1)
        ax[2].loglog(k3d*k_corr, p3d*p_corr, color='black', lw=0.5, alpha=0.1)
        ax[3].loglog(k3d*k_corr, p3d/meanP3T, color='black', lw=0.5, alpha=0.1)
    for p1d0, p1d1, p3d0, p3d1, k1d0, k1d1, k3d0, k3d1 in zip(data0['p1d'], data1['p1d'], data0['p3d'], data1['p3d'], data0['k1d'], data1['k1d'], data0['k3d'], data1['k3d']):
        lP, = ax[0].loglog(k1d0, 0.5*(p1d0 + p1d1), color='red', lw=0.5, alpha=0.1)
        ax[1].loglog(k1d0, 0.5*(p1d0 +  p1d1)/meanP1T, color='red', lw=0.1, alpha=0.25)
        ax[2].loglog(k3d0*k_corr, 0.5*(p3d0 + p3d1)*p_corr, color='red', lw=0.5, alpha=0.1)
        ax[3].loglog(k3d0*k_corr, 0.5*(p3d0 + p3d1)/meanP3T, color='red', lw=0.5, alpha=0.25)

    ax[0].set_xlim(xlim1d)
    ax[2].set_xlim(xlim3d)
    axmf.scatter(dataT['z'], dataT['lnmeanflux'], edgecolors='black', facecolors='none', alpha=0.5)
    axmf.scatter(data0['z'], data0['lnmeanflux'], edgecolors='red', facecolors='none', alpha=0.5)
    axmf.scatter(data1['z'], data1['lnmeanflux'], edgecolors='blue', facecolors='none', alpha=0.5)
    
    
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
    legend = [lT, lP]
    labels = ['traditional', 'paired']
    fig.legend(legend, labels,
               ncol=len(labels), frameon=False, mode="expand", borderaxespad=0.2, bbox_to_anchor=(0., 0.95, 0.95, 0.))#, loc=3
    ax[0].set_title('z={0:0.2f}  spectral resolution {1}km/s\n'.format(dataT['z'][0], spectral_res), y=1.25)
    fig.tight_layout()
    fig.savefig('ps_{0:03d}_ngrid{1:03d}_specres{2:03d}.pdf'.format(sn, grid_width, int(spectral_res)))
    #colors = [l.get_c() for l in legend]

zz = np.linspace(1, 5, 100)
axmf.plot(zz, lnMeanFlux(zz), lw=2, color='black')
axmf.set_ylim(-2.5, 0)
axmf.set_xlim(4.5, 1.5)
axmf.set_xlabel('redshift')
axmf.set_ylabel('ln mean flux')
figmf.savefig('meanFlux_ngrid{0:03d}_specres{1:03d}.pdf'.format(grid_width, int(spectral_res)))
print('saved meanflux plot for gridwidth {0} and spectral resolution {1}'.format(grid_width, int(spectral_res)))
