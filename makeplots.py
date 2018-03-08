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

snap_nums = [0, 1, 2]
zz = [4, 3, 2]
boxsize = 20. #Mpc/h
spectral_res = 10
grid_width = 200
xlim1d = (0.3, 30)
xlim3d = (0.3, 30)
ylim_avg = (0.3, 2)

figmf, axmf = plt.subplots()
savenpz_pre = ['T', 'NCV_0', 'NCV_1']
snap_pre = ['', 'NCV_0_', 'NCV_1_']

grid_width = 400
boxsize = 40
fixK = 2.

alpha_line = 0.1
alpha_fill = 0.75

for sn, z in zip(snap_nums, zz):
    fig1d, ax1d = plt.subplots(3, figsize=(6, 6)) #len(snap_nums))
    fig3d, ax3d = plt.subplots(3, figsize=(6, 6))
    figvar, axvar = plt.subplots()
    dataT = np.load('/Users/landerson/lyalpha/spec{0}_T_{1}.npz'.format(grid_width, sn))
    data0 = np.load('/Users/landerson/lyalpha/spec{0}_NCV_0_{1}.npz'.format(grid_width, sn))
    data1 = np.load('/Users/landerson/lyalpha/spec{0}_NCV_1_{1}.npz'.format(grid_width, sn))

    #PkFranciscoFile = 'Pk_m_mean_z={0}.txt'.format(z)
    #PkFrancisco = ascii.read(PkFranciscoFile, names=['k', 'mean', 'var'])

    #PkFranciscoFile = 'Pk_m_mean_NCV_z={0}.txt'.format(z)
    #PkFranciscoNCV = ascii.read(PkFranciscoFile, names=['k', 'mean', 'var'])

    paired_p1d = 0.5*(data0['p1d'] + data1['p1d'])
    #paired_p3d = 0.5*(data0['p3d'] + data1['p3d'])

    meanP1T = np.mean(np.vstack(dataT['p1d']), axis=0)
    #meanP3T = np.mean(np.vstack(dataT['p3d']), axis=0)

    meanP1P = np.mean(np.vstack(paired_p1d), axis=0)
    #meanP3P = np.mean(np.vstack(paired_p3d), axis=0)

    stdP1P = np.sqrt(np.sum((paired_p1d - meanP1P)**2., axis=0)/len(data0['p1d']))
    #stdP3P = np.sqrt(np.sum((paired_p3d - meanP3P)**2., axis=0)/len(data0['p3d']))

    stdP1T = np.sqrt(np.sum((dataT['p1d'] - meanP1T)**2., axis=0)/len(dataT['p1d']))
    #stdP3T = np.sqrt(np.sum((dataT['p3d'] - meanP3T)**2., axis=0)/len(dataT['p3d']))
    #import pdb; pdb.set_trace()
    axvar.semilogx(dataT['k1d'][0]/fixK, stdP1T/(stdP1P*np.sqrt(2)), label='1D')
    #axvar.semilogx(PkFrancisco['k'], PkFrancisco['var']/(PkFranciscoNCV['var']*np.sqrt(2)), label='3D')
    axvar.legend(loc='upper left')
    axvar.set_xlim(0.1, 40)
    axvar.set_ylim(0.9, 6)
    axvar.set_ylabel('$\sigma_T/\sigma_P$')
    axvar.axhline(1, linestyle='--', color='black')
    axvar.set_xlabel('k [h/Mpc]')
    axvar.set_title('z={0:0.2f}  1D spectral resolution {1}km/s\n'.format(dataT['z'][0], spectral_res))
    figvar.savefig('varRatio_{0:03d}_{1}Mpc.pdf'.format(sn, boxsize))
    #rmsP1T = np.sqrt(np.sum((dataT['p1d']/meanP1T - 1.)**2., axis=0)/len(dataT['p1d']))
    #rmsP1P = np.sqrt(np.sum((0.5*(data0['p1d'][0:25] + data1['p1d'][0:25])/meanP1T - 1.)**2., axis=0)/len(data0['p1d'][0:25]))

    #rmsP3T = np.sqrt(np.sum((dataT['p3d']*p_corr/meanP3T - 1.)**2., axis=0)/len(dataT['p3d']))
    #rmsP3P = np.sqrt(np.sum((0.5*(data0['p3d'][0:25]*p_corr + data1['p3d'][0:25]*p_corr)/meanP3T - 1.)**2., axis=0)/len(data0['p3d'][0:25]))

    lT1 = ax1d[2].fill_between(dataT['k1d'][0]/fixK, (meanP1T + stdP1T)/meanP1T, (meanP1T - stdP1T)/meanP1T, color='black', alpha=alpha_fill-0.2)
    lP1 = ax1d[2].fill_between(data0['k1d'][0]/fixK, (meanP1P + stdP1P)/meanP1T, (meanP1P - stdP1P)/meanP1T, color='red', alpha=alpha_fill-0.3)
    #lT3 = ax3d[2].fill_between(dataT['k3d'][0], (meanP3T + stdP3T)/meanP3T, (meanP3T - stdP3T)/meanP3T, color='black', alpha=alpha_fill-0.2)
    #lP3 = ax3d[2].fill_between(data0['k3d'][0], (meanP3P + stdP3P*np.sqrt(2))/meanP3T, (meanP3P - stdP3P*np.sqrt(2))/meanP3T, color='red', alpha=alpha_fill-0.3)

    #print(np.min(dataT['k3d'][0]), np.max(dataT['k3d'][0]))
    #print(np.min(PkFrancisco['k']), np.max(PkFrancisco['k']))

    #interpP3D = interp1d(dataT['k3d'][0], meanP3T, bounds_error=False, fill_value='extrapolate', kind='quadratic')


    #ax3d[2].plot(PkFrancisco['k'], (PkFrancisco['mean'] + PkFrancisco['var'])/PkFrancisco['mean'], color='blue')
    #ax3d[2].plot(PkFrancisco['k'], (PkFrancisco['mean'] - PkFrancisco['var'])/PkFrancisco['mean'], color='blue')

    #ax3d[2].plot(PkFranciscoNCV['k'], (PkFranciscoNCV['mean'] + PkFranciscoNCV['var']*np.sqrt(2))/PkFrancisco['mean'], color='green')
    #ax3d[2].plot(PkFranciscoNCV['k'], (PkFranciscoNCV['mean'] - PkFranciscoNCV['var']*np.sqrt(2))/PkFrancisco['mean'], color='green')
    #for p3d, k3d in zip(dataT['p3d'],dataT['k3d']):
    #    ax3d[0].loglog(k3d, p3d, color='black', lw=0.5, alpha=alpha_line)
    #    ax3d[1].semilogx(k3d, (p3d - meanP3T)/stdP3T, color='black', lw=0.5, alpha=alpha_line+0.2)
    for p1d, k1d in zip(paired_p1d, data0['k1d']/fixK):
        ax1d[0].loglog(k1d, p1d, color='red', lw=0.5, alpha=alpha_line)
        ax1d[1].semilogx(k1d, (p1d - meanP1T)/stdP1T, color='red', lw=0.5, alpha=alpha_line+0.2)

    #for p3d, k3d in zip(dataT['p3d'], dataT['k3d']):
    #    ax3d[0].loglog(k3d, p3d, color='black', lw=0.5, alpha=alpha_line)
    #    ax3d[1].semilogx(k3d, (p3d - meanP3T)/stdP3T, color='black', lw=0.5, alpha=alpha_line+0.2)
    for p1d, k1d in zip(dataT['p1d'], dataT['k1d']/fixK):
        ax1d[0].loglog(k1d, p1d, color='black', lw=0.5, alpha=alpha_line)
        ax1d[1].semilogx(k1d, (p1d - meanP1T)/stdP1T, color='black', lw=0.5, alpha=alpha_line+0.2)

    #for p3d, k3d in zip(paired_p3d, data0['k3d']):
    #    ax3d[0].loglog(k3d, p3d, color='red', lw=0.5, alpha=alpha_line)
    #    ax3d[1].semilogx(k3d, (p3d - meanP3T)/stdP3T, color='red', lw=0.5, alpha=alpha_line+0.2)

    for p1d, k1d in zip(paired_p1d, data0['k1d']/fixK):
        ax1d[0].loglog(k1d, p1d, color='red', lw=0.5, alpha=alpha_line)
        ax1d[1].semilogx(k1d, (p1d - meanP1T)/stdP1T, color='red', lw=0.5, alpha=alpha_line+0.2)

    ax1d[1].set_ylim(-2, 2)
    ax3d[1].set_ylim(-2, 2)
    ax1d[2].set_ylim(0.8, 1.2)
    ax3d[2].set_ylim(0.5, 1.5)
    ax1d[2].set_xscale('log')
    ax3d[2].set_xscale('log')
    axmf.scatter(dataT['z'], dataT['lnmeanflux'], edgecolors='black', facecolors='none', alpha=0.5)
    axmf.scatter(data0['z'], data0['lnmeanflux'], edgecolors='red', facecolors='none', alpha=0.5)
    axmf.scatter(data1['z'], data1['lnmeanflux'], edgecolors='blue', facecolors='none', alpha=0.5)
    ax1d[0].set_xlim(0.1, 30)
    ax1d[1].set_xlim(0.1, 30)
    ax1d[2].set_xlim(0.1, 30)
    ax3d[0].set_xlim(0.1, 150)
    ax3d[1].set_xlim(0.1, 150)
    ax3d[2].set_xlim(0.1, 150)

    for ax in [ax1d, ax3d]:
        ax[2].axhline(1.0, linestyle='--', alpha=0.5, color='black')
        ax[1].axhline(0.0, linestyle='--', alpha=0.5, color='black')
        ax[0].set_title('z={0:0.2f}  spectral resolution {1}km/s\n'.format(dataT['z'][0], spectral_res), y=1.25)

        #if i == 3: ax[i].set_ylim(ylim_avg)
        #if i == 1: ax[i].set_ylim(0.8, 1.2)
        #ax[i].set_xlim(xlim)
    ax1d[2].set_xlabel('k [h/Mpc]')
    ax3d[2].set_xlabel('k [h/Mpc]')
    ax1d[0].set_ylabel('LyA 1D P')
    ax1d[1].set_ylabel('$\Delta/\sigma_T$')
    ax1d[2].set_ylabel('($<P>\pm\sigma)/<P_T>$')
    ax3d[0].set_ylabel('3D Matter P')
    ax3d[1].set_ylabel('$\Delta/\sigma_T$')
    ax3d[2].set_ylabel('($<P>\pm\sigma)/<P_T>$')
    #ax[0].set_xlim(xlim)
    #ax[2].set_xlim(xlim)
    #ax[0].set_xlim(xlim1d)
    #ax[2].set_xlim(xlim3d)
    #ax[1].set_xlim(xlim1d)
    #ax[3].set_xlim(xlim3d)
    labels = ['traditional', 'paired']

    lT1 = mpatches.Patch(color='red', alpha=0.5, label='Traditional')
    lT3 = lT1
    lP1 = mpatches.Patch(color='black', alpha=0.5, label='Paired')
    lP3 = lP1
    legend1 = [lT1, lP1]
    legend3 = [lT3, lP3]

    for fig, ax, legend in zip([fig1d, fig3d], [ax1d[2], ax3d[2]], [legend1, legend3]):
        #ax.legend() #handles=legend)#,
               #ncol=2, frameon=False, mode="expand", borderaxespad=0.2, bbox_to_anchor=(0., 0.95, 0.95, 0.))#, loc=3
        fig.tight_layout()
    fig1d.savefig('ps1d_{0:03d}_{1}Mpc.pdf'.format(sn, boxsize))
    fig3d.savefig('ps3d_{0:03d}_{1}Mpc.pdf'.format(sn, boxsize))
    #colors = [l.get_c() for l in legend]

zz = np.linspace(1, 5, 100)
axmf.plot(zz, lnMeanFlux(zz), lw=2, color='black')
axmf.set_ylim(-1.2, 0)
axmf.set_xlim(4.5, 1.5)
axmf.set_xlabel('redshift')
axmf.set_ylabel('ln mean flux')
figmf.savefig('meanFlux_ngrid{0:03d}_specres{1:03d}_{2}Mpc.pdf'.format(grid_width, int(spectral_res), boxsize))
print('saved meanflux plot for gridwidth {0} and spectral resolution {1}'.format(grid_width, int(spectral_res)))
