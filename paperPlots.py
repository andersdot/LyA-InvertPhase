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

def plotAllPS(power, k, ax, color='black', alpha=0.1, label=None):
    for i, p, in enumerate(power):
        if i == 0: label = label
        else: label = None
        ax.loglog(k, p, color=color, lw=0.5, alpha=alpha, label=label)


def variance(power, nPaired=25., nTrad=50.):
    mean = np.mean(power, axis=0)
    var = np.sum((power - mean)**2., axis=0)/(len(power) - 1.)
    return var
    #meanT = np.mean(np.vstack(dataT[pkey]), axis=0)
    #meanP = np.mean(np.vstack(paired_p), axis=0)
    #meanT = np.mean(powerT, axis=0)
    #meanP = np.mean(powerP, axis=0)
    #varP = np.sum((powerP - meanP)**2., axis=0)/len(powerP)
    #varT = np.sum((powerT - meanT)**2., axis=0)/len(powerT)



def bias(k, powerT, powerP, nPaired=25., nTrad = 50.):
    meanT = np.mean(powerT, axis=0)
    meanP = np.mean(powerP, axis=0)
    varP = variance(powerP)
    varT = variance(powerT)
    bias = (meanT - meanP)/np.sqrt(varP/nPaired + varT/nTrad)
    return bias

snap_nums = [0, 1, 2]

#snap_nums = [0, 1, 2]
zz = [4, 3, 2]#, 3, 2]
colors = ['C0', 'C1', 'C2']
spectral_res = 10
xlim1d = (0.3, 30)
xlim3d = (0.3, 30)
ylim_avg = (0.3, 2)

savenpz_pre = ['T', 'NCV_0', 'NCV_1']
snap_pre = ['', 'NCV_0_', 'NCV_1_']

bigbox = False
dimension = '3d'


if bigbox:
    grid_width = 400
    boxsize = 40 #Mpc/h
    filenamePostfix = '_1mubin.npz'
else:
    grid_width = 200
    boxsize = 20
    filenamePostfix = '_1mubin_20Mpc.npz'


if dimension == '1d':
    normalization = boxsize
    plotPostFix = ''
else:
    normalization = boxsize**3.
    plotPostFix = '_1mubin'


alpha_line = 0.75
alpha_fill = 0.75

tradColor = 'C0'
pairColor = 'C1'

fig, ax = plt.subplots(3, sharex=True) #, figsize=(6, 6))

for sn, z, c in zip(snap_nums, zz, colors):
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

    dataT = np.load('/Users/landerson/LyA-InvertPhase/goodspec{0}/spec{0}{1}_{2}'.format(dimension, grid_width, sn) + filenamePostfix)
    data0 = np.load('/Users/landerson/LyA-InvertPhase/goodspec{0}/spec{0}{1}_NCV_0_{2}'.format(dimension, grid_width, sn) + filenamePostfix)
    data1 = np.load('/Users/landerson/LyA-InvertPhase/goodspec{0}/spec{0}{1}_NCV_1_{2}'.format(dimension, grid_width, sn) + filenamePostfix)

    shapeT = dataT['power'].shape
    powerT = dataT['power'].reshape(shapeT[0:2])*normalization
    nTrad = shapeT[0]
    if dimension == '1d':
        kT = np.arange(shapeT[1])[np.newaxis, :]*2*np.pi/boxsize
    else:
        kT = dataT['k'].reshape(shapeT[0:2])


    shapeP = data0['power'].shape
    powerP = 0.5*(data0['power'].reshape(shapeP[0:2]) + data1['power'].reshape(shapeP[0:2]))*normalization
    if dimension == '1d':
        kP = np.arange(shapeP[1])[np.newaxis, :]*2*np.pi/boxsize
    else:
        kP = data0['k'].reshape(shapeP[0:2])
    nPaired = shapeP[0]

    print(nTrad, nPaired)

    biasValues = bias(kT, powerT, powerP, nPaired=nPaired, nTrad = nTrad)
    ax[1].semilogx(kP[0], biasValues, color=c, lw=2, label='z='+str(z))

    varT = variance(powerT) #  nPaired=25., nTrad=50.):
    varP = variance(powerP) #,  nPaired=25., nTrad=50.)


    #import pdb; pdb.set_trace()
    uncert = 1/2.*np.sqrt(varT/varP)*np.sqrt(2./nTrad + 2./nPaired)
    yp = varT/varP + uncert
    ym = varT/varP - uncert
    print(yp, ym)
    ax[2].semilogx(kT[0], varT/varP, color=c, lw=2)
    ax[2].fill_between(kT[0], ym,  yp, color=c, alpha=0.75)
    #ax[2].fill_between(kT[0][:,0], ym[:,0],  yp[:,0], color=c, alpha=0.5)

    plotAllPS(powerT, kT[0], ax[0], color=c, alpha=0.1)
    plotAllPS(powerP, kP[0], ax[0], color=c, alpha=0.1)

#ax[0].set_title('spectral resolution {0}km/s\n'.format(spectral_res), y=1.25)

if dimension == '1d':
    ax[0].set_ylabel('$\mathrm{P_F(k) \; [h/Mpc]}$', fontsize=12)
else:
    ax[0].set_ylabel('$\mathrm{3D \; P_F(k) \;[h/Mpc]^3}$', fontsize=12)
ax[1].set_ylabel(r'$\mathrm{(\overline{P_T} - \overline{P_P})/\sigma_{\bar{P_T} - \bar{P_P}}}$', fontsize=11)
ax[2].set_ylabel('$\mathrm{\sigma_T^2/\sigma_P^2}$', fontsize=12)
ax[2].set_xlabel('$\mathrm{k_{\parallel} \; [h/Mpc]}$', fontsize=12)

ax[1].axhline(0, linestyle='--', color='black')
ax[1].axhline(-3, linestyle='--', color='black', lw=0.5)
ax[1].axhline(3, linestyle='--', color='black', lw=0.5)

ax[2].axhline(1, linestyle='--', color='black')

ax[0].set_xlim(0.2, 30)
if dimension == '1d':
    ax[0].set_ylim(1e-6, 1e0)
    ax[2].set_xlabel('$\mathrm{k_{\parallel} \; [h/Mpc]}$', fontsize=12)

else:
    ax[0].set_ylim(1e-6, 1e2)
    ax[2].set_xlabel('$\mathrm{k \;[h/Mpc]}$', fontsize=12)

ax[1].set_xlim(0.2, 30)
ax[1].set_ylim(-5, 5)
ax[2].set_xlim(0.2, 30)
ax[2].set_ylim(0.1, 300)
ax[2].set_yscale('log')

#lgd = pylab.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
lgd = [fig.legend(loc=9, bbox_to_anchor=(0.5, 1.02), ncol=4)]

#labels = ['traditional', 'paired']
#lT3 = mpatches.Patch(color='red', alpha=0.5, label='Traditional')
#lP3 = mpatches.Patch(color='black', alpha=0.5, label='Paired')
#legend3 = [lT3, lP3]

fig.tight_layout()
fig.savefig('varRatio{0}_{1}Mpc{2}.pdf'.format(dimension, boxsize, plotPostFix), additional_artists=lgd,
bbox_inches="tight")
