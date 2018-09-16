from pylab import *
import numpy as np
from matplotlib.ticker import ScalarFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LogNorm
from matplotlib.patches import Ellipse
rcParams["mathtext.fontset"]='cm'
rcParams['axes.linewidth'] = 2 #set the value globally


import matplotlib as mpl
#mpl.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import ascii
from scipy.interpolate import interp1d

def create_mean_std(bigbox, dimension, sn, fixed=False):
    if bigbox:
        grid_width = 400
        boxsize = 40 #Mpc/h
        filenamePostfix = '_1mubin.npz'
    else:
        grid_width = 200
        boxsize = 20
        filenamePostfix = '_1mubin_20Mpc.npz'

    dataT = np.load('/Users/landerson/LyA-InvertPhase/goodspec{0}/spec{0}{1}_{2}'.format(dimension, grid_width, sn) + filenamePostfix)
    data0 = np.load('/Users/landerson/LyA-InvertPhase/goodspec{0}/spec{0}{1}_NCV_0_{2}'.format(dimension, grid_width, sn) + filenamePostfix)
    data1 = np.load('/Users/landerson/LyA-InvertPhase/goodspec{0}/spec{0}{1}_NCV_1_{2}'.format(dimension, grid_width, sn) + filenamePostfix)


    if dimension == '1d':
        normalization = boxsize
    else:
        normalization = boxsize**3.

    shapeT = dataT['power'].shape
    powerT = dataT['power'].reshape(shapeT[0:2])*normalization

    nTrad = shapeT[0]
    if dimension == '1d':
        kT = np.arange(shapeT[1])[np.newaxis, :]*2*np.pi/boxsize
    else:
        kT = dataT['k'].reshape(shapeT[0:2])
    shapeP = data0['power'].shape
    powerP = 0.5*(data0['power'].reshape(shapeP[0:2]) + data1['power'].reshape(shapeP[0:2]))*normalization
    if fixed: powerP = data1['power'].reshape(shapeP[0:2])*normalization
    if dimension == '1d':
        #import pdb; pdb.set_trace()
        kP = np.arange(shapeP[1])[np.newaxis, :]*2*np.pi/boxsize
    else:
        kP = data0['k'].reshape(shapeP[0:2])
    #nPaired = shapeP[0]
    #import pdb; pdb.set_trace()

    meanT, varT = variance(powerT) #  nPaired=25., nTrad=50.):
    meanP, varP = variance(powerP) #,  nPaired=25., nTrad=50.)
    import pdb; pdb.set_trace()
    return kT[0], meanT, np.std(powerT, axis=0), kP[0], meanP, np.std(powerP, axis=0)


def make_plot(power, zz, colors, snaps, fout,f_1,f_2,fr, bigbox, dimension,
              nstandard,  npaired, x_min, x_max, y_min_1, y_max_1, y_lim_2, y_lim_3, y_max_4, ratio_switch,
              avoid_z=False):
    plt.clf()
    fig=figure(figsize=(9/1.4,15/1.4))

    ########### create plots and subplots ############
    gs = gridspec.GridSpec(3,1,height_ratios=[4,2,4])
    ax1=plt.subplot(gs[0]);  ax2=plt.subplot(gs[1])
    ax3=plt.subplot(gs[2])
    gs.update(hspace=0.0,wspace=0.0,bottom=0.6,top=1.05)

    ax  = [ax1,ax2,ax3]

    for a in [ax1, ax2]:
        a.xaxis.set_visible(False)

    if power == 'lya':
        if dimension == '1d':
            ax1.set_ylabel('$\mathrm{P_F(k) \; [h/Mpc]}$', fontsize=12)
        else:
            ax1.set_ylabel('$\mathrm{3D \; P_F(k) \;[h/Mpc]^3}$', fontsize=12)
        #ax[1].set_ylabel()

        #ax[1].set_ylabel(r'$\mathrm{(\overline{P_T} - \overline{P_P})/\overline{P_T}} \; [\%]$', fontsize=11)
        #ax[2].set_ylabel(r'$\mathrm{(\overline{P_T} - \overline{P_P})/\sigma_{\bar{P_T} - \bar{P_P}}}$', fontsize=11)
        #ax[2].set_ylabel(r'$\mathrm{(\overline{P_T} - \overline{P_P})/\sigma}$', fontsize=11)
        ax3.set_ylabel('$\mathrm{\sigma_S^2/\sigma_P^2}$', fontsize=12)
        ax3.set_xlabel('$\mathrm{k_{\parallel} \; [h/Mpc]}$', fontsize=12)

        ax2.axhline(0, linestyle='--', color='black')
        #ax[1].axhline(-3, linestyle='--', color='black', lw=0.5)
        #ax[1].axhline(3, linestyle='--', color='black', lw=0.5)

        ax3.axhline(1, linestyle='--', color='black')

        if dimension == '1d':
            ax3.set_xlabel('$\mathrm{k_{\parallel} \; [h/Mpc]}$', fontsize=12)

        else:
            ax3.set_xlabel('$\mathrm{k \;[h/Mpc]}$', fontsize=12)

    ax2.set_ylim(-y_lim_2, y_lim_2)
    ax3.set_ylim(0.5, y_max_4)
    ax3.set_yscale('log')

    if power == 'matter':
        ax1.set_ylabel('$\mathrm{3D \; P_M(k) \;[h/Mpc]^3}$', fontsize=12)
        #ax[1].set_ylabel(r'$\mathrm{(\overline{P_T} - \overline{P_P})/\sigma_{\bar{P_T} - \bar{P_P}}}$', fontsize=11)
        ax3.set_ylabel('$\mathrm{\sigma_S^2/\sigma_P^2}$', fontsize=12)
        ax3.set_xlabel('$\mathrm{k_{\parallel} \; [h/Mpc]}$', fontsize=12)

        ax3.set_xlabel('$\mathrm{k \;[h/Mpc]}$', fontsize=12)

    ax2.set_ylabel(r'$\mathrm{\Delta \overline{P}/\overline{P}_S\; [\%]}$')

    ax3.set_yscale('log')
    ax1.set_yscale('log')
    ax1.set_ylim(y_min_1, y_max_1)

    ax2.axhline(0, linestyle='--', color='black')
    #ax[1].axhline(-3, linestyle='--', color='black', lw=0.5)
    #ax[1].axhline(3, linestyle='--', color='black', lw=0.5)

    ax3.axhline(1, linestyle='--', color='black')
    ls = []
    hs = []
    for z, sn, c in zip(zz, snaps, colors):

        ##################################################

        f_out = fout+'.pdf'
        if power == 'matter':
            f1    = f_1 +'_z=%s.txt'%z
            f2    = f_2 +'_z=%s.txt'%z
            f_r   = fr  +'_z=%s.txt'%z
            # read data file
            k1,Pk1,dPk1 = np.loadtxt(f1,unpack=True)
            k2,Pk2,dPk2 = np.loadtxt(f2,unpack=True)
        if power == 'lya':
            k1, Pk1, dPk1, k2, Pk2, dPk2 = create_mean_std(bigbox, dimension, sn)
            _, _, _, kF, PkF, dPkF = create_mean_std(bigbox, dimension, sn, fixed=True)


        #### upper panel ####
        ax1.errorbar(k1, Pk1, yerr=dPk1, lw=2,fmt='o',ms=2, elinewidth=2, capsize=5,
                        linestyle='-',c=c)
        ax1.plot(k2, Pk2, lw=2,
                        linestyle='-',c=c)

        p1 = ax1.errorbar(k1, np.zeros(len(k1))+ 1e5, yerr=1., lw=2, fmt='o', ms=2, elinewidth=2, capsize=5,
                        linestyle='-', c='k')
        p2, = ax1.plot(k1, np.zeros(len(k1))+ 1e5, lw=2,
                        linestyle='-', c='k')

        #### bias panel ####
        pm, = ax2.plot(k1,(Pk1-Pk2)/Pk1*100., c=c,
                 linestyle='-',lw=2, label='z='+str(z))
        sigma = np.sqrt(dPk1**2./nstandard + dPk2**2./npaired)
        print(sigma[-5:])
        if z == 3: ax2.fill_between(k1, -sigma/Pk1*100., sigma/Pk1*100., color='black', alpha=0.2, edgecolor=None, linewidth=0.0)

        #import pdb; pdb.set_trace()
        #pm, = ax2.plot(k1,(Pk2-Pk1)/np.sqrt(dPk1**2/nstandard + dPk2**2/npaired), c=c,
        #         linestyle='-',lw=2, label='z='+str(z))
        ls.append(pm)
        hs.append('z='+str(z))
        #ax2.fill_between([x_min,x_max], -2, 2, color='grey',alpha=0.4)
        #ax2.plot([x_min,x_max],[0,0],linestyle='--',c='k')


        #### variance ratio panel ####
        # y = A/B ---> dy = y*sqrt((dA/A)^2 + (dB/B)^2)
        ratio = (dPk1/dPk2)**2./2.0
        f = interp1d(k1, ratio)
        kpaper = [np.min(k1[~np.isnan(k1) & (k1!=0.0)]), 2.]

        if bigbox: boxsize = 40
        else: boxsize = 20
        if z == 3: print('VarRatio at z={7} for {0} {1} {2} Mpc box at k = {3:0.2f}, {4:0.2f} is {5:0.4f}, {6:0.4f}'.format(power, dimension, boxsize, kpaper[0], kpaper[1], f(kpaper[0]), f(kpaper[1]), z))
        #dratio = ratio*np.sqrt(1.0/standard + 1.0/paired) #old wrong formula
        dratio = 0.5*ratio*(2.0/nstandard + 2.0/npaired) #old wrong formula
        #uncert = 1/2.*np.sqrt(varT/varP)*np.sqrt(2./nTrad + 2./nPaired)
        ax3.plot(k1, ratio, c=c,lw=2)

        if power == 'lya':
            ax3.plot(k1, (dPk1/dPkF)**2., c=c, lw=2, linestyle='--')
            p5,=ax3.plot(k1, np.zeros(len(k1))+1e5, c='k', lw=2, linestyle='--')
        p4,=ax3.plot(k1, np.zeros(len(k1))+1e5, c='k', lw=2)
        ax3.fill_between(k1, ratio+dratio, ratio-dratio, color=c,alpha=0.5)
        if power == 'matter':
            k_r,r,s1,s2 = np.loadtxt(f_r,unpack=True)
            ratio_lazy = np.zeros(len(k_r))
            if ratio_switch:
                for i in xrange(len(k_r)):
                    index = np.where(k1==k_r[i])[0]
                    ratio_lazy[i] = dPk1[index]/(0.5*(s1[i]+s2[i]))
            else:
                ratio_lazy = dPk1/(0.5*(s1+s2))
            ax3.plot(k_r, ratio_lazy, c=c,lw=2, linestyle='--')
            p5,=ax3.plot(k_r, np.zeros(len(k_r))+1e5, c='k', lw=2, linestyle='--')
            indexes = ~np.isnan(ratio)
            ratio_new = ratio[indexes]
            ratio_new = ratio_new[~np.isinf(ratio_new)]
            #print(fout,z,np.max(ratio_lazy),np.max(ratio_new))

        ax3.plot([x_min, x_max],[1.0, 1.0], c='k',linestyle=':')
    for a in ax:
        a.xaxis.set_major_formatter( NullFormatter() )   #unset x label
        a.set_xscale('log')
        a.set_xlim(x_min, x_max)

    # legend and title
    ax1.legend([p1,p2],
               ["%d standard      simulations"%nstandard,
                "%d paired fixed simulations"%npaired],
               loc=0,prop={'size':10.5},ncol=1,frameon=False)
    ax3.legend([p4,p5],
               ["paired fixed",
                "fixed"],
               loc=0,prop={'size':10.5},ncol=2,frameon=False)

    ax2.legend(ls, hs, ncol=3, bbox_to_anchor=(0.6, 1.7))
    fig.savefig(f_out, bbox_inches="tight")
    plt.close(fig)


def lnMeanFlux(z):
    return np.log(0.8)*((1. + z)/3.25)**3.2

def rms(x):
    np.sqrt(np.sum(x**2.)/len(x))

def plotAllPS(power, k, ax, color='black', alpha=0.1, label=None):
    for i, p, in enumerate(power):
        if i == 0: label = label
        else: label = None
        ax.loglog(k, p, color=color, lw=0.5, alpha=alpha, label=label)


def variance(power):

    mean = np.mean(power, axis=0)
    var = np.sum((power - mean)**2., axis=0)/(len(power) - 1.)
    return mean, var




######################### GENERIC ####################
y_min2 = -3.3;  y_max2 = 3.3
x_min  = 0.1;  x_max  = 10
y_lim_2 = 12
######################################################

zs = [4, 3, 2]
snaps = [0, 1, 2]
colors = ['C0', 'C1', 'C2']
"""
make_plots_lya(zz, colors, True, '3d')
make_plots_lya(zz, colors, True, '1d')
make_plots_lya(zz, colors, False, '3d')
make_plots_lya(zz, colors, False, '1d')
"""

###########################################################################
################################ 20 Mpc - hydro ##########################
###########################################################################
nstandard, npaired = 50, 50
root = '20Mpc_hydro'
ratio_switch = False


############################### Pk matter #################################
f_out = '../LyA-paper/Pk_m_20Mpc_hydro' #'%s/Pk_mm_20Mpc_hydro'%root
#f1    = '%s/mean/Pk_m_mean'%root
#f2    = '%s/mean/Pk_m_mean_NCV'%root
f1 = '%s/mean_Pk_m'%root
f2 = '%s/mean_Pk_m_NCV'%root
f_r   = '%s/r_NCV_Pk_m'%root

dimension = None
bigbox = None
power = 'matter'

y_min_1 = 2e-1
y_max_1 = 2e2

#y_lim_2 = 5

y_lim_3 = 5

y_max_4 = 200

make_plot(power, zs, colors, snaps, f_out, f1, f2, f_r, bigbox, dimension,
          nstandard,  npaired, x_min,  x_max, y_min_1, y_max_1, y_lim_2, y_lim_3, y_max_4, ratio_switch)
###########################################################################

############################### Pk Lya 3d #################################
f_out = '../LyA-paper/Pk_lya3d_20Mpc_hydro' #'%s/Pk_mm_20Mpc_hydro'%root
f1    = None
f2    = None
f_r   = None

dimension = '3d'
bigbox = False
power = 'lya'

y_min_1 = 2e-4
y_max_1 = 100

#y_lim_2 = 8

y_lim_3 = 5

y_max_4 = 80


make_plot(power, zs, colors, snaps, f_out, f1, f2, f_r, bigbox, dimension,
          nstandard,  npaired, x_min,  x_max, y_min_1, y_max_1, y_lim_2, y_lim_3, y_max_4, ratio_switch)
###########################################################################

############################### Pk Lya 1d #################################
f_out = '../LyA-paper/Pk_lya1d_20Mpc_hydro' #'%s/Pk_mm_20Mpc_hydro'%root
f1    = None
f2    = None
f_r   = None

dimension = '1d'
bigbox = False
power = 'lya'

y_min_1 = 2e-4
y_max_1 = 1

#y_lim_2 = 2

y_lim_3 = 5

y_max_4 = 30


make_plot(power, zs, colors, snaps, f_out, f1, f2, f_r, bigbox, dimension,
          nstandard,  npaired, x_min,  x_max, y_min_1, y_max_1, y_lim_2, y_lim_3, y_max_4, ratio_switch)
###########################################################################


###########################################################################
################################ 40 Mpc - hydro ##########################
###########################################################################
nstandard, npaired = 50, 25
root = '40Mpc_hydro'
ratio_switch = False


############################### Pk matter #################################
f_out = '../LyA-paper/Pk_m_40Mpc_hydro' #'%s/Pk_mm_20Mpc_hydro'%root
#f1    = '%s/mean/Pk_m_mean'%root
#f2    = '%s/mean/Pk_m_mean_NCV'%root

f1 = '%s/mean_Pk_m'%root
f2 = '%s/mean_Pk_m_NCV'%root
f_r   = '%s/r_NCV_Pk_m'%root

dimension = None
bigbox = None
power = 'matter'

y_min_1 = 2e-1
y_max_1 = 5e2

#y_lim_2 = 3

y_lim_3 = 2

y_max_4 = 400


make_plot(power, zs, colors, snaps, f_out, f1, f2, f_r, bigbox, dimension,
          nstandard,  npaired, x_min,  x_max, y_min_1, y_max_1, y_lim_2, y_lim_3, y_max_4, ratio_switch)
###########################################################################

############################### Pk Lya 3d #################################
f_out = '../LyA-paper/Pk_lya3d_40Mpc_hydro' #'%s/Pk_mm_20Mpc_hydro'%root
f1    = None
f2    = None
f_r   = None

dimension = '3d'
bigbox = True
power = 'lya'
y_min_1 = 2e-4
y_max_1 = 100


#y_lim_2 = 10

y_lim_3 = 3

y_max_4 = 150


make_plot(power, zs, colors, snaps, f_out, f1, f2, f_r, bigbox, dimension,
          nstandard,  npaired, x_min,  x_max, y_min_1, y_max_1, y_lim_2, y_lim_3, y_max_4, ratio_switch)
###########################################################################

############################### Pk Lya 1d #################################
f_out = '../LyA-paper/Pk_lya1d_40Mpc_hydro' #'%s/Pk_mm_20Mpc_hydro'%root
f1    = None
f2    = None
f_r   = None

dimension = '1d'
bigbox = True
power = 'lya'

y_min_1 = 5e-4
y_max_1 = 1

#y_lim_2 = 2

y_lim_3 = 3

y_max_4 = 80


make_plot(power, zs, colors, snaps, f_out, f1, f2, f_r, bigbox, dimension,
          nstandard,  npaired, x_min,  x_max, y_min_1, y_max_1, y_lim_2, y_lim_3, y_max_4, ratio_switch)
###########################################################################
