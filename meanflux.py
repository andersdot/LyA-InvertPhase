import glob
import numpy as np
import matplotlib.pyplot as plt
z = [4, 3, 2]
shifts = [-0.1, 0.1]
filenames = glob.glob('spec3d400*npz')
filenames.sort()

meanflux = {}
varflux = {}

fig, ax = plt.subplots(1, 3, figsize=(6, 3))

for file in filenames:
    bits = file.split('_')
    print(len(bits))
    if len(bits) == 5:
        label = bits[1] + bits[2]
        redshift = z[int(bits[3])]
        if bits[2] == '0': shift = -0.1
        if bits[2] == '1': shift =  0.1
    if len(bits) == 3:
        shift = 0.0
        label = 'S'
        redshift = z[int(bits[1])]
    data = np.load(file)
    meanflux[label+str(redshift)] = np.mean(data['meanflux'])
    varflux[label+str(redshift)] = np.std(data['meanflux'])

    ax[0].errorbar(redshift+shift, np.mean(data['meanflux']), yerr=3.*np.std(data['meanflux']), label=label, linewidth=3, color='k')

for axes, dict in zip(ax[1:3], [meanflux, varflux]):
    for redshift in [2, 3, 4]:
        ratio0 = dict['S' + str(redshift)]/dict['NCV0' + str(redshift)]
        ratio1 = dict['S' + str(redshift)]/dict['NCV1' + str(redshift)]
        print(ratio0, ratio1, redshift)
        axes.scatter(redshift-0.1, ratio0, label='NCV0', color='k')
        axes.scatter(redshift+0.1, ratio1, label='NCV1', color='k')

for axes in ax:
    axes.set_xlabel('z')

ax[0].set_ylabel('<F>')
ax[1].set_ylabel(r'$\overline{<F>}_S/\overline{<F>}_{FP}$')
ax[2].set_ylabel(r'$\sigma_{<F>S}/\sigma_{<F>FP}$')
print(meanflux, varflux)
#plt.legend()
plt.tight_layout()
plt.savefig('meanfluxNCV.pdf')
