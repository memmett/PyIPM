"""Alberta boreal PSP statistics."""


import numpy as np
import cPickle as pickle
import matplotlib.pylab as plt

from exp_abb_kernels import ignore
from kernels.abb import ABBSW, ABBAW
from utils.stats import ks1


with open('out/abb_with_comp.pkl', 'r') as f:
    plots = pickle.load(f)


npsp  = 0
nmeas = 0
nyears = []
dbh  = { 'sw': [], 'aw': [] }
sph  = { 'sw': [], 'aw': [], 'tot': [] }

for plotname in plots:
    if plotname in ignore:
        continue

    k = ABBSW()
    k.measurements('kernels/abb/%s.csv' % plotname, plotname)
    for year in k.years:
        dbh['sw'].extend(k.meas['SW'][year])
        dbh['aw'].extend(k.meas['AW'][year])

    sph['sw'].append(len(k.meas['SW'][k.first_year]) / k.plot_size * 1e4)
    sph['aw'].append(len(k.meas['AW'][k.first_year]) / k.plot_size * 1e4)
    sph['tot'].append(sph['sw'][-1] + sph['aw'][-1])

    years = np.asarray(k.years)
    nyears.extend(list(years[1:] - years[:-1]))

    nmeas = nmeas + len(k.years)
    npsp  = npsp + 1

print 'npsp', npsp
print 'nmeas', nmeas
print 'avg years', np.mean(nyears), np.std(nyears)

import pylab as plt

# dbh histograms
fig, ax = plt.subplots(1, 2, sharey=True)

ax[0].hist(dbh['sw'], 50, color='0.5')
ax[0].set_title('Spruce')
ax[0].set_ylabel('frequency')
ax[0].set_xlabel('dbh (mm)')

ax[1].hist(dbh['aw'], 50, color='0.5')
ax[1].set_title('Aspen')
ax[1].set_xlabel('dbh (mm)')

fig.savefig('plots/pspdbhhist.png')

# sph plots
plt.figure()
sph['sw'] = np.asarray(sph['sw'])
sph['aw'] = np.asarray(sph['aw'])

idx = np.argsort(sph['sw'])
x   = np.arange(len(idx))

width = 0.35

plt.bar(x, sph['sw'][idx], width, color='black', label='spruce')
plt.bar(x+width, sph['aw'][idx], width, color='white', hatch='/', label='aspen')

plt.xlabel('psp')
plt.ylabel('sph (#/ha)')
plt.legend(loc='best')

plt.savefig('plots/pspsph.png')

plt.show()



