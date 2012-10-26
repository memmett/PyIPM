"""Alberta boreal SW and AW competition experiment."""


import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW
from kernels.abb import first_year, last_year, measurements

import matplotlib.pylab as plt


###############################################################################
# load results

globals().update(np.load('abb.npz'))


###############################################################################
# init

method = MidPoint()

sw = ABBSW()
aw = ABBAW()

# setup kernels and connect them to each other
for k in [ sw, aw ]:
    k.L, k.U = L, U
    k.setup(method, N)


###############################################################################
# main

for j, t in enumerate(T):
    plt.figure()
    plt.subplot(211)

    if t in measurements['SW']:
        plt.hist(measurements['SW'][t], bins=N/4,
                 histtype='stepfilled', color='blue', alpha=0.6, label='spruce meas')

    plt.plot(sw.x, nsw[j], '-b', label='spruce', linewidth=2)

    plt.legend(loc='best')
    plt.ylabel('stems per hectare')
    plt.title('year ' + str(t))

    plt.subplot(212)

    if t in measurements['AW']:
        plt.hist(measurements['AW'][t], bins=N/4,
                 histtype='stepfilled', color='red', alpha=0.6, label='aspen meas')

    plt.plot(aw.x, naw[j], '-r', label='aspen', linewidth=2)

    plt.legend(loc='best')
    plt.xlabel('dbh (mm)')
    plt.ylabel('stems per hectare')
    plt.savefig('plots/projection_%02d.pdf' % (j))



psw = [ sw.population(nsw[j]) for j in range(len(T)) ]
paw = [ aw.population(naw[j]) for j in range(len(T)) ]

psw[0] = sum(nsw[0])
paw[0] = sum(naw[0])

Tmeas   = [ t for t in measurements['SW'] ]
pswmeas = np.asarray([ len(measurements['SW'][t]) for t in Tmeas ]) / sw.plot_size * 1e4
pawmeas = np.asarray([ len(measurements['AW'][t]) for t in Tmeas ]) / aw.plot_size * 1e4

plt.figure()
plt.plot(T, psw, '-ob', label='spruce')
plt.plot(T, paw, '-or', label='aspen')
plt.plot(Tmeas, pswmeas, 'sc', label='spruce (meas)')
plt.plot(Tmeas, pawmeas, 'sm', label='aspen (meas)')
plt.legend(loc='best')
plt.xlabel('year')
plt.ylabel('population')
plt.title('population vs time')
plt.savefig('plots/population.pdf')

plt.show()
