"""Alberta boreal SW and AW competition experiment."""


import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

import cPickle as pickle

import matplotlib.pylab as plt


###############################################################################
# load results

with open('out/abb_with_comp.pkl', 'r') as f:
    plots = pickle.load(f)


###############################################################################
# main

for plotname in plots:

    plot = plots[plotname]

    # init

    sw = ABBSW()
    aw = ABBAW()

    L, U, N, T, x, plotfile = plot['attrs']
    nsw, naw = plot['nsw'], plot['naw']

    dx = float(U - L) / N

    plotname = plotfile.split('/')[-1].split('.')[0]

    # setup kernels and connect them to each other
    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(MidPoint(), N)
        k.measurements(plotfile, plotname)


    for j, t in enumerate(T):
        if t in sw.years and j > 0:
            plt.figure()

            plt.subplot(211)

            meas    = np.asarray(sw.meas['SW'][t])
            weights = 1e4 / sw.plot_size * (dx/4) * np.ones(meas.shape)
            plt.hist(meas, bins=N/4, weights=weights,
                     histtype='stepfilled', color='0.5', edgecolor='0.5', label='spruce meas')

            plt.plot(sw.x, nsw[j], '-b', label='spruce', linewidth=2)

            plt.legend(loc='best')
            plt.ylabel('stems per hectare')
            plt.title(plotname + ' year ' + str(t))

            plt.subplot(212)

            meas    = np.asarray(sw.meas['AW'][t])
            weights = 1e4 / sw.plot_size * (dx/4) * np.ones(meas.shape)
            plt.hist(meas, bins=N/4, weights=weights,
                     histtype='stepfilled', color='0.5', edgecolor='0.5', label='aspen meas')

            plt.plot(aw.x, naw[j], '-r', label='aspen', linewidth=2)

            plt.legend(loc='best')
            plt.xlabel('dbh (mm)')
            plt.ylabel('stems per hectare')
            plt.savefig('plots/proj_%s_%d.png' % (plotname, t))

            plt.close()


    psw = [ sw.population(nsw[j]) for j in range(len(T)) ]
    paw = [ aw.population(naw[j]) for j in range(len(T)) ]

    psw[0] = sum(nsw[0])
    paw[0] = sum(naw[0])

    # Tmeas   = [ t for t in sw.meas['SW'] ]
    Tmeas   = sw.years
    pswmeas = np.asarray([ len(sw.meas['SW'][t]) for t in Tmeas ]) / sw.plot_size * 1e4
    pawmeas = np.asarray([ len(aw.meas['AW'][t]) for t in Tmeas ]) / aw.plot_size * 1e4

    plt.figure()
    plt.plot(T[1:], psw[1:], 'ob', label='spruce')
    plt.plot(T[1:], paw[1:], 'or', label='aspen')
    plt.plot(Tmeas, pswmeas, 'sc', label='spruce (meas)')
    plt.plot(Tmeas, pawmeas, 'sm', label='aspen (meas)')
    plt.legend(loc='best')
    plt.xlabel('year')
    plt.ylabel('population')
    plt.title('population vs time ' + plotname)
    plt.savefig('plots/pop_%s.png' % plotname)
    plt.close()

#plt.show()

