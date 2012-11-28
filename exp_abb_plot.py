"""Alberta boreal SW and AW competition experiment."""


import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

import cPickle as pickle

import matplotlib.pylab as plt

from ks_1samp import ks1


###############################################################################
# load results

with open('abb.pkl', 'r') as f:
    plots = pickle.load(f)


###############################################################################
# main

for plotname in plots:

    plot = plots[plotname]

    # init

    method = MidPoint()

    sw = ABBSW()
    aw = ABBAW()

    L, U, N, T, x, plotfile = plot['attrs']
    nsw, naw = plot['nsw'], plot['naw']

    plotname = plotfile.split('/')[-1].split('.')[0]

    # setup kernels and connect them to each other
    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(method, N)
        k.measurements(plotfile, plotname)


    for j, t in enumerate(T):
        plt.figure()
        plt.subplot(211)

        if t in sw.years:
            if j > 0:
                print ks1(sw.meas['SW'][t], nsw[j], x, 10000)
            plt.hist(sw.meas['SW'][t], bins=N/4,
                     histtype='stepfilled', color='blue', alpha=0.6, label='spruce meas')

        if j > 0:
            plt.plot(sw.x, nsw[j], '-b', label='spruce', linewidth=2)

        plt.legend(loc='best')
        plt.ylabel('stems per hectare')
        plt.title('year ' + str(t))

        plt.subplot(212)

        if t in aw.years:
            plt.hist(aw.meas['AW'][t], bins=N/4,
                     histtype='stepfilled', color='red', alpha=0.6, label='aspen meas')

        if j > 0:
            plt.plot(aw.x, naw[j], '-r', label='aspen', linewidth=2)

        plt.legend(loc='best')
        plt.xlabel('dbh (mm)')
        plt.ylabel('stems per hectare')
        plt.savefig('plots/%s_projection_%02d.png' % (plotname, j))



    psw = [ sw.population(nsw[j]) for j in range(len(T)) ]
    paw = [ aw.population(naw[j]) for j in range(len(T)) ]

    psw[0] = sum(nsw[0])
    paw[0] = sum(naw[0])

    Tmeas   = [ t for t in sw.meas['SW'] ]
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
    plt.title('population vs time')
    plt.savefig('plots/%s_population.png' % plotname)

plt.show()
