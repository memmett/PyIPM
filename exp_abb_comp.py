"""Alberta boreal SW and AW competition experiment."""


import numpy as np
import cPickle as pickle
import pylab as plt

from exp_abb_kernels import abb_init_kernels


###############################################################################
# load results

dsets = [ 'abb_with_nocomp.pkl', 'abb_with_comp.pkl' ]

with open(dsets[0], 'r') as f:
    plots1 = pickle.load(f)

with open(dsets[1], 'r') as f:
    plots2 = pickle.load(f)


###############################################################################
# main

for plotname in plots_comp:

    plot1 = plots1[plotname]
    plot2 = plots2[plotname]

    L, U, N, T, x, plotfile = plot1['attrs']
    dx = float(U - L) / N

    nsw1, naw1 = plot1['nsw'], plot1['naw']
    nsw2, naw2 = plot2['nsw'], plot2['naw']

    sw1, aw1 = abb_init_kernels(L, U, N, 'with_comp',   plotname)
    sw2, aw2 = abb_init_kernels(L, U, N, 'with_nocomp', plotname)

    for j, t in enumerate(T):
        if t in sw1.years and j > 0:
            plt.figure()

            # spruce
            plt.subplot(211)
            meas    = np.asarray(sw1.meas['SW'][t])
            weights = 1e4 / sw1.plot_size * (dx/4) * np.ones(meas.shape)
            plt.hist(meas, bins=N/4, weights=weights,
                     histtype='stepfilled', color='blue', alpha=0.6, label='spruce meas')

            plt.plot(sw1.x, nsw1[j], '-b', label='spruce (comp)',    linewidth=2)
            plt.plot(sw2.x, nsw2[j], '-c', label='spruce (no comp)', linewidth=2)

            plt.legend(loc='best')
            plt.ylabel('stems per hectare')
            plt.title(plotname + ' year ' + str(t))

            plt.subplot(212)

            meas    = np.asarray(aw1.meas['AW'][t])
            weights = 1e4 / aw1.plot_size * (dx/4) * np.ones(meas.shape)
            plt.hist(meas, bins=N/4, weights=weights,
                     histtype='stepfilled', color='red', alpha=0.6, label='aspen meas')

            plt.plot(aw1.x, naw1[j], '-r', label='aspen (comp)',    linewidth=2)
            plt.plot(aw2.x, naw2[j], '-m', label='aspen (no comp)', linewidth=2)

            plt.legend(loc='best')
            plt.xlabel('dbh (mm)')
            plt.ylabel('stems per hectare')
            plt.savefig('plots/%s_projection_comp_%d.png' % (plotname, t))

            plt.close()


    psw1 = [ sw1.population(nsw1[j]) for j in range(len(T)) ]
    psw2 = [ sw2.population(nsw2[j]) for j in range(len(T)) ]

    paw1 = [ aw1.population(naw1[j]) for j in range(len(T)) ]

    Tmeas   = sw1.years
    pswmeas = np.asarray([ len(sw1.meas['SW'][t]) for t in Tmeas ]) / sw1.plot_size * 1e4
    pawmeas = np.asarray([ len(aw1.meas['AW'][t]) for t in Tmeas ]) / aw1.plot_size * 1e4

    plt.figure()
    plt.plot(T[1:], psw1[1:], 'ob', label='spruce (comp)')
    plt.plot(T[1:], psw2[1:], 'sc', label='spruce (no comp)')

    plt.plot(T[1:], paw1[1:], 'or', label='aspen')

    plt.plot(Tmeas, pswmeas, 'vk', label='spruce (meas)')
    plt.plot(Tmeas, pawmeas, '^k', label='aspen (meas)')
    plt.legend(loc='best')
    plt.xlabel('year')
    plt.ylabel('population')
    plt.title('population vs time ' + plotname)
    plt.savefig('plots/%s_population_comp.png' % plotname)
    plt.close()
