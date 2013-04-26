"""Alberta boreal SW and AW competition experiment.

Generate plots to compare with/with out competition in the mortality
model.

"""

import matplotlib
matplotlib.rcParams['font.size'] = 24
matplotlib.rcParams['legend.fontsize'] = 18


import numpy as np
import cPickle as pickle
import pylab as plt

from exp_abb_kernels import abb_init_kernels

from utils.stats import ks1


def comp_ks_plots(plots1, plots2, flavour):

    for plotname in plots1:

        if plotname != '237':
            continue

        plot1 = plots1[plotname]
        plot2 = plots2[plotname]

        L, U, N, T, x, plotfile = plot1['attrs']
        dx = float(U - L) / N

        nsw1, naw1 = plot1['nsw'], plot1['naw']
        nsw2, naw2 = plot2['nsw'], plot2['naw']

        sw1, aw1 = abb_init_kernels(L, U, N, flavour, plotname)
        sw2, aw2 = abb_init_kernels(L, U, N, flavour, plotname)

        for j, t in enumerate(T):
            if t in sw1.years and j > 0:

                if t != 2000:
                    continue

                #fig, ax = plt.subplots(ncols=2)
                plt.figure(1)

                meassw = np.asarray(sw1.meas['SW'][t])
                # measaw = np.asarray(aw1.meas['AW'][t])

                # ytop = max(nsw1[j].max(), nsw2[j].max(), naw1[j].max(), naw2[j].max()) * 1.2

                # spruce
                weights = 1e4 / sw1.plot_size * (dx/4) * np.ones(meassw.shape)
                plt.hist(meassw, bins=N/4, weights=weights,
                             histtype='stepfilled', color='0.5', edgecolor='0.5')

                plt.plot(sw1.x, nsw1[j], color='r', linestyle='-', linewidth=1)
                plt.plot(sw2.x, nsw2[j], color='b', linestyle='-', linewidth=1)

                plt.ylabel('stems per hectare')
                plt.xlabel('size')
                plt.xlim(L, 500)

                plt.savefig('plots/projkscompswonly_%s_%d_sph.png' % (plotname, t))

                # spruce ks
                ccdf, dcdf, d, p = ks1(sw1.meas['SW'][t], nsw1[j], x)
                z = np.sort(sw1.meas['SW'][t])

                plt.figure(2)

                plt.plot(z, dcdf, color='0.5', linestyle='-', linewidth=2)
                plt.plot(z, ccdf, '-r')

                ccdf, dcdf, d, p = ks1(sw1.meas['SW'][t], nsw2[j], x)
                plt.plot(z, ccdf, '-b')

                plt.xlabel('size')
                plt.ylabel('cumulative distribution function')

                plt.legend(['measurements', 'with comp', 'without comp'], loc='lower right')

                plt.xlim(L, 500)
                plt.ylim(0, 1)

                plt.savefig('plots/projkscompswonly_%s_%d_cdf.png' % (plotname, t))
                plt.close()




        # psw1 = [ sw1.population(nsw1[j]) for j in range(len(T)) ]
        # psw2 = [ sw2.population(nsw2[j]) for j in range(len(T)) ]

        # paw1 = [ aw1.population(naw1[j]) for j in range(len(T)) ]

        # Tmeas   = sw1.years
        # pswmeas = np.asarray([ len(sw1.meas['SW'][t]) for t in Tmeas ]) / sw1.plot_size * 1e4
        # pawmeas = np.asarray([ len(aw1.meas['AW'][t]) for t in Tmeas ]) / aw1.plot_size * 1e4

        # plt.figure()
        # plt.plot(T[1:], psw1[1:], 'ob', label='spruce (comp)')
        # plt.plot(T[1:], psw2[1:], 'sc', label='spruce (no comp)')

        # plt.plot(T[1:], paw1[1:], 'or', label='aspen')

        # plt.plot(Tmeas, pswmeas, 'vk', label='spruce (meas)')
        # plt.plot(Tmeas, pawmeas, '^k', label='aspen (meas)')
        # plt.legend(loc='best')
        # plt.xlabel('year')
        # plt.ylabel('population')
        # plt.title('population vs time ' + plotname)
        # plt.savefig('plots/popcomp_%s.png' % plotname)
        # plt.close()


if __name__ == '__main__':

    with open('out/abb_with_comp.pkl', 'r') as f:
        plots1 = pickle.load(f)

    with open('out/abb_with_nocomp.pkl', 'r') as f:
        plots2 = pickle.load(f)

    comp_ks_plots(plots1, plots2, 'with_comp')

