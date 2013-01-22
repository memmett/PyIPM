"""Alberta boreal SW and AW competition experiment.

Generate plots to compare with/with out competition in the mortality
model.

"""


import numpy as np
import cPickle as pickle
import pylab as plt

from exp_abb_kernels import abb_init_kernels

from utils.stats import ks1



def comp_ks_plots(plots1, plots2, flavour):

    for plotname in plots1:

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

                fig, ax = plt.subplots(2, 2)

                meassw = np.asarray(sw1.meas['SW'][t])
                measaw = np.asarray(aw1.meas['AW'][t])

                ytop = max(nsw1[j].max(), nsw2[j].max(), naw1[j].max(), naw2[j].max()) * 1.2

                # spruce
                weights = 1e4 / sw1.plot_size * (dx/4) * np.ones(meassw.shape)
                ax[0,0].hist(meassw, bins=N/4, weights=weights,
                             histtype='stepfilled', color='0.5', edgecolor='0.5')

                ax[0,0].plot(sw1.x, nsw1[j], color='r', linestyle='-', linewidth=1)
                ax[0,0].plot(sw2.x, nsw2[j], color='b', linestyle='-', linewidth=1)

                # ax[0,0].legend(loc='best')
                ax[0,0].set_ylabel('stems per hectare')
                ax[0,0].set_title('spruce')

                ax[0,0].set_ylim(0, ytop)
                plt.setp(ax[0,0].get_xticklabels(), visible=False)

                # aspen
                weights = 1e4 / aw1.plot_size * (dx/4) * np.ones(measaw.shape)
                ax[0,1].hist(measaw, bins=N/4, weights=weights,
                             histtype='stepfilled', color='0.5', edgecolor='0.5')

                ax[0,1].plot(aw1.x, naw1[j], color='r', linestyle='-', linewidth=1)
                ax[0,1].plot(aw2.x, naw2[j], color='b', linestyle='-', linewidth=1)

                ax[0,1].set_ylim(0, ytop)
                ax[0,1].set_title('aspen')

                plt.setp(ax[0,1].get_xticklabels(), visible=False)
                plt.setp(ax[0,1].get_yticklabels(), visible=False)

                # spruce ks
                ccdf, dcdf, d, p = ks1(sw1.meas['SW'][t], nsw1[j], x)
                z = np.sort(sw1.meas['SW'][t])
                ax[1,0].plot(z, dcdf, color='0.5', linestyle='--', linewidth=2)
                ax[1,0].plot(z, ccdf, '-r')

                ccdf, dcdf, d, p = ks1(sw1.meas['SW'][t], nsw2[j], x)
                ax[1,0].plot(z, ccdf, '-b')

                ax[1,0].set_xlabel('dbh')
                ax[1,0].set_ylabel('cdf')

                ax[1,0].set_xlim(L, U)
                ax[1,0].set_ylim(0, 1)

                # aspen ks
                ccdf, dcdf, d, p = ks1(aw1.meas['AW'][t], naw1[j], x)
                z = np.sort(aw1.meas['AW'][t])

                ax[1,1].plot(z, dcdf, color='0.5', linestyle='--', linewidth=2, label='meas.')
                ax[1,1].plot(z, ccdf, '-r', label='w/ comp.')

                ccdf, dcdf, d, p = ks1(aw1.meas['AW'][t], naw2[j], x)

                ax[1,1].plot(z, ccdf, '-b', label='w/o comp.')

                ax[1,1].set_xlabel('dbh')
                ax[1,1].set_ylabel('cdf')

                ax[1,1].set_xlim(L, U)
                ax[1,1].set_ylim(0, 1)

                ax[1,1].legend(loc='lower right', prop={'size': 10})

                plt.setp(ax[1,1].get_yticklabels(), visible=False)


                fig.savefig('plots/projkscomp_%s_%d.png' % (plotname, t))
                # plt.show()
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

