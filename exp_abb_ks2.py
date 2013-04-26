"""Alberta boreal SW and AW competition experiment.

Generate all KS related plots for the with_comp and with_nocomp
flavours.

"""

import matplotlib
matplotlib.rcParams['font.size'] = 14


import numpy as np
import cPickle as pickle
import matplotlib.pylab as plt

from exp_abb_kernels import abb_init_kernels, ignore
from utils.stats import ks1


##############################################################################
# ks plots

def ks_plots(plots, flavour, species):

    pvals = []
    for plotname in plots:

        if plotname in ignore:
            continue

        plot = plots[plotname]

        L, U, N, T, x, plotfile = plot['attrs']
        sw, aw = abb_init_kernels(L, U, N, flavour, plotname)

        nsw, naw = plot['nsw'], plot['naw']

        for j, t in enumerate(T):
            if t in sw.years:
                if j > 0:
                    if species == 'sw':
                        ccdf, dcdf, d, p = ks1(sw.meas['SW'][t], nsw[j], x)
                        z = np.sort(sw.meas['SW'][t])
                    elif species == 'aw':
                        ccdf, dcdf, d, p = ks1(aw.meas['AW'][t], naw[j], x)
                        z = np.sort(aw.meas['AW'][t])
                    pvals.append(p)


    # # histogram of pvals
    # plt.figure()
    # plt.hist(pvals, bins=21)
    # plt.xlabel('p-value')
    # plt.ylabel('frequency')
    # plt.savefig('plots/kshist_%s_%s.png' % (flavour, species))
    # plt.close()


    # sorted scatter plot
    plt.figure()
    x = range(len(pvals))
    plt.plot(x, sorted(pvals), '.k')
    plt.plot(x, len(x)*[0.05], '-r')
    plt.xlabel('index')
    plt.ylabel('p-value')
    plt.savefig('plots/ksscatter_%s_%s.png' % (flavour, species))
    plt.close()

    pvals = np.asarray(pvals)

    n    = len(pvals)
    nsig = np.count_nonzero(pvals > 0.05)
    print flavour, species, n, nsig, n - nsig

    return pvals


def ks_compare(plots1, plots2, flavour, species):

    pvals1 = []
    pvals2 = []

    for plotname in plots1:

        if plotname in ignore:
            continue

        plot1 = plots1[plotname]
        plot2 = plots2[plotname]

        L, U, N, T, x, plotfile = plot1['attrs']
        sw, aw = abb_init_kernels(L, U, N, flavour, plotname)

        nsw1, naw1 = plot1['nsw'], plot1['naw']
        nsw2, naw2 = plot2['nsw'], plot2['naw']

        for j, t in enumerate(T):
            if t in sw.years:
                if j > 0:
                    if species == 'sw':
                        ccdf, dcdf, d, p = ks1(sw.meas['SW'][t], nsw1[j], x)
                        # pvals1.append(p)
                        pvals1.append(d)

                        ccdf, dcdf, d, p = ks1(sw.meas['SW'][t], nsw2[j], x)
                        # pvals2.append(p)
                        pvals2.append(d)

                    elif species == 'aw':
                        ccdf, dcdf, d, p = ks1(aw.meas['AW'][t], naw1[j], x)
                        # pvals1.append(p)
                        pvals1.append(d)

                        ccdf, dcdf, d, p = ks1(aw.meas['AW'][t], naw2[j], x)
                        # pvals2.append(p)
                        pvals2.append(d)

    # pmax = 1.0
    pmax = 0.6

    # p-val vs p-val (actually d-val vs d-val for now...)
    plt.figure()
    plt.plot(pvals1, pvals2, 'ok', alpha=0.5)
    # plt.axhline(y=0.05, color='black', linestyle='--')
    # plt.axvline(x=0.05, color='black', linestyle='--')
    plt.plot([0.0, pmax], [0.0, pmax], color='black')
    # plt.xlabel('p-value (without competition)')
    # plt.ylabel('p-value (with competition)')
    plt.xlabel('Kolmogorov-Smirnov statistic (without competition)')
    plt.ylabel('Kolmogorov-Smirnov statistic (with competition)')

    dvals1 = np.asarray(pvals1)
    dvals2 = np.asarray(pvals2)

    print 'species', species
    print 'npoints', len(dvals1)
    # print 'nbetter', np.count_nonzero(dvals1 > dvals2)

    plt.savefig('plots/ksscatter_%s_comp.png' % species)
    plt.close()


##############################################################################
# main

if __name__ == '__main__':

    with open('out/abb_with_nocomp.pkl', 'r') as f:
        plots1 = pickle.load(f)

    # pvals = ks_plots(plots1, 'with_nocomp', 'aw')
    # pvals = ks_plots(plots1, 'with_nocomp', 'sw')

    with open('out/abb_with_comp.pkl', 'r') as f:
        plots2 = pickle.load(f)

    # pvals = ks_plots(plots2, 'with_comp', 'aw')
    # pvals = ks_plots(plots2, 'with_comp', 'sw')

    ks_compare(plots1, plots2, 'with_comp', 'sw')
    ks_compare(plots1, plots2, 'with_comp', 'aw')
