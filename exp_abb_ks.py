"""Alberta boreal SW and AW competition experiment.

Generate all KS related plots for the with_comp and with_nocomp
flavours.

"""


import numpy as np
import cPickle as pickle
import matplotlib.pylab as plt

from exp_abb_kernels import abb_init_kernels
from utils.stats import ks1


##############################################################################
# ks plots

def ks_plots(plots, flavour):

    pvals = []

    for plotname in plots:

        plot = plots[plotname]

        L, U, N, T, x, plotfile = plot['attrs']
        sw, aw = abb_init_kernels(L, U, N, flavour, plotname)

        nsw, naw = plot['nsw'], plot['naw']

        for j, t in enumerate(T):
            if t in sw.years:
                if j > 0:
                    ccdf, dcdf, d, p = ks1(sw.meas['SW'][t], nsw[j], x)
                    # print plotname, t, p
                    pvals.append(p)

                    plt.figure()
                    z = np.sort(sw.meas['SW'][t])
                    plt.plot(z, ccdf, '-k', label='model')
                    plt.plot(z, dcdf, '-r', label='meas')
                    plt.xlabel('dbh')
                    plt.ylabel('cdf')
                    plt.title('plot: %s, year: %d' % (plotname, t))
                    plt.legend(loc='best')
                    plt.savefig('plots/kscdf_%s_%s_%d.png' % (flavour, plotname, t))
                    plt.close()


    # histogram of pvals
    plt.figure()
    plt.hist(pvals, bins=21)
    plt.xlabel('p-value')
    plt.ylabel('frequency')
    plt.savefig('plots/kshist_%s.png' % flavour)
    plt.close()

    # sorted scatter plot
    plt.figure()
    x = range(len(pvals))
    plt.plot(x, sorted(pvals), '.k')
    plt.plot(x, len(x)*[0.05], '-r')
    plt.xlabel('index')
    plt.ylabel('p-value')
    plt.savefig('plots/ksscatter_%s.png' % flavour)
    plt.close()

    # plt.figure()
    # plt.semilogy(x, sorted(pvals))
    # plt.semilogy(x, len(x)*[0.05], '-r')
    # plt.title('ks pval')
    # plt.xlabel('plot')
    # plt.ylabel('log(pval)')
    # plt.savefig('plots/ks_log.png')
    # plt.close()

    return pvals


##############################################################################
# main

if __name__ == '__main__':

    with open('out/abb_with_nocomp.pkl', 'r') as f:
        plots = pickle.load(f)

    pvals = ks_plots(plots, 'with_nocomp')


    with open('out/abb_with_comp.pkl', 'r') as f:
        plots = pickle.load(f)

    pvals = ks_plots(plots, 'with_comp')
