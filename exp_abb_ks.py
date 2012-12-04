"""Alberta boreal SW and AW competition experiment."""


import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

import cPickle as pickle

import matplotlib.pylab as plt

from utils.stats import ks1


###############################################################################
# load results

with open('abb.pkl', 'r') as f:
    plots = pickle.load(f)


plot_cdf = True

###############################################################################
# main

pvals = []

for plotname in plots:

    plot = plots[plotname]

    # init

    sw = ABBSW()
    aw = ABBAW()

    L, U, N, T, x, plotfile = plot['attrs']
    nsw, naw = plot['nsw'], plot['naw']

    # setup kernels and connect them to each other
    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(MidPoint(), N)
        k.measurements(plotfile, plotname)


    for j, t in enumerate(T):
        if t in sw.years:
            if j > 0:
                ccdf, dcdf, d, p = ks1(sw.meas['SW'][t], nsw[j], x)
                # print plotname, t, p
                pvals.append(p)

                if plot_cdf:
                  plt.figure()
                  z = np.sort(sw.meas['SW'][t])
                  plt.plot(z, ccdf, '-k', label='model')
                  plt.plot(z, dcdf, '-r', label='meas')
                  plt.xlabel('dbh')
                  plt.ylabel('cdf')
                  plt.title('plot: %s, year: %d' % (plotname, t))
                  plt.legend(loc='best')
                  plt.savefig('plots/%s_ks_%d.png' % (plotname, t))
                  plt.close()

        # if t in aw.years:
        #     plt.hist(aw.meas['AW'][t], bins=N/4,
        #              histtype='stepfilled', color='red', alpha=0.6, label='aspen meas')


plt.figure()
x = range(len(pvals))
plt.plot(x, sorted(pvals))
plt.plot(x, len(x)*[0.05], '-r')
plt.title('ks pval')
plt.xlabel('plot')
plt.ylabel('pval')
plt.savefig('plots/ks.png')
plt.close()

plt.figure()
plt.semilogy(x, sorted(pvals))
plt.semilogy(x, len(x)*[0.05], '-r')
plt.title('ks pval')
plt.xlabel('plot')
plt.ylabel('log(pval)')
plt.savefig('plots/ks_log.png')
plt.close()
