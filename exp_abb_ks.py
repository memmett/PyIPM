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

    # setup kernels and connect them to each other
    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(method, N)
        k.measurements(plotfile)


    for j, t in enumerate(T):
        if t in sw.years:
            if j > 0:
                plt.figure()
                ccdf, dcdf, d, p = ks1(sw.meas['SW'][t], nsw[j], x)
                print d, p
                z = np.sort(sw.meas['SW'][t])
                plt.plot(z, ccdf, '-k')
                plt.plot(z, dcdf, '-r')

        # if t in aw.years:
        #     plt.hist(aw.meas['AW'][t], bins=N/4,
        #              histtype='stepfilled', color='red', alpha=0.6, label='aspen meas')


plt.show()
