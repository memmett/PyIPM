"""Alberta boreal SW and AW competition experiment."""


import numpy as np
import cPickle as pickle
import matplotlib.pylab as plt

from exp_abb_kernels import abb_init_kernels


###############################################################################
# load results

with open('out/abb_eviction.pkl', 'r') as f:
    plots = pickle.load(f)


###############################################################################
# main

loss = []
init = []

for plotname in plots:

    plot = plots[plotname]

    L, U, N, T, x, plotfile = plot['attrs']
    dx = float(U - L) / N

    nsw, naw = plot['nsw'], plot['naw']

    sw, aw = abb_init_kernels(L, U, N, 'eviction', plotname)

    for i, t in enumerate(sw.years[:-1]):

        j0 = T.index(sw.years[i]) + 1
        j1 = T.index(sw.years[i+1])

        psw0 = sw.population(nsw[j0])
        psw1 = sw.population(nsw[j1])

        loss.append(abs(psw0 - psw1) / abs(psw0))


        j0 = T.index(sw.years[i])
        j1 = j0 + 1

        psw0 = len(sw.meas['SW'][T[j0]]) / sw.plot_size * 1e4
        psw1 = sw.population(nsw[j1])

        init.append(abs(psw0 - psw1) / abs(psw0))

print 'maximum empirical init error (percent):', max(init)*100.0
print 'maximum tree loss (percent):', max(loss)*100.0

