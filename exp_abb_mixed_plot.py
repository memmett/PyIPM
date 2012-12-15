"""Alberta boreal SW and AW competition experiment."""


import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

import cPickle as pickle

import matplotlib.pylab as plt


###############################################################################
# load results

with open('out/abb_mixed.pkl', 'r') as f:
    plots = pickle.load(f)


###############################################################################
# main

pens = {
    'C':  { 'color': 'b', 'marker': 's' },
    'D':  { 'color': 'r', 'marker': 'o' },
    'CD': { 'color': 'c', 'marker': '^' },
    'DC': { 'color': 'm', 'marker': 'v' },
}

for plotname in plots:

    plot = plots[plotname]

    # init

    sw = ABBSW()
    aw = ABBAW()

    L, U, N, T, x, name = plot['attrs']
    nsw, naw = plot['nsw'], plot['naw']

    dx = float(U - L) / N

    # setup kernels and connect them to each other
    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(MidPoint(), N)

    mrates = []
    for j in range(1, nsw.shape[0]-1):
        mr = sw.population(nsw[j+1]) / sw.population(nsw[j])
        mrates.append(mr)
    mrates = np.asarray(mrates)
    print name, np.mean(mrates), np.var(mrates)

    # plt.figure()
    # for j, t in enumerate(T):
    #     plt.plot(sw.x, nsw[j], '-b')
    #     plt.plot(aw.x, naw[j], '-r')

    psw = [ sw.population(nsw[j]) for j in range(len(T)) ]
    paw = [ aw.population(naw[j]) for j in range(len(T)) ]

    psw[0] = sum(nsw[0])
    paw[0] = sum(naw[0])

    plt.figure(1)
    plt.plot(T[1:], psw[1:], label=plotname, **pens[plotname])

    plt.figure(2)
    plt.plot(T[1:], paw[1:], label=plotname, **pens[plotname])

for f in [ 1, 2 ]:
    plt.figure(f)
    plt.legend(loc='best')
    plt.xlabel('year')
    plt.ylabel('population')
    #plt.title('population vs time ' + plotname)

plt.figure(1)
plt.savefig('plots/popmixed_sw.png')

plt.figure(2)
plt.savefig('plots/popmixed_aw.png')

plt.show()

