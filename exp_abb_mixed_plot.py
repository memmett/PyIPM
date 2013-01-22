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

figsph, sph = plt.subplots(1, 3, sharey=True)
figba,  ba  = plt.subplots(1, 3, sharey=True)


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

    psw = np.asarray([ sw.population(nsw[j]) for j in range(len(T)) ])
    paw = np.asarray([ aw.population(naw[j]) for j in range(len(T)) ])

    psw[0] = sum(nsw[0])
    paw[0] = sum(naw[0])

    sph[0].plot(T[1:], psw[1:], label=plotname, **pens[plotname])
    sph[1].plot(T[1:], paw[1:], label=plotname, **pens[plotname])
    sph[2].plot(T[1:], psw[1:] + paw[1:], label=plotname, **pens[plotname])

    basal_area = lambda k, s: np.pi * np.dot(k.method.P, s * (k.x/2)**2) / 1e6
    swba = np.asarray([ basal_area(sw, nsw[j]) for j in range(len(T)) ])
    awba = np.asarray([ basal_area(aw, naw[j]) for j in range(len(T)) ])

    ba[0].plot(T[1:], swba[1:], label=plotname, **pens[plotname])
    ba[1].plot(T[1:], awba[1:], label=plotname, **pens[plotname])
    ba[2].plot(T[1:], swba[1:] + awba[1:], label=plotname, **pens[plotname])


# sph
sph[2].legend(loc='best')

sph[0].set_ylabel('stems per hectare')
sph[1].set_xlabel('year')

sph[0].set_title('Spruce')
sph[1].set_title('Aspen')
sph[2].set_title('Total')

figsph.savefig('plots/sphmixed.png')


# ba
ba[0].legend(loc='best')

ba[0].set_ylabel('basal area')
ba[1].set_xlabel('year')

ba[0].set_title('Spruce')
ba[1].set_title('Aspen')
ba[2].set_title('Total')

figba.savefig('plots/bamixed.png')


plt.show()

