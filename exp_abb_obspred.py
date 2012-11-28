"""Alberta boreal SW and AW competition experiment."""


import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

import cPickle as pickle

import pylab as plt


###############################################################################
# load results

with open('abb.pkl', 'r') as f:
    plots = pickle.load(f)


###############################################################################
# main

pops = {
    'obs': { 'sw': [], 'aw': [] },
    'prd': { 'sw': [], 'aw': [] },
}

for plotname in plots:

    plot = plots[plotname]

    L, U, N, T, x, plotfile = plot['attrs']
    nsw, naw = plot['nsw'], plot['naw']

    method = MidPoint()
    sw = ABBSW()
    aw = ABBAW()

    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(method, N)
        k.measurements(plotfile, plotname)

    for obs_year in sw.years[1:]:
        pops['obs']['sw'].append(len(sw.meas['SW'][obs_year]) / sw.plot_size * 1e4)
        pops['obs']['aw'].append(len(aw.meas['AW'][obs_year]) / sw.plot_size * 1e4)

        j = T.index(obs_year)
        pops['prd']['sw'].append(sw.population(nsw[j]))
        pops['prd']['aw'].append(aw.population(naw[j]))


for species in [ 'sw', 'aw' ]:
    plt.figure()
    plt.plot(pops['obs'][species], pops['prd'][species], 'ok')
    plt.xlabel('observed')
    plt.ylabel('predicted')
    plt.title({'sw': 'Spruce', 'aw': 'Aspen'}[species])

plt.show()