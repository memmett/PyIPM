"""Alberta boreal SW and AW mixed stand experiment."""

import logging
import time

import numpy as np
import cPickle as pickle

from kernels.abb import ABBSW, ABBAW
from methods.midpoint import MidPoint

from utils.stats import dttnorm

###############################################################################
# config

N = 1601
L = 0.0
U = 400.0


###############################################################################
# init

# create root logger (debug messages to abb.log)
logging.basicConfig(filename='abb.log',
                    format='%(levelname)s: %(message)s',
                    level=logging.DEBUG)


# add stderr info logger
logger = logging.getLogger()
info   = logging.StreamHandler()
info.setLevel(logging.INFO)
logger.addHandler(info)

logging.info("START: %s", time.asctime())


###############################################################################
# initial conditions

initial_populations = []

# instantiate an SW kernel to grab an 'x' array...
sw = ABBSW()
sw.L, sw.U = L, U
sw.setup(MidPoint(), N)


# dens: # / 10 m^2
scenarios = {
    'C': {
        'sw': { 'dens': 8.846, 'mu': 25.7, 'sd': 17.1, 'min':  3.8, 'max': 53.5 },
        'aw': { 'dens': 3.308, 'mu': 30.5, 'sd': 22.9, 'min':  9.4, 'max': 76.1 },
    },
    'CD': {
        'sw': { 'dens': 3.778, 'mu': 25.7, 'sd': 17.1, 'min':  3.8, 'max':  53.5 },
        'aw': { 'dens': 7.278, 'mu': 63.8, 'sd': 28.5, 'min': 23.7, 'max': 109.2 },
    },
    'D': {
        'sw': { 'dens': 0.036, 'mu': 25.7, 'sd': 17.1, 'min':  3.8, 'max':  53.5 },
        'aw': { 'dens': 6.964, 'mu': 52.5, 'sd': 27.3, 'min': 11.1, 'max': 112.5 },
    },
    'DC': {
        'sw': { 'dens': 2.529, 'mu': 25.7, 'sd': 17.1, 'min':  3.8, 'max': 53.5 },
        'aw': { 'dens': 9.118, 'mu': 29.3, 'sd': 19.4, 'min':  3.3, 'max': 68.6 },
    },
}

for scenario in scenarios:
    swp = scenarios[scenario]['sw']
    awp = scenarios[scenario]['aw']

    sw0 = swp['dens'] * dttnorm(sw.x, **swp) * 1000
    aw0 = awp['dens'] * dttnorm(sw.x, **awp) * 1000

    # import pylab as plt
    # plt.figure()
    # plt.plot(sw.x, sw0, '-b')
    # plt.plot(sw.x, aw0, '-r')
    # plt.title(scenario)

    initial_populations.append((sw0, aw0, scenario))

# plt.show()


###############################################################################
# main

plots = {}

for sw0, aw0, mix in initial_populations:

    logging.info("MIX: %s", mix)

    sw = ABBSW()
    aw = ABBAW()

    for k in [ sw, aw ]:
        k.L, k.U = L, U
        k.setup(MidPoint(), N)

    T = range(30)
    n = np.zeros((len(T)+1, 2, len(sw.x)))

    n[0, 0] = sw0
    n[0, 1] = aw0

    for j, t in enumerate(T):

        for k in [ sw, aw ]:
            k.update(n[j, 0], n[j, 1], t)

        n[j+1, 0] = sw.project(n[j, 0])
        n[j+1, 1] = aw.project(n[j, 1])


    plots[mix] = {}
    plots[mix]['nsw'] = np.asarray(n[:, 0])
    plots[mix]['naw'] = np.asarray(n[:, 1])
    plots[mix]['attrs'] = (L, U, N, T, sw.x, mix)


with open('out/abb_mixed.pkl', 'w') as f:
  pickle.dump(plots, f)
