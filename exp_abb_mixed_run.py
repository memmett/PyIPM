"""Alberta boreal SW and AW competition experiment."""

import logging
import time

import numpy as np
import cPickle as pickle

from kernels.abb import ABBSW, ABBAW
from methods.midpoint import MidPoint


###############################################################################
# config

N = 801
L = 0.0
U = 800.0


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

# some spruce, no aspen
sw0 = np.ones(sw.x.shape)
aw0 = np.zeros(sw.x.shape)
initial_populations.append((sw0, aw0, 'spruce-only'))

# some spruce, some aspen
sw0 = np.ones(sw.x.shape)
aw0 = np.ones(sw.x.shape)
initial_populations.append((sw0, aw0, 'mixed'))

# no spruce, some aspen
sw0 = np.zeros(sw.x.shape)
aw0 = np.ones(sw.x.shape)
initial_populations.append((sw0, aw0, 'aspen-only'))


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


with open('abb_mixed.pkl', 'w') as f:
  pickle.dump(plots, f)
