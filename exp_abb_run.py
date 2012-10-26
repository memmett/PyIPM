"""Alberta boreal SW and AW competition experiment."""

import logging
import time

import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW
from kernels.abb import first_year, last_year, measurements


###############################################################################
# config

N = 400
L = 10.0
U = 500.0

T = range(first_year, last_year+10)

###############################################################################
# init

method = MidPoint()

sw = ABBSW()
aw = ABBAW()

# setup kernels
for k in [ sw, aw ]:
    k.L, k.U = L, U
    k.setup(method, N)

# create root logger (debug messages to abb.log)
logging.basicConfig(filename='abb.log',
                    format='%(levelname)s: %(message)s',
                    level=logging.DEBUG)


# add stderr info logger
logger = logging.getLogger()
info   = logging.StreamHandler()
info.setLevel(logging.INFO)
logger.addHandler(info)


###############################################################################
# main

logging.info("start:  %s", time.asctime())
logging.info("method: %s", method.name)

n = np.zeros((len(T), 2, len(sw.x)))

for j, t in enumerate(T[1:]):

    logging.info('year: %d', t)

    # project
    if j == 0:
        for k in [ sw, aw ]:
            k.sw0 = sw.n0
            k.aw0 = aw.n0

        n[j, 0] = sw.method.histogram(sw.n0) / sw.plot_size * 1e4 # convert spp to sph
        n[j, 1] = aw.method.histogram(aw.n0) / aw.plot_size * 1e4 # convert spp to sph

        n[j+1, 0] = sw.first_projection() / sw.plot_size * 1e4    # convert spp to sph
        n[j+1, 1] = aw.first_projection() / aw.plot_size * 1e4    # convert spp to sph

    else:
        # set current populations and resample
        for k in [ sw, aw ]:
            k.update(n[j, 0], n[j, 1], t)

        n[j+1, 0] = sw.project(n[j, 0])
        n[j+1, 1] = aw.project(n[j, 1])


np.savez('abb.npz', x=sw.x, nsw=n[:, 0], naw=n[:, 1], L=L, U=U, N=N, T=T)
