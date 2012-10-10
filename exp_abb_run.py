"""Alberta boreal SW and AW competition experiment."""

import cPickle as pickle
import collections
import logging
import time

import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW


###############################################################################
# config

N = 400
L = 0.0
U = 300.0

T = range(0, 10)

###############################################################################
# init

method = MidPoint()

sw = ABBSW()
aw = ABBAW()

# setup kernels and connect them to each other
for k in [ sw, aw ]:
    k.L, k.U = L, U
    k.setup(method, N)

    k.kernel_sw = sw
    k.kernel_aw = aw    

# create root logger (debug messages to project.log)
logging.basicConfig(filename='abb.log',
                    format='%(levelname)s: %(message)s',
                    level=logging.DEBUG)

# add stderr info logger
logger = logging.getLogger()
info   = logging.StreamHandler()
info.setLevel(logging.INFO)
logger.addHandler(info)


###############################################################################
# test

# n = np.zeros((2, len(sw.x)))
# n[0] = np.exp(-(sw.x - U/2)**2/(U/8)**2)
# n[1] = 0.5 * np.exp(-(sw.x - U/2)**2/(U/8)**2)

# import matplotlib.pylab as plt
# plt.plot(sw.x, n[0], '-r')
# plt.plot(sw.x, n[1], '-b')
# plt.show()


# set current populations and resample
# for k in [ sw, aw ]:
#     k.update(n[0], n[1], 0)

# print sw.sw_pct(np.array([100.0, 101.0]))


###############################################################################
# main

logging.info("RUN STARTED AT: %s", time.asctime())
logging.info("method: %s", method.name)

n = np.zeros((len(T), 2, len(sw.x)))

# XXX: initial condition

for j, t in enumerate(T[1:]):

    # set current populations and resample
    for k in [ sw, aw ]:
        k.update(n[j, 0], n[j, 1], t)

    # project
    n[j+1, 0] = sw.project(n[j, 0])
    n[j+1, 1] = aw.project(n[j, 1])

