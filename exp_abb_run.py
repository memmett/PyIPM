"""Alberta boreal SW and AW competition experiment."""

import cPickle as pickle
import collections
import logging
import time

import numpy as np

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

import matplotlib.pylab as plt


###############################################################################
# config

N = 400
L = 10.0
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
# main

logging.info("start:  %s", time.asctime())
logging.info("method: %s", method.name)

n = np.zeros((len(T), 2, len(sw.x)))

# XXX: initial condition
M = 0.5 * (L + U)
n[0, 0] = np.exp(-(sw.x - M)**2/(M/4)**2)
n[0, 1] = 0.5 * np.exp(-(sw.x - M)**2/(M/4)**2)

# for k in [ sw, aw ]:
#     k.update(n[0, 0], n[0, 1], 0)
# print sw.sw_pct(np.array([0.0]))
# print sw.sw_pct(np.array([M]))
# raise SystemExit



plt.plot(n[0, 0], '-b')
plt.plot(n[0, 1], '-r')
plt.title(str(0))
plt.figure()
plt.show()

for j, t in enumerate(T[1:]):

    logging.info('step: %d (%s)', t, time.asctime())

    # set current populations and resample
    for k in [ sw, aw ]:
        k.update(n[j, 0], n[j, 1], t)

    # project
    n[j+1, 0] = sw.project(n[j, 0])
    n[j+1, 1] = aw.project(n[j, 1])


    plt.plot(n[j+1, 0], '-b')
    plt.plot(n[j+1, 1], '-r')
    plt.title(str(j+1))
    plt.figure()
