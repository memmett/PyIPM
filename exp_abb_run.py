"""Alberta boreal SW and AW competition experiment."""


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
U = 500.0

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

# # XXX: initial condition
# M = 0.5 * (L + U)
# n[0, 0] = 5 * np.exp(-(sw.x - M)**2/(M/4)**2)
# n[0, 1] = 10 * np.exp(-(sw.x - M)**2/(M/4)**2)

# plt.figure()
# plt.plot(n[0, 0], '-b')
# plt.plot(n[0, 1], '-r')
# plt.title(str(0))

for j, t in enumerate(T[1:]):

    logging.info('step: %d (%s)', t, time.asctime())

    # project
    if j == 0:
        for k in [ sw, aw ]:
            k.sw0 = sw.n0
            k.aw0 = aw.n0

        n[j, 0] = sw.method.histogram(sw.n0) / sw.plot_size * 1e4 # convert spp to sph
        n[j, 1] = aw.method.histogram(aw.n0) / aw.plot_size * 1e4 # convert spp to sph

        n[j+1, 0] = sw.first_projection() / sw.plot_size * 1e4 # convert spp to sph
        n[j+1, 1] = aw.first_projection() / aw.plot_size * 1e4 # convert spp to sph

    else:
        # set current populations and resample
        for k in [ sw, aw ]:
            k.update(n[j, 0], n[j, 1], t)

        n[j+1, 0] = sw.project(n[j, 0])
        n[j+1, 1] = aw.project(n[j, 1])


    plt.figure()
    plt.plot(sw.x, n[j+1, 0], '-b', label='spruce')
    plt.plot(aw.x, n[j+1, 1], '-r', label='aspen')
    plt.xlabel('dbh (mm)')
    plt.ylabel('stems per hectare')
    plt.title('year ' + str(j+1))
    plt.savefig('plots/projection_%02d.pdf' % (j+1))



psw = [ sw.population(n[j, 0]) for j in range(len(T)) ]
paw = [ aw.population(n[j, 1]) for j in range(len(T)) ]

plt.figure()
plt.plot(T, psw, '-ob', label='spruce')
plt.plot(T, paw, '-or', label='aspen')
plt.legend(loc='best')
plt.xlabel('year')
plt.ylabel('population')
plt.title('population vs time')
plt.savefig('plots/population.pdf')

plt.show()
