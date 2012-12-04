"""Alberta boreal SW and AW competition experiment."""

import logging
import time

import numpy as np
import cPickle as pickle

from methods import MidPoint
from kernels.abb import ABBSW, ABBAW

from glob import glob


np.set_printoptions(linewidth=200)


###############################################################################
# config

N = 801
L = 0.0
U = 800.0


###############################################################################
# init

sw = ABBSW()
aw = ABBAW()

# setup kernels
sw.L, sw.U = L, U
sw.setup(MidPoint(), N)

aw.L, aw.U = L, U
aw.setup(MidPoint(), N)


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
logging.info("method: %s", sw.method.name)

plots = {}

for plotfile in glob('kernels/abb/*.csv'):
    if plotfile == 'kernels/abb/plotsizes.csv':
        continue

    plotname = plotfile.split('/')[-1].split('.')[0]
    logging.info('plot: %s', plotname)

    for k in [ sw, aw ]:
      k.measurements(plotfile, plotname)

    logging.info('plot size: %f', sw.plot_size)

    T = range(sw.first_year, sw.last_year+10)
    n = np.zeros((len(T)+1, 2, len(sw.x)))

    for j, t in enumerate(T):

        logging.debug('year: %d', t)

        if t in sw.years:
            from_data = True
            for k in [ sw, aw ]:
                k.set_n0(t)
        else:
            from_data = False


        # project
        if from_data:
            for k in [ sw, aw ]:
                k.sw0 = sw.n0
                k.aw0 = aw.n0

            n[j+1, 0] = sw.first_projection() / sw.plot_size * 1e4
            n[j+1, 1] = aw.first_projection() / aw.plot_size * 1e4

        else:
            # set current populations and resample
            for k in [ sw, aw ]:
                k.update(n[j, 0], n[j, 1], t)

            n[j+1, 0] = sw.project(n[j, 0])
            n[j+1, 1] = aw.project(n[j, 1])

            # dx = (U - L) / N

            # print 'A'
            # print sw.method.A[:8,:8]
            # # print 'k'
            # # print sw.kernel(sw.x[:5], sw.x[0], 0.0, ix=np.asarray(range(10))) * dx

            # import pylab as plt
            # plt.plot(sw.x, n[j, 0], '-b')
            # plt.plot(sw.x, n[j+1, 0], '-r')
            # plt.show()


    plots[plotname] = {}
    plots[plotname]['nsw'] = np.asarray(n[:, 0])
    plots[plotname]['naw'] = np.asarray(n[:, 1])
    plots[plotname]['attrs'] = (L, U, N, T, sw.x, plotfile)



with open('abb.pkl', 'w') as f:
  pickle.dump(plots, f)
