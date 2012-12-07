"""Alberta boreal SW and AW competition experiment.

This runs all flavours of kernels as defined in exp_abb_kernels.py.
"""

import logging
import time

import numpy as np
import cPickle as pickle

from glob import glob
from exp_abb_kernels import flavours, abb_init_kernels


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


###############################################################################
# main

logging.info("START: %s", time.asctime())

for flavour in flavours:

    logging.info("FLAVOUR: %s", flavour)

    plots = {}
    for plotfile in glob('kernels/abb/*.csv'):
        if plotfile == 'kernels/abb/plotsizes.csv':
            continue
        plotname = plotfile.split('/')[-1].split('.')[0]

        sw, aw = abb_init_kernels(L, U, N, flavour, plotname)

        logging.info('PLOT: %s %f', plotname, sw.plot_size)

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


        plots[plotname] = {}
        plots[plotname]['nsw'] = np.asarray(n[:, 0])
        plots[plotname]['naw'] = np.asarray(n[:, 1])
        plots[plotname]['attrs'] = (L, U, N, T, sw.x, plotfile)



    with open('out/abb_%s.pkl' % flavour, 'w') as f:
      pickle.dump(plots, f)
