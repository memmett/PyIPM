"""ARTR experiment."""

import logging
import time
from itertools import product
from collections import defaultdict

import numpy as np
import matplotlib.pylab as plt
import scipy.stats

from scipy.stats.mstats import mquantiles

from kernels.artr import ARTTRI
from methods.midpoint import MidPoint


###############################################################################
# config

iterations = 1

kernel = ARTTRI()
method = MidPoint()
N      = 400

mortality_types = [ 'noexp' ]            # see kernels/artr.py 
fecundity_types = [ 'uniform_noexp' ]    # see kernels/artr.py

adult_cutoff = 0.0

covariats = [ 'pptLag1', 'pptWin2', 'TmeanSum1', 'recpptApr', 'recFebMarSnow' ]

elasticity_from_data = True


###############################################################################
# create loggers

# create root logger (debug messages to project.log)
logging.basicConfig(filename='artr.log',
                    format='%(levelname)s: %(message)s',
                    level=logging.DEBUG)

# add stderr info logger
logger = logging.getLogger()
info = logging.StreamHandler()
info.setLevel(logging.INFO)
logger.addHandler(info)


###############################################################################
# compute and plot projections

logging.info("ARTR: RUN STARTED AT: %s", time.asctime())
logging.info("ARTR: kernel: %s", kernel.name)
logging.info("ARTR: method: %s", method.name)
logging.info("ARTR: mesh size: %d", N)

#### compute projections

method = MidPoint()

kernel.reset_count()
kernel.setup(method, N)

juvenile_mask = method.x <  adult_cutoff
adult_mask    = method.x >= adult_cutoff

T = kernel.T

for mortality, fecundity in product(mortality_types, fecundity_types):

    kernel.mortality_type = mortality
    kernel.fecundity_type = fecundity

    elasticity  = defaultdict(list)
    projections = np.zeros((iterations, len(T)+1, N))
    populations = np.zeros((iterations, len(T)+1))
    juveniles   = np.zeros((iterations, len(T)+1))
    adults      = np.zeros((iterations, len(T)+1))

    for k in range(iterations):

        logging.info("ARTR: iteration: %d", k)

        n0 = method.histogram(kernel.n0)

        projections[k, 0] = n0
        populations[k, 0] = sum(n0)
        juveniles[k, 0]   = sum(n0 * juvenile_mask)
        adults[k, 0]      = sum(n0 * adult_mask)

        logging.debug("ARTR: iter: %d, step: 0, population: %f",
                      k, populations[k, 0])

        for j, t in enumerate(T):

            logging.info("ARTR: time: %d (%d)", t, j)

            kernel.update(t)
            if iterations > 1:
                kernel.draw_parameters(t)

            if j == 0:
                n1 = kernel.first_projection()
            else:
                n1 = kernel.project(n0)

            if fecundity == 'population':
                # kernel.fecundity_type = 'uniform'
                # y = kernel.x
                # f = np.zeros((N, N))
                # for kk, x in enumerate(kernel.x):
                #     f[kk, :] = kernel.fecundity(x, y, t)
                # n1 += np.dot(f, n0)
                # kernel.fecundity_type = 'population'
                
                n1 += kernel.population_fecundity(n0)

            if iterations == 1:
                for covariat in covariats:

                    if elasticity_from_data:
                        c0 = method.histogram(kernel.measurements(t))
                        if not np.any(c0):
                            elasticity[covariat].append(np.nan)
                            continue
                    else:
                        c0 = n0

                    c1 = kernel.project(c0)

                    kernel.tweak(covariat, 1.1)
                    kernel.update(t)
                    c2 = kernel.project(c0)

                    gr1 = method.growth_rate(c0, c1)
                    gr2 = method.growth_rate(c0, c2)

                    logging.debug("ELASTICITY: step: %d, "
                                  + "covariat: %s, lambda relative: %f", 
                                  j, covariat, gr2/gr1)
                    elasticity[covariat].append(gr2/gr1)

                kernel.tweak(None, 1.0)

            projections[k, j+1] = n1
            populations[k, j+1] = method.total_population(n1)
            juveniles[k, j+1]   = method.total_population(n1 * juvenile_mask)
            adults[k, j+1]      = method.total_population(n1 * adult_mask)

            logging.debug("ARTR: iter: %d, step: %d, population: %f",
                          k, j+1, populations[k, j+1])

            n0 = n1


    #### plots

    ## plot projections with predictive intervals
    # logging.info('PLOTTING: projections')

    # for j, t in enumerate(T):
    #     plt.figure()
    #     plt.title('%s: %s (%d), %d' % (kernel.name, method.name, N, t))

    #     q = mquantiles(projections[:, j, :], prob=[0.05, 0.5, 0.95], axis=0)
    #     plt.plot(method.x, q[0], '-b')
    #     plt.plot(method.x, q[1], '-r')
    #     plt.plot(method.x, q[2], '-b')
    #     # try:
    #     #     plt.hist(kernel.measurements(t), method.X, facecolor='none')
    #     #     plt.legend([ '0.05', '0.50', '0.95', 'meas' ])
    #     # except KeyError:
    #     #     plt.legend([ '0.05', '0.50', '0.95' ])
    #     plt.legend([ '0.05', '0.50', '0.95' ])

    #     plt.ylim(ymax=6.0)
    #     plt.savefig('projection-%d-%s-%s.pdf' % (t, mortality, fecundity))

    ## plot population vs time
    logging.info('PLOTTING: population')

    plt.figure()
    plt.title('Population: %s: %s (%d); %s' % (kernel.name, method.name, N, fecundity))

    q = mquantiles(populations, prob=[0.05, 0.5, 0.95], axis=0)
    t = list(T) + [ T[-1]+1 ]

    # quantiles
    plt.plot(t, q[0], '-b')
    plt.plot(t, q[1], '-r')
    plt.plot(t, q[2], '-b')

    # measurements
    plt.plot(T, [ len(kernel.measurements(t)) for t in T ], 's-k')

    plt.legend([ '0.05', '0.50', '0.95', 'meas' ])
    plt.savefig('population-%s-%s.pdf' % (mortality, fecundity))


    ## plot juvenile and adult population vs time
    logging.info('PLOTTING: juvenile/adult')

    plt.figure()
    plt.title('Juvenile/adult population: %s: %s (%d); %s' % 
              (kernel.name, method.name, N, fecundity))

    q = mquantiles(juveniles, prob=[0.05, 0.5, 0.95], axis=0)
    t = list(T) + [ T[-1]+1 ]

    plt.plot(t, q[0], '-m')
    plt.plot(t, q[1], '-g')
    plt.plot(t, q[2], '-m')


    q = mquantiles(adults, prob=[0.05, 0.5, 0.95], axis=0)
    t = list(T) + [ T[-1]+1 ]

    plt.plot(t, q[0], '-b')
    plt.plot(t, q[1], '-r')
    plt.plot(t, q[2], '-b')

    # measurements
    juv  = []
    adlt = []
    for t in T:
        meas = kernel.measurements(t)
        juv.append(len([ x for x in meas if x < adult_cutoff ]))
        adlt.append(len([ x for x in meas if x >= adult_cutoff ]))

    plt.plot(T, juv, 's-c')
    plt.plot(T, adlt, 's-k')

    plt.legend([ 'J 0.05', 'J 0.50', 'J 0.95', 
                 'A 0.05', 'A 0.50', 'A 0.95', 
                 'J meas', 'A meas' ])
    plt.savefig('juvadlt-%s-%s.pdf' % (mortality, fecundity))


    ## annual growth rate
    logging.info('PLOTTING: growth rates')

    plt.figure()
    plt.title('Annual growth rate: %s: %s (%d); %s'
              % (kernel.name, method.name, N, fecundity))

    growth_rates = np.zeros((iterations, len(T)))
    for j in range(len(T)):
        growth_rates[:, j] = populations[:, j+1] / populations[:, j]

    q = mquantiles(growth_rates, prob=[0.05, 0.5, 0.95], axis=0)
    plt.plot(T, q[0], '-b')
    plt.plot(T, q[1], '-r')
    plt.plot(T, q[2], '-b')

    q = mquantiles(populations[:, -1] / populations[:, 0], 
                   prob=[0.05, 0.5, 0.95])
    q = q**(1/float(len(T)))

    # print 'Long term growth rate:', q
    # print 'Eigen values:', sorted(np.linalg.eigvals(method.A), reverse=True)[:2]

    measured_growth_rates = []
    for t in T[:-1]:
        try:
            gr = ( float(len(kernel.measurements(t+1))) 
                   / len(kernel.measurements(t)) )
        except ZeroDivisionError:
            gr = 0.0
        measured_growth_rates.append(gr)
    plt.plot(T[:-1], measured_growth_rates, 's-k')

    plt.legend([ '0.05', '0.50', '0.95', 'meas' ])
    plt.savefig('growth_rate-%s-%s.pdf' % (mortality, fecundity))

    ## elasticity
    if iterations == 1:
        logging.info('PLOTTING: elasticity')

        plt.figure()
        for covariat in covariats:
            plt.plot(T, elasticity[covariat], label=covariat)
        plt.legend(loc='best')
        plt.xlabel('year')
        plt.ylabel('ratio of growth rates')
        plt.title('Change in growth rate due to 10% change in covariat')


#### done

# save projections an r data file
R = False
try:
    import rpy2.robjects as ro
    import rpy2.robjects.numpy2ri 
    R = True
except:
    pass

if R:
    rpy2.robjects.numpy2ri.activate()

    ro.r.assign('projections', projections)
    ro.r.assign('x', method.x)
    ro.r.save('projections', 'x', file='projections.rdata')
    logging.info('OUTPUT: projections saved to projections.rdata')

    for c in covariats:
        ro.r.assign('elasticity_%s' % c, np.asarray(elasticity[c]))
    ro.r.assign('t', np.asarray(T))
    save = [ 't' ] + [ 'elasticity_%s' % c for c in covariats ]
    ro.r.save(*save, file='elasticity.rdata')
    logging.info('OUTPUT: elasticity saved to elasticity.rdata')


plt.show()
