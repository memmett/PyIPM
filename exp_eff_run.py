"""Projection experiment."""

import cPickle as pickle
import collections
import logging
import time

import numpy as np

from time import clock

from methods import *
from kernels import *

from exp_eff_mesh import mesh_sizes

###############################################################################
# config

norm_order  = np.inf
growth_rate = True

kernels = [
    Exact,
    EED,
    Zuidema,
    ]

methods = [
    (MidPoint, {}),
    (MidPointZuidema,  {}),
    (INTClark, {}),
    (INTClark, { 'k': 3 }),
    (INTClark, { 'k': 5 }),
    (INTClark, { 'k': 9 }),
    (GENClark, {}),
    (GENClark, { 'k': 3 }),
    (GENClark, { 'k': 5 }),
    (GENClark, { 'k': 9 }),
    (GENClark, { 'k': 3, 'qtype': 'ClenshawCurtis' }),
    (GENClark, { 'k': 5, 'qtype': 'ClenshawCurtis' }),
    (GENClark, { 'k': 9, 'qtype': 'ClenshawCurtis' }),
    (GaussQuad, {}),
    (GaussQuad, { 'k': 3 }),
    (GaussQuad, { 'k': 5 }),
    (GaussQuad, { 'k': 9 }),
    (GaussQuad, { 'k': 3, 'adjust': True }),
    (GaussQuad, { 'k': 5, 'adjust': True }),
    (GaussQuad, { 'k': 9, 'adjust': True }),
    (GaussQuad, { 'qtype': 'ClenshawCurtis' }),
    (GaussQuad, { 'qtype': 'ClenshawCurtis', 'k': 3 }),
    (GaussQuad, { 'qtype': 'ClenshawCurtis', 'k': 5 }),
    (GaussQuad, { 'qtype': 'ClenshawCurtis', 'k': 9 }),
    ]

###############################################################################
# compute and plot projections

def compute_projections(kernel, method, mesh_size, growth_rate=False):
    """Compute projections and populations."""

    logging.info("PROJECT: %s, %s, %d", kernel.name, method.name, mesh_size)

    T = kernel.T
    kernel.reset_count()

    t1 = clock()
    kernel.setup(method, mesh_size)
    kernel.update(0)
    runtime = clock() - t1

    projections = np.zeros((len(T), len(kernel.x)))
    populations = np.zeros(len(T))

    logging.debug("PROJECT: first projection")

    n0 = kernel.first_projection()
    projections[0] = n0
    populations[0] = kernel.population(n0)

    for j, t in enumerate(T[1:]):

        n1 = kernel.project(n0)
        projections[j+1] = n1
        populations[j+1] = kernel.population(n1)

        n0 = n1

    if growth_rate:
        import scipy.sparse.linalg.eigen

        logging.debug("PROJECT: computing dominant eigenvalue")
        try:
            evals = scipy.sparse.linalg.eigen.eigs(
                kernel.projection_matrix, k=1, return_eigenvectors=False)
            growth_rate = evals[0]
        except:
            growth_rate = 0.0

    else:
        growth_rate = 0.0

    return projections, populations, growth_rate, runtime


###############################################################################
# create loggers

# create root logger (debug messages to project.log)
logging.basicConfig(filename='project.log',
                    format='%(levelname)s: %(message)s',
                    level=logging.DEBUG)

# add stderr info logger
logger = logging.getLogger()
info   = logging.StreamHandler()
info.setLevel(logging.INFO)
logger.addHandler(info)


###############################################################################
# main

logging.info("PROJECT: RUN STARTED AT: %s", time.asctime())

measurements = collections.defaultdict(dict)

for Kernel in kernels:

    kernel = Kernel()

    # compute a reference solution
    logging.info("REFERENCE: computing reference solution")

    method = GaussQuad('GaussLegendre', k=13)
    proj, pop, gr, time = compute_projections(kernel, method, 800, growth_rate)

    ref_pop = pop[-1]
    ref_gr  = gr

    # cycle through methods and project!
    for Method, kwargs in methods:

        method = Method(**kwargs)

        logging.info("PROJECT: method: %s", method.name)

        errors = []
        counts = []
        times  = []

        # do a throw-away run to prime it the method...
        kernel.reset_count()
        kernel.setup(method, 100)
        kernel.update(0)

        for N in mesh_sizes.get((kernel.name, method.name), []):

            proj, pop, gr, time = compute_projections(kernel, method, N, growth_rate)

            errors.append((abs(pop[-1] - ref_pop), abs(gr - ref_gr)))

            counts.append(kernel.count)
            times.append(time)

        if errors:
            measurements[kernel.name][method.name] = (errors, counts, times)


with open('measurements.pkl', 'w') as f:
    pickle.dump(measurements, f)
