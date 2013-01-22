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
    (Exact, {}),
    (EED, {}),
    (EED, { 'growth_type': 'slow' }),
    (Zuidema, {}),
]

methods = [
    (MidPoint, {}),
    (MidPointZuidema, {}),
    (INTClark, { 'k': 9 }),
    (GENClark, { 'k': 9 }),
    (GaussQuad, { 'k': 9 }),
    (GaussQuad, { 'k': 9, 'adjust': True }),
    (GaussQuad, { 'qtype': 'ClenshawCurtis', 'k': 9 }),

    # disabled below after receiving comments from reviewers
  # (INTClark, {}),
  # (INTClark, { 'k': 3 }),
  # (INTClark, { 'k': 5 }),
  # (GENClark, {}),
  # (GENClark, { 'k': 3 }),
  # (GENClark, { 'k': 5 }),
  # (GENClark, { 'k': 3, 'qtype': 'ClenshawCurtis' }),
  # (GENClark, { 'k': 5, 'qtype': 'ClenshawCurtis' }),
  # (GENClark, { 'k': 9, 'qtype': 'ClenshawCurtis' }),
  # (GaussQuad, {}),
  # (GaussQuad, { 'k': 3 }),
  # (GaussQuad, { 'k': 5 }),
  # (GaussQuad, { 'k': 3, 'adjust': True }),
  # (GaussQuad, { 'k': 5, 'adjust': True }),
  # (GaussQuad, { 'qtype': 'ClenshawCurtis' }),
  # (GaussQuad, { 'qtype': 'ClenshawCurtis', 'k': 3 }),
  # (GaussQuad, { 'qtype': 'ClenshawCurtis', 'k': 5 }),
]

###############################################################################
# compute and plot projections

def compute_projections(kernel, method, mesh_size,
                        growth_rate=False, plot_kernel=False):
    """Compute projections and populations."""

    logging.info("PROJECT: %s, %s, %d", kernel.name, method.name, mesh_size)

    T = kernel.T
    kernel.reset_count()

    t1 = clock()
    kernel.setup(method, mesh_size)
    kernel.update(0)
    runtime = clock() - t1

    if plot_kernel:
      import pylab
      pylab.imshow(kernel.method.A)
      pylab.show()

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
        from scipy.sparse.linalg.eigen import eigs, ArpackNoConvergence

        logging.debug("PROJECT: computing dominant eigenvalue")
        try:
            evals = eigs(kernel.projection_matrix,
                         k=1, tol=1e-14, maxiter=500, sigma=1.0, which='LM',
                         return_eigenvectors=False)
            gr = evals[0]

        except ArpackNoConvergence as err:
            logging.info("WARNING: EVALS: ARPACK didn't converge")
            gr = np.nan

        # import scipy.linalg
        # evals = scipy.linalg.eig(kernel.projection_matrix, left=False, right=False)
        # gr = abs(evals).max()

    else:
        gr = 0.0

    return projections, populations, gr, runtime


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

for Kernel, kernel_init_args in kernels:

    kernel = Kernel(**kernel_init_args)

    # compute a reference solution
    logging.info("REFERENCE: computing reference solution")

    method = GaussQuad('GaussLegendre', k=13)
    proj, pop, gr, time = compute_projections(kernel, method, 800, growth_rate)

    ref_pop = pop[-1]
    ref_gr  = gr

    # cycle through methods and project!
    for Method, method_init_args in methods:

        method = Method(**method_init_args)

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

            errors.append((pop[-1] - ref_pop, gr - ref_gr))

            counts.append(kernel.count)
            times.append(time)

        if errors:
            measurements[kernel.name][method.name] = (errors, counts, times)


with open('measurements.pkl', 'w') as f:
    pickle.dump(measurements, f)
