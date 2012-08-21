"""Clark's method.

There are several variants that can be constructed.  First, we
distinguish between kernels that can be integrated symbolically with
respect to y over a grid cell, and denote the method used to project
these by INTClark (integrable).  Otherwise, we denote the method by
GENClark.  Next, we postfix INTClark or GENClark with the numerical
integration scheme used to perform (any remaining) integrals:

* QP - QuadPack ('integrate' in R)

* GL(X) - Fixed point Gauss-Legendre with X nodes.

"""

import logging

import numpy as np
import base

from stats import pnorm
from math import sqrt

import scipy.integrate


###############################################################################

class INTClark(base.Method):
    """Clark method for kernels that can be symbolically integrate wrt y."""

    mesh_type = 'cell'

    def __init__(self, k=None):

        if k:
            name = 'INTClarkGL(' + str(k) + ')'
        else:
            name = 'INTClarkQP'

        self.k = k
        self.name = name
        # self.fixed_quad = fquad


    def sample(self, kernel, t):

        N = self.N
        U = self.U
        L = self.L
        X = self.X

        A = np.zeros((N, N))
        for i in range(N):
            logging.debug('%s: sampling (%d): %d', self.name, N, i)

            for j in range(N):

                F = kernel.F(X[j], X[j+1], t)

                def dx(x):
                    s  = kernel.s(x, t)
                    r  = kernel.r(x, t)
                    mu = kernel.mu(x, t)
                    sd = kernel.sd(x, t)

                    p2 = pnorm(X[j+1], mu, sd)
                    p1 = pnorm(X[j],   mu, sd)

                    return s * (p2 - p1) + r * F

                if self.k:
                    v, _ = scipy.integrate.fixed_quad(dx, X[i], X[i+1], n=self.k)
                else:
                    v, _ = scipy.integrate.quad(dx, X[i], X[i+1])

                dx = (U - L) / N
                A[j, i] = v / dx

        self.A = A 

        logging.debug("%s: update: kernel sampled", self.name)


###############################################################################

class GENClark(base.Method):
    """Clark method for general kernels."""

    mesh_type = 'cell'


    def __init__(self, k=None, qtype='GaussLegendre'):

        if k:
            if qtype == 'GaussLegendre':
                name = 'GENClarkGL(' + str(k) + ')'
                from scipy.special.orthogonal import p_roots
                self.xi, self.w = p_roots(k)
            elif qtype == 'ClenshawCurtis':
                name = 'GENClarkCC(' + str(k) + ')'
                import nodes
                self.xi = np.asarray(map(float, nodes.clenshaw_curtis_nodes(k-1)))
                self.w  = np.asarray(map(float, nodes.clenshaw_curtis_weights(k-1)))
        else:
            name = 'GENClarkQP'

        self.k = k
        self.name = name

    
    def sample(self, kernel, t):

        N = self.N
        U = self.U
        L = self.L
        X = self.X

        h = (U - L) / N
        s = 0.5 * h
        A = np.zeros((N, N))
        for i in range(N):
            logging.debug('%s: sampling (%d): %d', self.name, N, i)

            if self.k:
                xp = 0.5*(X[i+1] + X[i]) + 0.5*(X[i+1] - X[i]) * self.xi


            for j in range(N):

                if self.k:
                    yp = 0.5*(X[j+1] + X[j]) + 0.5*(X[j+1] - X[j]) * self.xi

                    v = 0.0
                    for m, x in enumerate(xp):
                        v += s**2 * self.w[m] * np.dot(self.w, kernel.kernel(x, yp, t))

                else:

                    dxdy = lambda x, y: kernel.kernel(x, y, t)
                    v, _ = scipy.integrate.dblquad(dxdy, X[j], X[j+1], 
                                                   lambda x: X[i], lambda x: X[i+1])
                A[j, i] = v / h

        self.A = A 

        logging.debug("%s: update: kernel sampled", self.name)
