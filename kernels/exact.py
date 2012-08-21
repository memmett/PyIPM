"""Exact kernel (test kernel with a known solution)."""

from math import sqrt
from stats import dnorm, pnorm
from scipy.integrate import quad

import numpy as np

import base

class Exact(base.Kernel):

    def __init__(self):

        self.L = -10.0
        self.U =  30.0
        self.T = range(5)

        self.name = 'Exact'
        self.time_dependent = False

        self._ss = 0.1          # sigma^2
        self._mu = 2.05         # mean offset
        self._theta = 0.02      # survival

        self.mu0 = 0.4
        self.ss0 = 0.1


    def n0(self, x):

        return dnorm(x, mu=self.mu0, sd=sqrt(self.ss0))

    
    #### general form

    def kernel(self, x, y, t):

        self.increment_count(x, y)

        return np.exp(-self._theta * x**2) * dnorm(y, mu=x+self._mu, sd=sqrt(self._ss))

    #### y-integrable form

    def s(self, x, t):

        return np.exp(-self._theta * x**2)


    def mu(self, x, t):

        self.increment_count(x, [])

        return x + self._mu


    def sd(self, x, t):

        return sqrt(self._ss)


    def r(self, x, t):

        return 0.0


    def F(self, y1, y2, t):

        return 0.0





    # def exact(self, t, method):

    #     pdf = lambda x: self.theta**t * dnorm(
    #         x, mu=self.mu0 + t*self.mu, sd=sqrt(self.ss0+t*self.ss))

    #     cdf = lambda x: self.theta**t * pnorm(
    #         x, mu=self.mu0 + t*self.mu, sd=sqrt(self.ss0+t*self.ss))

    #     n = np.zeros(method.N)
    #     if method.mesh_type == 'cell':
    #         for i in range(method.N):
    #             n[i] = cdf(method.X[i+1]) - cdf(method.X[i])
    #     else:
    #         n = np.zeros(method.N)
    #         for i in range(method.N):
    #             n[i] = pdf(method.x[i])

    #     return n

