"""Exact kernel (test kernel with a known solution)."""

from numpy import sqrt, exp
from scipy.integrate import quad
from utils.stats import dnorm, pnorm
from base import Kernel


class Exact(Kernel):

    def __init__(self, **kwargs):

        self.L = 0.0
        self.U = 20.0
        self.T = range(5)

        self.name = 'Exact'
        self.time_dependent = False

        self._ss = 2.0
        self._mu = 1.0
        self._theta = 0.02

        self.mu0 = 5.0
        self.ss0 = 1.0


    def n0(self, x):

        return dnorm(x, mu=self.mu0, sd=sqrt(self.ss0))


    #### general form

    def kernel(self, x, y, t, **kwargs):

        self.increment_count(x, y)

        return exp(-self._theta * x**2) * dnorm(y, mu=x+self._mu, sd=sqrt(self._ss))


    #### y-integrable form

    def s(self, x, t):

        return exp(-self._theta * x**2)


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

