"""Easterling, Ellner and Dixon IPM kernel."""

from stats import dnorm
from numpy import sqrt, exp

import base


class EED(base.Kernel):

    def __init__(self):

        base.Kernel.__init__(self)

        self.L =  0.0
        self.U = 10.0
        self.T = range(0, 5)

        self.name = 'EED'
        self.time_independent = True

        self.survival_params  = [ 1.34, 0.92 ]
        self.growth_params    = [ 0.37, 0.73, 0.127, 0.23 ]
        self.fecundity_params = [ 0.034, 0.038 ]

        self.discontinuities  = [ 0.15, 0.25 ]
 

    def n0(self, x):

        return dnorm(x, mu=5.0, sd=sqrt(1.0))


    #### general form

    def survival(self, x, t):

        a0, a1 = self.survival_params

        return 1.0/(1.0 + exp(-(a0 + a1 * x)))


    def growth(self, x, y, t):

        a0, a1, b0, b1 = self.growth_params

        m = a0 + a1 * x
        v = b0 + b1 * x

        return dnorm(y, mu=m, sd=sqrt(v))


    def fecundity(self, x, y, t):

        a0, a1 = self.fecundity_params

        f1 = a0 + a1 * x
        f2 = 1.0 * (0.15 <= y) * (y <= 0.25) / (0.25 - 0.15)

        return f1 * f2


    def kernel(self, x, y, t):

        self.increment_count(x, y)

        s = self.survival(x, t)
        g = self.growth(x, y, t)
        f = self.fecundity(x, y, t)

        return s*g + f


    #### y-integrable form

    def s(self, x, t):

        return self.survival(x, t)


    def mu(self, x, t):

        self.increment_count(x, [])

        a0, a1, b0, b1 = self.growth_params

        return a0 + a1 * x


    def sd(self, x, t):

        a0, a1, b0, b1 = self.growth_params

        return sqrt(b0 + b1 * x)


    def r(self, x, t):

        a0, a1 = self.fecundity_params
        return a0 + a1 * x


    def F(self, y1, y2, t):

        if y2 < 0.15 or 0.25 < y1:
            return 0.0

        if y1 < 0.15:
            y1 = 0.15

        if y2 > 0.25:
            y2 = 0.25
        
        return (y2 - y1) / (0.25 - 0.15)
