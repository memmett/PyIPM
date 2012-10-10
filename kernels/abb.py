"""IPM kernels for the Alberta boreal SW and AW species."""

import logging
import numpy as np

from numpy import pi, exp, dot, ones, array
from utils.stats import dnorm

from base import Kernel


class ABB(Kernel):


    def __init__(self):
        Kernel.__init__(self)

        self.sw_pct  = np.frompyfunc(self._sw_pct,  1, 1)
        self.sw_psd  = np.frompyfunc(self._sw_psd,  1, 1)
        self.aw_pba  = np.frompyfunc(self._aw_pba,  1, 1)
        self.net_pba = np.frompyfunc(self._net_pba, 1, 1)

        
    def update(self, nsw, naw, t, **kwargs):
        self.nsw = nsw
        self.naw = naw
        self.method.sample(self, t)


    def _sw_pct(self, x):
        """Return number of Spruce with dbh greater than *x*."""

        idx = self.x >= x
        return dot(self.method.P[idx], self.nsw[idx])


    def _sw_psd(self, x):
        """Return Spruce sum-of-diameters for Spruce with dbh greater than *x*."""

        idx = self.x >= x
        return dot(self.method.P[idx], self.nsw[idx] * self.x[idx])



    def _aw_pba(self, x):
        """Return Aspen basal area for Aspen with dbh greater than *x*."""

        idx = self.x >= x
        return dot(self.method.P[idx], self.naw[idx] * pi * (self.x[idx]/2.0)**2)


    def _net_pba(self, x):
        """Return total basal area for trees with dbh greater than *x*."""

        idx = self.x >= x
        return dot(self.method.P[idx], (self.naw[idx] + self.nsw[idx]) * pi * (self.x[idx]/2.0)**2)



class ABBSW(ABB):


    def __init__(self):
        ABB.__init__(self)

        self.growth_params = array(
            [ 1.358768, 1.025693, -5.700839e-5, -2.256733e-1, -1.882815e-2 ])

        self.survival_params = array(
            [ 3.3668568, 0.023043, -0.183418, -0.041333 ])


    def kernel(self, x, y, t):

        one    = ones(len(x))
        dbh    = x
        sw_pct = self.sw_pct(x)
        aw_pba = self.aw_pba(x)

        # growth
        xi = array([ one, dbh, dbh**2, sw_pct, aw_pba ])
        mu = dot(self.growth_params, xi)

        g = dnorm(y, mu)

        # survival
        xi = array([ one, dbh, sw_pct, aw_pba ])
        mu = dot(self.survival_params, xi)

        s = exp(mu) / (1.0 + exp(mu))

        return g * s


class ABBAW(ABB):


    def __init__(self):
        ABB.__init__(self)

        self.growth_params = array(
            [ 3.6741201, 0.9954211, -0.0075639, -0.0139588 ])

        self.survival_params = array(
            [ 2.20, 0.207, -0.00286, -.659, -0.0000778 ])


    def kernel(self, x, y, t):

        one     = ones(len(x))
        dbh     = x
        di      = y - x
        sw_psd  = self.sw_psd(x)
        aw_pba  = self.aw_pba(x)
        net_pba = self.net_pba(x)

        # growth
        xi = array([ one, dbh, aw_pba, sw_psd ])
        mu = dot(self.growth_params, xi)

        g = dnorm(y, mu)

        # survival
        xi = array([ one, dbh**2, di, dbh**2 * net_pba ])
        mu = dot(self.survival_params, xi)

        s = exp(mu) / (1.0 + exp(mu))

        return g * s



