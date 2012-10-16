"""IPM kernels for the Alberta boreal SW and AW species.

Note that throughout these kernels the units are:

  * nsw, naw:    # / ha
  * dbh (x):     mm
  * sw_pct:      # / 10m^2
  * sw_psd:      m / ha
  * aw_pba:      m^2 / ha
  * net_pba:     m^2 / ha
  * di:          cm / y
  *

"""

import logging
import numpy as np

from numpy import pi, exp, sqrt
from utils.stats import dnorm

from base import Kernel


class ABB(Kernel):


    def update(self, nsw, naw, t, **kwargs):
        self.nsw = nsw
        self.naw = naw
        self.method.sample(self, t)


    def sw_pct(self, x):
        """Return number of Spruce with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i in xrange(len(x)):
            idx = self.x >= x[i]
            r[i] = np.dot(self.method.P[idx], self.nsw[idx]) 

        return r / 1e3          # convert # / ha to # / 10 m^2


    def sw_psd(self, x):
        """Return Spruce sum-of-diameters for Spruce with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i in xrange(len(x)):
            idx = self.x >= x[i]
            r[i] = np.dot(self.method.P[idx], self.nsw[idx] * self.x[idx])

        return r / 1e3          # convert mm / ha to m / ha


    def aw_pba(self, x):
        """Return Aspen basal area for Aspen with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i in xrange(len(x)):
            idx = self.x >= x[i]
            r[i] = np.dot(self.method.P[idx], self.naw[idx] * pi * (self.x[idx]/2.0)**2)

        return r / 1e6          # convert mm^2 / ha to m^2 / ha


    def net_pba(self, x):
        """Return total basal area for trees with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i in xrange(len(x)):
            idx = self.x >= x[i]
            r[i] = np.dot(self.method.P[idx], (self.naw[idx] + self.nsw[idx]) * pi * (self.x[idx]/2.0)**2)

        return r / 1e6          # convert mm^2 / ha to m^2 / ha



class ABBSW(ABB):


    def __init__(self):
        ABB.__init__(self)

        self.growth_params = np.array(
            [ 1.358768, 1.025693, -5.700839e-5, -2.256733e-1, -1.882815e-2 ])

        self.survival_params = np.array(
            [ 3.3668568, 0.023043, -0.183418, -0.041333 ])

        self.sd = sqrt(3.687244)



    def kernel(self, x, y, t):

        one    = np.ones(len(x))
        dbh    = x
        sw_pct = self.sw_pct(x)
        aw_pba = self.aw_pba(x)

        # growth
        xi = np.vstack([ one, dbh, dbh**2, sw_pct, aw_pba ])
        mu = np.dot(self.growth_params, xi)

        g = dnorm(y, mu, sd=self.sd)

        # survival
        xi = np.vstack([ one, dbh, sw_pct, aw_pba ])
        mu = np.dot(self.survival_params, xi)

        s = exp(mu) / (1.0 + exp(mu))

        return g * s


class ABBAW(ABB):


    def __init__(self):
        ABB.__init__(self)

        self.growth_params = np.array(
            [ 3.6741201, 0.9954211, -0.0075639, -0.0139588 ])

        self.survival_params = np.array(
            [ 2.20, 0.0207, -0.0000286, -0.0659, -0.000000778 ])

        self.sd = sqrt(3.687244)


    def kernel(self, x, y, t):

        one     = np.ones(len(x))
        dbh     = x
        di      = y - x
        sw_psd  = self.sw_psd(x)
        aw_pba  = self.aw_pba(x)
        net_pba = self.net_pba(x)

        # growth
        xi = np.vstack([ one, dbh, aw_pba, sw_psd ])
        mu = np.dot(self.growth_params, xi)

        g = dnorm(y, mu, sd=self.sd)

        # survival
        xi = np.vstack([ one, dbh, dbh**2, di, dbh**2 * net_pba ])
        mu = np.dot(self.survival_params, xi)

        s  = exp(mu) / (1.0 + exp(mu))

        return g * s



