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
from utils.io import read_csv

from base import Kernel


class ABB(Kernel):


    def update(self, nsw, naw, t, **kwargs):
        self.nsw = nsw
        self.naw = naw

        # cache sw_pct, sw_psd, aw_pba etc...
        ix = np.array(range(len(self.x)))
        self._sw_pct  = self.sw_pct(self.x, ix)
        self._sw_psd  = self.sw_psd(self.x, ix)
        self._aw_pba  = self.aw_pba(self.x, ix)
        self._net_pba = self.net_pba(self.x, ix)

        # at this point, eg, _sw_pct[i] is the number of spruce with
        # dbh greater than x[i] where x corresponds to the methods'
        # grid

        self.method.sample(self, t)


    def sw_pct(self, x, ix):
        """Return number of Spruce with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i, j in enumerate(ix):
            r[i] = np.dot(self.method.P[j:], self.nsw[j:])

        return r / 1e3          # convert # / ha to # / 10 m^2


    def sw_psd(self, x, ix):
        """Return Spruce sum-of-diameters for Spruce with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i, j in enumerate(ix):
            r[i] = np.dot(self.method.P[j:], self.nsw[j:] * self.x[j:])

        return r / 1e3          # convert mm / ha to m / ha


    def aw_pba(self, x, ix):
        """Return Aspen basal area for Aspen with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i, j in enumerate(ix):
            r[i] = np.dot(self.method.P[j:], self.naw[j:] * pi * (self.x[j:]/2.0)**2)

        return r / 1e6          # convert mm^2 / ha to m^2 / ha


    def net_pba(self, x, ix):
        """Return total basal area for trees with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i, j in enumerate(ix):
            r[i] = np.dot(self.method.P[j:], (self.naw[j:] + self.nsw[j:]) * pi * (self.x[j:]/2.0)**2)

        return r / 1e6          # convert mm^2 / ha to m^2 / ha


    # XXX: units for all of the below?

    def sw_pct_from_meas(self, x):
        """Return number of Spruce with dbh greater than *x*."""

        r = [ d for d in self.sw0 if d > x ]

        return float(len(r)) / 1e3          # convert # / ha to # / 10 m^2


    def sw_psd_from_meas(self, x):
        """Return Spruce sum-of-diameters for Spruce with dbh greater than *x*."""

        r = [ d for d in self.sw0 if d > x ]

        return float(sum(r)) / 1e3          # convert mm / ha to m / ha


    def aw_pba_from_meas(self, x):
        """Return Aspen basal area for Aspen with dbh greater than *x*."""

        r = [ pi * (d/2.0)**2 for d in self.aw0 if d > x ]

        return sum(r) / 1e6          # convert mm^2 / ha to m^2 / ha


    def net_pba_from_meas(self, x):
        """Return total basal area for trees with dbh greater than *x*."""

        r = [ pi * (d/2.0)**2 for d in self.sw0 if d > x ]
        r.extend([ pi * (d/2.0)**2 for d in self.aw0 if d > x ])

        return sum(r) / 1e6          # convert mm^2 / ha to m^2 / ha



class ABBSW(ABB):


    def __init__(self):
        ABB.__init__(self)

        self.growth_params = np.array(
            [ 1.358768, 1.025693, -5.700839e-5, -2.256733e-1, -1.882815e-2 ])

        self.survival_params = np.array(
            [ 3.3668568, 0.023043, -0.183418, -0.041333 ])

        self.sd = sqrt(3.687244)

        measurements = read_csv('kernels/abb/sw.csv', header=['dbh'])
        self.n0 = np.asarray([ x.dbh for x in measurements ], dtype=np.float64)


    def kernel(self, x, y, t, ix=None, **kwargs):

        if ix is not None:
            one    = np.ones(len(x))
            dbh    = x
            sw_pct = self._sw_pct[ix]
            aw_pba = self._aw_pba[ix]
        else:
            one    = np.ones(len(y))
            dbh    = np.array(len(y) * [x])
            sw_pct = np.array(len(y) * [self.sw_pct_from_meas(x)])
            aw_pba = np.array(len(y) * [self.aw_pba_from_meas(x)])

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

        measurements = read_csv('kernels/abb/aw.csv', header=['dbh'])
        self.n0 = np.asarray([ x.dbh for x in measurements ], dtype=np.float64)


    def kernel(self, x, y, t, ix=None, **kwargs):

        if ix is not None:
            one     = np.ones(len(x))
            dbh     = x
            di      = y - x
            sw_psd  = self._sw_psd[ix]
            aw_pba  = self._aw_pba[ix]
            net_pba = self._net_pba[ix]
        else:
            one     = np.ones(len(y))
            dbh     = np.array(len(y) * [x])
            di      = y - x
            sw_psd  = np.array(len(y) * [self.sw_psd_from_meas(x)])
            aw_pba  = np.array(len(y) * [self.aw_pba_from_meas(x)])
            net_pba = np.array(len(y) * [self.net_pba_from_meas(x)])


        # growth
        xi = np.vstack([ one, dbh, aw_pba, sw_psd ])
        mu = np.dot(self.growth_params, xi)

        g = dnorm(y, mu, sd=self.sd)

        # survival
        xi = np.vstack([ one, dbh, dbh**2, di, dbh**2 * net_pba ])
        mu = np.dot(self.survival_params, xi)

        s  = exp(mu) / (1.0 + exp(mu))

        return g * s



