"""IPM kernels for the Alberta boreal SW and AW species.

Note that throughout these kernels the units are:

  * nsw, naw:    # / ha
  * dbh (x):     mm
  * sw_pct:      # / 10m^2
  * sw_psd:      m / ha
  * aw_pba:      m^2 / ha
  * net_pba:     m^2 / ha
  * di:          cm / y

If you change the initial conditions (in abb/sw.csv and abb/aw.csv)
then you should also update the plot sizes.

The plot size is assumed to be in units of m^2.

"""

import logging
import numpy as np

from collections import defaultdict
from numpy import pi, exp, sqrt
from utils.stats import dnorm, dtnorm
from utils.io import read_csv

from base import Kernel

plotsizes = read_csv('kernels/abb/plotsizes.csv')


class ABB(Kernel):

    def __init__(self):

        Kernel.__init__(self)

        self.sw_mort     = 'model'
        self.competition = True
        self.no_mort     = False


    def measurements(self, raw_meas, plotname):

        raw_meas = read_csv(raw_meas)
        plotsize = [ x.plotsize for x in plotsizes if x.plot == plotname ]
        if len(plotsize) == 1:
            self.plot_size = float(plotsize[0])
        else:
            raise ValueError('something funny with the plotsize')

        self.meas = {
            'SW': defaultdict(list),
            'AW': defaultdict(list),
        }

        for r in raw_meas:
            if r.spec in [ 'SW', 'AW' ]:
                self.meas[r.spec][int(r.year)].append(float(r.dbh))

        self.years      = sorted(list(set([ int(x.year) for x in raw_meas ])))
        self.first_year = self.years[0]
        self.last_year  = self.years[-1]


    def update(self, nsw, naw, t, **kwargs):
        self.nsw = nsw
        self.naw = naw

        # cache sw_psd, aw_pba etc...
        ix = np.array(range(len(self.x)))
        self._sw_psd  = self.sw_psd(self.x, ix)
        self._sw_pct  = self.sw_pct(self.x, ix)
        self._aw_pct  = self.aw_pct(self.x, ix)
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


    def aw_pct(self, x, ix):
        """Return number of Spruce with dbh greater than *x*."""

        r = np.empty(x.shape)

        for i, j in enumerate(ix):
            r[i] = np.dot(self.method.P[j:], self.naw[j:])

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


    def sw_pct_from_meas(self, x):
        """Return number of Spruce with dbh greater than *x*."""

        r = [ d for d in self.sw0 if d > x ]

        return float(len(r)) / (self.plot_size / 10.0)  # convert # / plot_size to # / 10 m^2


    def aw_pct_from_meas(self, x):
        """Return number of aspen with dbh greater than *x*."""

        r = [ d for d in self.aw0 if d > x ]

        return float(len(r)) / (self.plot_size / 10.0)  # convert # / plot_size to # / 10 m^2


    def sw_psd_from_meas(self, x):
        """Return Spruce sum-of-diameters for Spruce with dbh greater than *x*."""

        r = [ d for d in self.sw0 if d > x ]

        return float(sum(r)) / (self.plot_size / 10.0)  # convert mm / plot_size to m / ha


    def aw_pba_from_meas(self, x):
        """Return Aspen basal area for Aspen with dbh greater than *x*."""

        r = [ pi * (d/2.0)**2 for d in self.aw0 if d > x ]

        return sum(r) / (self.plot_size * 100.0)      # convert mm^2 / plot_size to m^2 / ha


    def net_pba_from_meas(self, x):
        """Return total basal area for trees with dbh greater than *x*."""

        r = [ pi * (d/2.0)**2 for d in self.sw0 if d > x ]
        r.extend([ pi * (d/2.0)**2 for d in self.aw0 if d > x ])

        return sum(r) / (self.plot_size * 100.0)      # convert mm^2 / plot_size to m^2 / ha



class ABBSW(ABB):


    def __init__(self):
        ABB.__init__(self)

        self.growth_params = np.array(
            [ 1.430, 1.023, -3.883e-3, -6.659e-2, -5.435e-5 ])

        self.growth_params_no_comp = np.array(
            [ 0.5905, 1.029, 0.0, 0.0, -0.00006619 ])

        self.survival_params = np.array(
            [ 4.780, 0.0342, -0.00738, -0.0299 ])

        self.sd = sqrt(3.550338)
        self.sd_no_comp = sqrt(1.003784)


    def set_n0(self, year):

        self.n0 = np.asarray(self.meas['SW'][year])


    def kernel(self, x, y, t, ix=None, **kwargs):

        if ix is not None:
            one    = np.ones(len(x))
            dbh    = x
            sw_psd = self._sw_psd[ix]
            sw_pct = self._sw_pct[ix]
            aw_pct = self._aw_pct[ix]
            aw_pba = self._aw_pba[ix]
        else:
            one    = np.ones(len(y))
            dbh    = np.array(len(y) * [x])
            sw_psd = np.array(len(y) * [self.sw_psd_from_meas(x)])
            sw_pct = np.array(len(y) * [self.sw_pct_from_meas(x)])
            aw_pct = np.array(len(y) * [self.aw_pct_from_meas(x)])
            aw_pba = np.array(len(y) * [self.aw_pba_from_meas(x)])

        # growth
        if self.competition:
            params = self.growth_params
            sd = self.sd
        else:
            params = self.growth_params_no_comp
            sd = self.sd_no_comp

        xi = np.vstack([ one, dbh, sw_psd, aw_pct + sw_pct, dbh**2 ])
        mu = np.dot(params, xi)
        mu = np.where(mu > dbh, mu, dbh)

        g = dnorm(y, mu, sd=sd)
        # g = dtnorm(y, mu=mu, a0=dbh, sd=self.sd)

        # survival
        if self.no_mort:
            s = 1.00
        elif self.sw_mort == 'model':
            xi = np.vstack([ one, dbh, sw_psd, aw_pba ])
            mu = np.dot(self.survival_params, xi)
            s  = exp(mu) / (1.0 + exp(mu))
        elif self.sw_mort == 'const':
            s = 0.99
        else:
            raise ValueError("invalid sw_mort, should be 'model' or 'const'")

        return s * g


class ABBAW(ABB):


    def __init__(self):
        ABB.__init__(self)

        self.growth_params = np.array(
            [ 3.674, 0.995, -0.0148, 0.0 ])

        self.growth_params_no_comp = np.array(
            [ 2.351, 1.004, 0.0, -1.336 ])

        self.sd = sqrt(3.551841)
        self.sd_no_comp = sqrt(1.000585)


    def set_n0(self, year):

        self.n0 = np.asarray(self.meas['AW'][year])


    def kernel(self, x, y, t, ix=None, **kwargs):

        if ix is not None:
            one     = np.ones(len(x))
            dbh     = x
            # di      = (y - x) / 10.0
            sw_psd  = self._sw_psd[ix]
            # aw_pba  = self._aw_pba[ix]
            # net_pba = self._net_pba[ix]
        else:
            one     = np.ones(len(y))
            dbh     = np.array(len(y) * [x])
            # di      = y - x
            sw_psd  = np.array(len(y) * [ self.sw_psd_from_meas(x)  ])
            # aw_pba  = np.array(len(y) * [ self.aw_pba_from_meas(x)  ])
            # net_pba = np.array(len(y) * [ self.net_pba_from_meas(x) ])


        # growth
        if self.competition:
            params = self.growth_params
            sd     = self.sd
        else:
            params = self.growth_params_no_comp
            sd     = self.sd_no_comp

        xi = np.vstack([ one, dbh, sw_psd, dbh**2 ])
        mu = np.dot(params, xi)
        mu = np.where(mu > dbh, mu, dbh)

        g = dnorm(y, mu, sd=sd)
        # g = dtnorm(y, mu=mu, a0=dbh, sd=self.sd)

        # # survival
        # xi = np.vstack([ one, dbh, dbh**2, di, dbh**2 * net_pba ])
        # mu = np.dot(self.survival_params, xi)

        # s  = exp(mu) / (1.0 + exp(mu))

        if self.no_mort:
            s = 1.00
        else:
            s = 0.99

        return s * g
