"""IPM kernel for the ART TRI data."""

import logging

import numpy as np

from numpy import sqrt, exp, log
from utils.io import read_csv
from utils.stats import dnorm

from base import Kernel

class ARTTRI(Kernel):

    L = -8.0
    U = 18.0
    # T = range(1932, 2006)
    T = range(1932, 1932+25+1)

    name = 'ARTR'

    def __init__(self):
        
        Kernel.__init__(self)

        self.time_dependent = True

        # read climate data and measurements (cm^2)
        climate      = read_csv('kernels/artr/climate.csv')
        measurements = read_csv('kernels/artr/survivalDataARTR.csv')
        years        = sorted(list(set([ int(x.year) for x in measurements ])))
        fecundity_measurements = read_csv('kernels/artr/fecdat.csv')

        log_areas = {}
        for year in years:
            log_areas[year] = [ float(x.logArea) for x in measurements 
                                                 if int(x.year) == year ]
            
        # set attributes
        self.n0        = log_areas[years[0]]
        self.log_areas = log_areas
        self.years     = years
        self.climate   = climate
        self.fec_meas  = fecundity_measurements

        self.mortality_type = 'noexp'
        self.fecundity_type = 'uniform_exp'

        self._tweak_covariat = None
        self._tweak_factor   = 1.0


        # set parameters
        self.params = {

            'growth': {

                'mean': np.array(
            [ -6.134753578,  0.765515413,  0.006574618,  0.362328659 ]),

                'covariance': np.array(
            [[  2.3003264050, -3.551420e-04, -9.540608e-04, -1.202835e-01 ],
             [ -0.0003551420,  3.699543e-04,  1.274024e-06, -6.374569e-05 ],
             [ -0.0009540608,  1.274024e-06,  3.545569e-06,  3.898729e-05 ],
             [ -0.1202834507, -6.374569e-05,  3.898729e-05,  6.350332e-03 ]])
                },

            'survival': {

                'mean': np.array(
            [ -1.482555395,  0.753372235,  0.003454076, -0.001100492 ]),

                'covariance':
                    np.array(
            [[  1.149244e-01, -1.422819e-02, -4.075111e-04,  5.042504e-05 ],
             [ -1.422819e-02,  1.058307e-02,  5.006554e-05, -3.693585e-05 ],
             [ -4.075111e-04,  5.006554e-05,  1.520893e-06, -1.896803e-07 ],
             [  5.042504e-05, -3.693585e-05, -1.896803e-07,  1.351324e-07 ]])
                },

            'fecundity': {

                'mean': np.array(
            [ 3.16420466, 0.06773083, 0.02146190 ]),

                'covariance': 
                    np.array(
            [[  2.900329e-02,  3.790527e-04,  2.113102e-05 ],
             [  3.790527e-04,  1.101385e-04,  8.486999e-06 ],
             [  2.113102e-05,  8.486999e-06,  1.079253e-05 ]])
                }
            }

        # default parameters
        self.growth_params    = self.params['growth']['mean']
        self.survival_params  = self.params['survival']['mean']
        self.fecundity_params = self.params['fecundity']['mean']


    def measurements(self, t):

        logging.debug("ARTR: measurements: %d", t)

        return self.log_areas.get(t, [])

    
    def setup(self, method, N):

        Kernel.setup(self, method, N)

        self.row = 1

        logging.debug("ARTR: setup: climate row: %d", self.row)


    def update(self, t):

        Kernel.update(self, t)

        # set row for looking up climate data
        self.row = max(0, t - int(self.years[0]) - 1)

        # set default parameters
        self.growth_params    = self.params['growth']['mean']
        self.survival_params  = self.params['survival']['mean']
        self.fecundity_params = self.params['fecundity']['mean']

        logging.debug("ARTR: update: t: %d, "
                      + "climate row: %d, "
                      + "params set to means", t, self.row)


    def draw_parameters(self, t):

        # growth
        self.growth_params = np.random.multivariate_normal(
            self.params['growth']['mean'], 
            self.params['growth']['covariance'])

        # survival
        self.survival_params = np.random.multivariate_normal(
            self.params['survival']['mean'], 
            self.params['survival']['covariance'])

        # fecundity
        self.fecundity_params = np.random.multivariate_normal(
            self.params['fecundity']['mean'], 
            self.params['fecundity']['covariance'])

        logging.debug("ARTR: growth params:    %s", self.growth_params)
        logging.debug("ARTR: survival params:  %s", self.survival_params)
        logging.debug("ARTR: fecundity params: %s", self.fecundity_params)

    
    def tweak(self, covariat, factor):

        self._tweak_covariat = covariat
        self._tweak_factor   = factor


    def covariat(self, name):

        v = float(getattr(self.climate[self.row], name))

        if name == self._tweak_covariat:
            v *= self._tweak_factor

        return v


    def survival(self, x, t):

        pptLag1 = self.covariat('pptLag1')

        a0, a1, a2, a3 = self.survival_params

        if self.mortality_type == 'noexp':

            return 1.0/(1.0 + exp(-(a0 + a1*x + a2*pptLag1 + a3*pptLag1*x)))

        elif self.mortality_type == 'exp':

            return 1.0/(1.0 + exp(-(a0 + a1*exp(x) 
                                    + a2*pptLag1 + a3*pptLag1*exp(x))))

        else:
            return self.mortality_type


    def growth(self, x, y, t):

        pptWin2   = self.covariat('pptWin2')
        TmeanSum1 = self.covariat('TmeanSum1')

        a0, a1, a2, a3 = self.growth_params

        m = a0 + a1 * x + a2 * pptWin2 + a3 * TmeanSum1
        v = 1.79

        return dnorm(y, mu=m, sd=sqrt(v))


    def fecundity(self, x, y, t):

        y = np.asarray(y)

        recpptApr     = self.covariat('recpptApr')
        recFebMarSnow = self.covariat('recFebMarSnow')

        a0, a1, a2 = self.fecundity_params

        if self.fecundity_type == 'uniform_exp':

            # this is the one that is in the paper and that Andria wants to use

            rpg = exp(a0 + a1*recpptApr + a2*recFebMarSnow)*exp(x)/40000.0
            fa  = 1.0/0.4 * (-1.5 <= y) * (y <= -1.1)

        elif self.fecundity_type == 'uniform_noexp':

            rpg = (a0 + a1*recpptApr + a2*recFebMarSnow)*exp(x)/40000.0
            fa  = 1.0/0.4 * (-1.5 <= y) * (y <= -1.1)

        elif self.fecundity_type == 'exp':

            rpg = (a0 + a1*recpptApr + a2*recFebMarSnow)*exp(x)/40000.0
            fa  = 0.25 * exp(-0.25*exp(y))

        return rpg * fa


    def kernel(self, x, y, t):

        self.increment_count(x, y)

        s = self.survival(x, t)
        g = self.growth(x, y, t)

        if self.fecundity_type == 'population':
            f = 0.0
        else:
            f = self.fecundity(x, y, t)

        return s*g + f


    def population_fecundity(self, n0):

        recpptApr     = self.covariat('recpptApr')
        recFebMarSnow = self.covariat('recFebMarSnow')
        a0, a1, a2    = self.fecundity_params

        # note: P is the quadrature vector that integrates over the
        # entire domain
        total_cover = np.dot(self.method.P, n0 * exp(self.x))
        recruits    = total_cover * exp(a0 + a1*recpptApr + a2*recFebMarSnow) / 40000.0

        fa = 1.0/0.4 * (-1.5 <= self.x) * (self.x <= -1.1)

        return  recruits * fa

        
        

        
        


