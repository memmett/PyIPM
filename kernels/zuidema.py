"""Zuidema Parashorea Chinensis kernel."""

from numpy import sqrt, exp, dot, zeros, asarray
from utils.stats import dnorm, dtnorm
from base import Kernel


class Zuidema(Kernel):

    name = 'Zuidema'
    time_dependent = False

    def __init__(self):

        Kernel.__init__(self)

        self.L =   1.0
        self.U = 150.0
        self.T = range(0, 5)

        # parashorea chinensis
        self.survival_params = [ 0.98 ]
        self.growth_params   = [ 42.1, 144.0, 2.258, 0.1054 ] # a, b, c, sd
        self.k_st_params     = [ -3.232, 0.146, 0.829 ]       # b, mu, sdl
        self.k_ts_params     = [ 0.118, 1.078, 0.369 ]        # a54, mu, sd

        self.k_ss = asarray([[ 0.700, 0.0,   0.0,   0.0   ],
                                [ 0.101, 0.704, 0.0,   0.0   ],
                                [ 0.0,   0.136, 0.793, 0.0   ],
                                [ 0.0,   0.0,   0.096, 0.819 ]])


    def n0(self, x):

        return dnorm(x, mu=50.0, sd=sqrt(10.0))


    #### continuous part

    def survival(self, x, t):

        s = self.survival_params

        return s[0]


    def growth(self, x, y, t):

        a, b, c, sd = self.growth_params

        m = b * c * x**(c-1.0) / ( b + x**c/a )**2

        return dnorm(y, mu=x+m, sd=sd)


    def kernel(self, x, y, t):

        self.increment_count(x, y)

        s = self.survival(x, t)
        g = self.growth(x, y, t)

        return s*g


    #### discrete parts

    def k_st(self, x):

        N = self.N
        A = zeros((4, N))

        s = self.survival_params[0]
        b0, b1, kids = self.k_st_params
        p = 1.0/(1.0 + exp(-(b0 + b1 * x)))

        A[0, :] = s * p * kids * self.method.P

        return A


    def k_ts(self, x):

        L = self.L
        N = self.N
        A = zeros((N, 4))

        a54, mu, sd = self.k_ts_params
        
        A[:, -1] = a54 * dtnorm(x, L, mu, sd)

        return A


    def setup(self, method, N):

        self.method = method
        self.method.setup(self.L, self.U, N)
        self.N = len(self.method.x)

        N = self.N
        self.x = zeros(4+N)
        self.x[4:] = self.method.x
        self.x[:4] = asarray([ -4, -3, -2, -1 ])

        self.method.sample(self, 0) # sample continuous part of myself

        self.A = zeros((N+4, N+4))
        self.A[:4, :4] = self.k_ss
        self.A[:4, 4:] = self.k_st(self.method.x)
        self.A[4:, :4] = self.k_ts(self.method.x)
        self.A[4:, 4:] = self.method.A


    def population(self, n):

        return sum(n[:4]) + dot(self.method.P, n[4:])


    def first_projection(self):

        n1 = zeros(4+self.N)

        if callable(self.n0):
            n0 = self.n0
        else:
            n0 = self.empirical_distribution(self.n0)

        n1[4:] = self.method.eval_or_integrate(n0)

        return n1
        

    def project(self, n0):
        return dot(self.A, n0)


    @property
    def projection_matrix(self):
        return self.A


    #### y-integrable form

    def s(self, x, t):

        return self.survival(x, t)


    def mu(self, x, t):

        self.increment_count(x, [])

        a, b, c, sd = self.growth_params

        m = b * c * x**(c-1.0) / ( b + x**c/a )**2

        return x + m


    def sd(self, x, t):

        a, b, c, sd = self.growth_params

        return sd


    def r(self, x, t):

        return 0.0


    def F(self, y1, y2, t):

        return 0.0
