"""Base class for kernels."""

import numpy as np


class Kernel(object):


    def __init__(self):
        self.count = 0
        self.method = None
        self.time_dependent = False


    def reset_count(self):
        self.count = 0


    def increment_count(self, x, y):
        if isinstance(x, float):
            x = [ x ]
        if isinstance(y, float):
            y = [ y ]
        self.count += len(x) + len(y)


    def setup(self, method, N, **kwargs):
        if hasattr(self, 'discontinuities'):
            method.discontinuities = self.discontinuities
        self.method = method
        self.method.setup(self.L, self.U, N, **kwargs)
        self.x = self.method.x


    def update(self, t, **kwargs):
        if self.time_dependent or self.method.A is None:
            self.method.sample(self, t)

    
    def first_projection(self):
        if callable(self.n0):
            n0 = self.n0
        else:
            n0 = self.empirical_distribution(self.n0)
        return self.method.eval_or_integrate(n0)


    def empirical_distribution(self, measurements):
        k = self.kernel
        f = lambda y: sum(map(lambda x: k(x, y, 0), measurements))
        return f


    def project(self, n0):
        return np.dot(self.method.A, n0)


    def population(self, n0):
        return np.dot(self.method.P, n0)


    def growth_rate(self, n1, n2):
        """Compute the growth rate between *n1* and *n2*."""

        p1 = self.population(n1)
        p2 = self.population(n2)

        return p2 / p1


    @property
    def projection_matrix(self):
        return self.method.A
