"""IPM integration method base class."""

import logging

import numpy as np

class Method(object):

    def eval_or_integrate(self, f):
        """Evaluate *f* at mesh points, or integrate *f* over mesh cells."""

        if self.mesh_type == 'cell':

            import scipy.integrate
            n = np.zeros(self.N)
            for i in range(self.N):
                # n[i], _ = scipy.integrate.quad(f, self.X[i], self.X[i+1])
                n[i], _ = scipy.integrate.fixed_quad(f, self.X[i], self.X[i+1])
                n[i] /= self.dx

        else:

            n = np.zeros(len(self.x))
            n[...] = f(self.x)

        return n


    def histogram(self, n0):
        """Compute histogram of *n0* based on cells."""

        assert self.mesh_type == 'cell'

        n, edges = np.histogram(n0, self.X)

        return n


    def total_population(self, n):
        """Compute the total population."""

        if self.mesh_type == 'cell':
            return sum(n) * self.dx

        return np.dot(self.P, n)


    def growth_rate(self, n1, n2):
        """Compute the growth rate between *n1* and *n2*."""

        p1 = self.total_population(n1)
        p2 = self.total_population(n2)

        return p2 / p1


    def setup(self, L, U, N, *args, **kwargs):

        self.L = L
        self.U = U
        self.N = N
        self.A = None
        self.P = None

        if self.mesh_type == 'cell':
            dx = (U - L) / N

            X = np.linspace(L, U, N+1)
            x = X[:-1] + dx/2

            self.X  = X
            self.x  = x
            self.dx = dx
            self.P  = dx * np.ones(N)


