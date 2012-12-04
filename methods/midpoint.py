"""Mid-point method."""

import logging

import numpy as np
import base

class MidPoint(base.Method):

    name = 'MidPoint'
    mesh_type = 'cell'

    def eval_or_integrate(self, f):

        return f(self.x) / self.dx


    def sample(self, kernel, t):

        N = self.N
        U = self.U
        L = self.L
        x = self.x

        logging.debug("MIDPOINT: sampling kernel")

        dx = (U - L) / N

        ix = np.array(range(N))

        A = np.zeros((N, N))
        for i in xrange(N):
            A[i, :] = dx * kernel.kernel(x, x[i], t, ix=ix, iy=i)

        self.A = A

