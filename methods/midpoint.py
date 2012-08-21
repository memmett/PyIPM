"""Mid-point method."""

import logging

import numpy as np
import base

class MidPoint(base.Method):

    name = 'MidPoint'
    mesh_type = 'cell'

    def sample(self, kernel, t):

        N = self.N
        U = self.U
        L = self.L
        x = self.x

        logging.debug("MIDPOINT: sampling kernel")

        dx = (U - L) / N

        A = np.zeros((N, N))
        for i in range(N):
            A[i, :] = dx * kernel.kernel(x, x[i], t)

        self.A = A

