"""Zuidema's method."""

import logging

import numpy as np
import base

class MidPointZuidema(base.Method):

    name = 'Zuidema'
    mesh_type = 'cell'
    sub_mesh_size = 200


    def sample(self, kernel, t):

        N = self.N
        U = self.U
        L = self.L
        M = self.sub_mesh_size
        X = self.X

        h = (U - L) / N / M

        A = np.zeros((N, N))
        for i in range(N):
            logging.debug('ZUIDEMA: sampling (%d): %d', N, i)
            y = np.arange(X[i]+h/2, X[i+1], h)
            for j in range(N):
                x = np.arange(X[j]+h/2, X[j+1], h)
                A[i, j] = np.sum(kernel.kernel(x, y, t))

        self.A = h * A

        logging.debug("ZUIDEMA: update: kernel sampled")
