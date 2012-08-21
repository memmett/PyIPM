"""Point-based Gaussian quadrature."""

import numpy as np
import base


def adjust(x, discontinuities):
    """Adjust linearly spaced cell edges contained in x to be
    consistent with the discontinuities."""

    x = np.asarray(x).copy()
    tie_pts = [ x[0] ] + discontinuities + [ x[-1] ]

    # some notes:
    # - i marks the current cell edge
    # - j marks the index of the next cell edge that is close to a
    #   discontinuity, or the end of the array

    i = 0
    for tie in tie_pts:
        j = i + np.argmin(abs(x[i:] - tie))

        if j == i:
            x[i] = tie
        else:
            i = i - 1
            slope    = (tie - x[i]) / (j - i)
            x[i:j+1] = x[i] + slope * np.arange(0, j-i+1, 1)

        i = j + 1

    return x


class GaussQuad(base.Method):

    mesh_type = 'point'

    def __init__(self, qtype='GaussLegendre', k=None, adjust=False):

        self.k = k
        self.qtype = qtype
        self.adjust = adjust

        if k:
            self.name = qtype + '(' + str(k) + ')'
        else:
            self.name = qtype

        if adjust:
            self.name = 'Adj' + self.name


    def setup(self, L, U, N, discontinuities=None, **kwargs):

        super(GaussQuad, self).setup(L, U, N, **kwargs)

        k = self.k
        L = self.L
        U = self.U
        N = self.N

        if k:
            # compute quad weights and points
            if self.qtype == 'GaussLegendre':
                from scipy.special.orthogonal import p_roots
                xi, w = p_roots(k)
            elif self.qtype == 'ClenshawCurtis':
                import nodes
                xi = np.asarray(map(float, nodes.clenshaw_curtis_nodes(k-1)))
                w  = np.asarray(map(float, nodes.clenshaw_curtis_weights(k-1)))

            # compute cell edges and adjust if necessary
            edges = np.linspace(L, U, N+1)
            jumps = getattr(self, 'discontinuities', discontinuities)
            if self.adjust and jumps:
                edges = adjust(edges, jumps)

            # compute quad weights and points
            x = np.zeros(k*N)
            P = np.zeros(k*N)

            for i in range(N):
                dx2 = 0.5 * (edges[i+1] - edges[i])
                for l in range(k):
                    x[i*k+l] = edges[i] + dx2 * (xi[l] + 1.0)
                    P[i*k+l] = dx2 * w[l]

            self.P = P
            self.x = x

        else:

            # compute quad weights and points
            if self.qtype == 'GaussLegendre':
                from scipy.special.orthogonal import p_roots
                xi, w = p_roots(N)
            elif self.qtype == 'ClenshawCurtis':
                import nodes
                xi = np.asarray(map(float, nodes.clenshaw_curtis_nodes(N-1)))
                w  = np.asarray(map(float, nodes.clenshaw_curtis_weights(N-1)))

            self.P = w * (U - L)/2.0
            self.x = (U + L)/2 + (U - L)/2 * xi


    def sample(self, kernel, t):
        """Update the integration matrix A."""

        N = len(self.P)
        A = np.empty((N, N))
        for j in range(len(self.x)):
            A[j, :] = self.P * kernel.kernel(self.x, self.x[j], t)

        self.A = A
