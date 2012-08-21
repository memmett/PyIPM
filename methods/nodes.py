"""SDC quadrature matrix generation routines.

:author: Matthew Emmett <matt@emmett.ca>

"""

# Copyright (c) 2011, 2012.  Matthew Emmett.  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#   1. Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#   2. Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


import numpy as np
import sympy
import sympy.mpmath as mpmath


###############################################################################
# polynomial generator, roots etc

def legendre_poly(n):
  """Return Legendre polynomial :math:`P_n(x)`.

  :param n: polynomial degree
  """

  x = sympy.var('x')
  p = (1.0*x**2 - 1.0)**n

  top = p.diff(x, n)
  bot = 2**n * 1.0*sympy.factorial(n)

  return (top / bot).as_poly()


def find_roots(p):
  """Return list of the roots of the polynomial *p*."""

  return sorted(p.nroots(n=100))


###############################################################################
# quadrature points

def gauss_legendre_nodes(n):
  """Return Gauss-Legendre nodes.

  Gauss-Legendre nodes are roots of: :math:`P_n(x)`.
  """

  p = legendre_poly(n)
  return find_roots(p)


def gauss_lobatto_nodes(n):
  """Return Gauss-Lobatto nodes.

  Gauss-Lobatto nodes are roots of: :math:`P'_{n-1}(x)`, and -1, 1.
  """

  x = sympy.var('x')
  p = legendre_poly(n-1).diff(x)
  r = find_roots(p)
  r = [mpmath.mpf('-1.0'), mpmath.mpf('1.0')] + r

  return sorted(r)


def gauss_radau_nodes(n):
  """Return Gauss-Radau nodes.

  Gauss-Radau nodes are: -1 times the roots of :math:`P_n(x) + P_{n-1}(x)`.
  """

  if n == 1:
    r = [ -1.0 ]
  else:
    p = legendre_poly(n) + legendre_poly(n-1)
    r = find_roots(p)

  r = [ -1.0*x for x in r ]
  return sorted(r)


def clenshaw_curtis_nodes(n):
  """Return Clenshaw-Curtis nodes (actually returns n+1 nodes)."""

  r = set()
  for k in range(n+1):
    r.add(mpmath.cos(k*mpmath.pi/mpmath.mpf(n)))

  return sorted(r)

def clenshaw_curtis_weights(n):
  """Return Clenshaw-Curtis weights (actually the n+1 weights)."""

  w = []
  for k in range(n+1):
    c_k = 1.0 if k % n == 0 else 2.0
    t_k = k*mpmath.pi/mpmath.mpf(n)

    acc = 1.0
    for j in range(1, n/2+1):
      b_j = 2.0 if j < n/2 else 1.0
      acc -= b_j / (4.0*j**2 - 1.0) * mpmath.cos(2*j*t_k)

    w.append(c_k / n * acc)

  return w
