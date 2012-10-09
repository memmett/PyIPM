
from numpy import sqrt, exp, pi
from scipy.special import erf

sqrt2       = sqrt(2.0)
inv_sqrt2pi = 1.0/sqrt(2.0*pi)


def dnorm(x, mu=0.0, sd=1.0):
    c = 1.0/sd
    x = c * (x - mu)
    return inv_sqrt2pi * c * exp(-x**2/2.0)


def pnorm(x, mu=0.0, sd=1.0):
    c = 1.0/(sqrt2*sd)
    x = c * (x - mu)
    return 0.5 * (1.0 + erf(x))


def dtnorm(x, a=0.0, mu=0.0, sd=1.0):
    c = 1.0/sd
    x = c * (x - mu)
    return dnorm(x)/(1.0 - pnorm((a - mu)/sd))/sd * (x >= a)

