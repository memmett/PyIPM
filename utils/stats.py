
import numpy as np
import scipy.stats

from numpy import sqrt, exp, pi
from scipy.special import erf, erfc, smirnov
from scipy import polyfit, polyval
from numpy.linalg import norm


sqrt2       = sqrt(2.0)
inv_sqrt2pi = 1.0/sqrt(2.0*pi)
sqrt2sqrtpi = sqrt(2.0/pi)

def dnorm(x, mu=0.0, sd=1.0):
    """Normal distribution."""
    c = 1.0/sd
    x = c * (x - mu)
    return inv_sqrt2pi * c * exp(-0.5 * x**2)


def pnorm(x, mu=0.0, sd=1.0):
    c = 1.0/(sqrt2*sd)
    x = c * (x - mu)
    return 0.5 * (1.0 + erf(x))


def dtnorm(x0, a0=0.0, mu=0.0, sd=1.0):
    """Truncated normal distribution."""
    c = 1.0/sd
    x = c * (x0 - mu)
    a = c * (a0 - mu) / sqrt2
    return c * sqrt2sqrtpi * exp(-0.5 * x**2) / erfc(a) * (x0 > a0)


def ks1(data, model, x, ord=np.inf):
    vals = np.sort(data)

    # build the continuous cdf
    total = np.sum(model)
    ccdf  = np.zeros(vals.shape)
    for i in xrange(len(vals)):
        idx = x <= vals[i]
        ccdf[i] = np.sum(model[idx]) / total

    # build the discrete cdf
    N = len(vals)
    dcdf = np.cumsum(np.ones(vals.shape)) / N

    d = norm(dcdf-ccdf, ord)
    p = smirnov(N, d)

    return ccdf, dcdf, d, p


def dent_blackie(obs, prd):

  obs = np.asarray(obs)
  prd = np.asarray(prd)

  n    = len(obs)
  b, a = polyfit(prd, obs, 1)
  yhat = polyval([b, a], prd)

  top = n * a**2 + 2 * a * (b - 1) * sum(prd) + (b - 1)**2 * sum(prd**2)
  bot = 2 * sum((obs-yhat)**2) / (n - 2)

  F = top / bot

  return F, scipy.stats.f.sf(F, n-2, 2), (a, b)


def theil(obs, prd):

  obs = np.asarray(obs)
  prd = np.asarray(prd)
  n   = len(obs)

  top = np.sqrt(1.0/n * sum((prd-obs)**2))
  bot = np.sqrt(1.0/n * sum(prd**2)) + np.sqrt(1.0/n * sum(obs**2))

  U = top / bot

  return U


def MSEP(obs, prd):

  from scipy.stats import pearsonr

  obs = np.asarray(obs)
  prd = np.asarray(prd)
  n   = len(obs)

  mu_x = np.mean(prd)
  sd_x = np.std(prd)

  mu_y = np.mean(obs)
  sd_y = np.std(obs)

  r, pval = pearsonr(obs, prd)

  msep = 1.0/n * sum((prd-obs)**2)

  mc = (mu_x - mu_y)**2 / msep
  sc = (sd_x - r * sd_y)**2 / msep
  rc = (1.0 - r**2) * sd_y**2 / msep

  return mc, sc, rc, msep


