"""Alberta boreal SW and AW competition experiment.

Using "sw_mort_const" and "sw_mort_model" flavoured runs, compare
effect of constant vs modeled spruce mortality.

"""


import numpy as np
import cPickle as pickle
import pylab as plt

from pprint import pprint

from exp_abb_kernels import abb_init_kernels
from utils.stats import dent_blackie


###############################################################################
# load runs

with open('out/abb_sw_mort_const.pkl', 'r') as f:
    plots_const = pickle.load(f)

with open('out/abb_sw_mort_model.pkl', 'r') as f:
    plots_model = pickle.load(f)


###############################################################################
# histogram of annual mortality rates (model only)

mrstats = { 'model': {}, 'const': {} }

mrates = []

for plotname, plot in plots_model.iteritems():
    L, U, N, T, x, plotfile = plot['attrs']

    sw,  aw  = abb_init_kernels(L, U, N, 'sw_mort_model', plotname)
    nsw, naw = plot['nsw'], plot['naw']

    for j in range(1, nsw.shape[0]-1):
        mr = sw.population(nsw[j+1]) / sw.population(nsw[j])
        if mr <= 1.0:
            mrates.append(mr)

mrates = np.asarray(mrates)

mrstats['model']['mean'] = np.mean(mrates)
mrstats['model']['var']  = np.var(mrates)

plt.figure()
plt.hist(mrates, bins=501)
plt.xlim([0.9, 1.0])
plt.xlabel('mortality rate')
plt.ylabel('frequency')
plt.savefig('plots/mr_hist.pdf')

print "histogram saved to: plots/mr_hist.pdf"


###############################################################################
# observed avg mortality rate vs predicted avg mortality rate for model

obs = []
prd = { 'model': [], 'const': [] }

for plotname in plots_const:

    plot1 = plots_const[plotname]
    plot2 = plots_model[plotname]

    L, U, N, T, x, plotfile = plot1['attrs']
    sw,  aw  = abb_init_kernels(L, U, N, 'sw_mort_model', plotname)

    nsw1, naw1 = plot1['nsw'], plot1['naw']
    nsw2, naw2 = plot2['nsw'], plot2['naw']

    for i, prev_obs_year in enumerate(sw.years[:-1]):
        obs_year = sw.years[i+1]
        nyears   = obs_year - prev_obs_year

        p0 = len(sw.meas['SW'][prev_obs_year]) / sw.plot_size * 1e4
        p1 = len(sw.meas['SW'][obs_year])      / sw.plot_size * 1e4

        # XXX: should really remove these plots
        if p1 > p0:
            print 'NEW TREES:', plotname, obs_year
            continue

        j = T.index(obs_year)
        n1 = sw.population(nsw1[j])
        n2 = sw.population(nsw2[j])

        obs.append(p1/p0)
        prd['const'].append(n1/p0)
        prd['model'].append(n2/p0)


for flavour in prd:
    F, pval, (a, b) = dent_blackie(obs, prd[flavour])
    mrstats[flavour]['F'] = (F, pval)

    x    = np.asarray(sorted(prd[flavour]))
    yhat = a + b*x
    yex  = x

    plt.figure()
    plt.plot(prd[flavour], obs, 'ok', label='data', alpha=0.5)
    plt.plot(x, yex,  '-r', label='exact')
    plt.plot(x, yhat, '-b', label='fit')
    plt.ylabel('observed')
    plt.xlabel('predicted')
    plt.legend(loc='best')
    plt.savefig('plots/mr_obsprd_%s.pdf' % flavour)
    print "obs vs prd saved to: plots/mr_obsprd_%s.pdf" % flavour


print "mortality stats:"
pprint(mrstats)

plt.show()