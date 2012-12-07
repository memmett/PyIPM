"""Alberta boreal SW and AW competition experiment."""


import numpy as np
import cPickle as pickle
import pylab as plt

from exp_abb_kernels import abb_init_kernels


with open('abb_sw_mort_const.pkl', 'r') as f:
    plots = pickle.load(f)

# with open('abb_sw_mort_model.pkl', 'r') as f:
#     plots = pickle.load(f)


###############################################################################
# histogram of annual mortality rates

mrates = []

for plotname in plots:
    plot = plots[plotname]
    L, U, N, T, x, plotfile = plot['attrs']

    sw,  aw  = abb_init_kernels(L, U, N, 'sw_mort_model', plotname)
    nsw, naw = plot['nsw'], plot['naw']

    for j in range(1, nsw.shape[0]-1):
        mr = sw.population(nsw[j+1]) / sw.population(nsw[j])
        if mr <= 1.0:
            mrates.append(mr)

mrates = np.asarray(mrates)
print np.mean(mrates)

plt.hist(mrates, bins=101)


###############################################################################
# observed avg mortality rate vs predicted avg mortality rate for model

obs = []
prd = []

for plotname in plots:

    plot = plots[plotname]
    L, U, N, T, x, plotfile = plot['attrs']

    sw,  aw  = abb_init_kernels(L, U, N, 'sw_mort_model', plotname)
    nsw, naw = plot['nsw'], plot['naw']

    # fy = sw.years[0]
    # if aw.meas['AW'][fy] > sw.meas['SW'][fy]:
    #   continue


    for i, prev_obs_year in enumerate(sw.years[:-1]):
        obs_year = sw.years[i+1]
        nyears   = obs_year - prev_obs_year

        p0 = len(sw.meas['SW'][prev_obs_year]) / sw.plot_size * 1e4
        p1 = len(sw.meas['SW'][obs_year])      / sw.plot_size * 1e4

        # print nyears, p0, p1, (p1/p0)**(1.0/nyears)

        if p1 > p0:
            print 'NEW TREES:', plotname, obs_year
            continue

        j = T.index(obs_year)
        p2 = sw.population(nsw[j])

        obs.append(p1/p0)
        prd.append(p2/p0)


from utils.stats import dent_blackie

F, pval, (a, b) = dent_blackie(obs, prd)
print '  F: ', F, pval

x    = np.asarray(sorted(prd))
yhat = a + b*x
yex  = x

plt.figure()
plt.plot(prd, obs, 'ok', label='data', alpha=0.5)
plt.plot(x,   yex, '-r', label='exact')
plt.plot(x,   yhat, '-b', label='fit')
plt.ylabel('observed')
plt.xlabel('predicted')
#plt.axis('equal')
plt.legend(loc='best')

plt.show()


