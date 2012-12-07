"""Alberta boreal SW and AW competition experiment."""


import numpy as np
import cPickle as pickle
import pylab as plt

from exp_abb_kernels import abb_init_kernels


###############################################################################
# load results

with open('abb_with_nocomp.pkl', 'r') as f:
    plots = pickle.load(f)


###############################################################################
# main

pops = {
    'obs': { 'sw': [], 'aw': [] },
    'prd': { 'sw': [], 'aw': [] },
}

for plotname in plots:

    plot = plots[plotname]

    L, U, N, T, x, plotfile = plot['attrs']
    nsw, naw = plot['nsw'], plot['naw']

    # sw = ABBSW()
    # aw = ABBAW()

    # for k in [ sw, aw ]:
    #     k.L, k.U = L, U
    #     k.setup(MidPoint(), N)
    #     k.measurements(plotfile, plotname)

    # sw, aw = abb_init_kernels(L, U, N, 'sw_mort_model', plotname)
    sw, aw = abb_init_kernels(L, U, N, 'with_comp', plotname)

    # XXX
    fy = sw.years[0]
    if aw.meas['AW'][fy] > sw.meas['SW'][fy]:
      continue

    for obs_year in sw.years[1:]:
        pops['obs']['sw'].append(len(sw.meas['SW'][obs_year]) / sw.plot_size * 1e4)
        pops['obs']['aw'].append(len(aw.meas['AW'][obs_year]) / sw.plot_size * 1e4)

        j = T.index(obs_year)
        pops['prd']['sw'].append(sw.population(nsw[j]))
        pops['prd']['aw'].append(aw.population(naw[j]))


for species in [ 'sw', 'aw' ]:

    from utils.stats import dent_blackie, theil, MSEP
    from scipy.stats import ttest_ind, mannwhitneyu, pearsonr

    print 'species:', species

    obs = pops['obs'][species]
    prd = pops['prd'][species]


    F, pval, (a, b) = dent_blackie(obs, prd)
    print '  F: ', F, pval

    U = theil(obs, prd)
    print '  U: ', U

    tstat, pval = ttest_ind(prd, obs)
    print '  t: ', tstat, pval

    u, pval = mannwhitneyu(prd, obs)
    print '  u: ', u, pval

    p, pval = pearsonr(prd, obs)
    print '  p: ', p, pval

    mc, sc, rc, msep = MSEP(obs, prd)
    print '  msep: ', msep
    print '    mc: ', mc
    print '    sc: ', sc
    print '    rc: ', rc
    print '   sum: ', mc + sc + rc


    x    = np.asarray(sorted(prd))
    yhat = a + b*x
    yex  = x

    plt.figure()
    plt.plot(prd, obs, 'ok', label='data', alpha=0.5)
    plt.plot(x,   yex, '-r', label='exact')
    plt.plot(x,   yhat, '-b', label='fit')
    plt.ylabel('observed')
    plt.xlabel('predicted')
    plt.legend(loc='best')
    plt.title({'sw': 'Spruce', 'aw': 'Aspen'}[species])
    plt.savefig('plots/obspred_%s.png' % species)



plt.show()