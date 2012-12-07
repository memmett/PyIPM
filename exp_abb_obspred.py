"""Alberta boreal SW and AW competition experiment.

Generate obs vs prd plots for a particular flavour.
"""


import numpy as np
import cPickle as pickle
import pylab as plt

from pprint import pprint

from exp_abb_kernels import abb_init_kernels


###############################################################################
# plotting routine

def obsprd_plot(plots, flavour):

    pops = {
        'obs': { 'sw': [], 'aw': [] },
        'prd': { 'sw': [], 'aw': [] },
    }

    for plotname in plots:

        plot = plots[plotname]

        L, U, N, T, x, plotfile = plot['attrs']
        nsw, naw = plot['nsw'], plot['naw']

        sw, aw = abb_init_kernels(L, U, N, flavour, plotname)

        for obs_year in sw.years[1:]:
            pops['obs']['sw'].append(len(sw.meas['SW'][obs_year]) / sw.plot_size * 1e4)
            pops['obs']['aw'].append(len(aw.meas['AW'][obs_year]) / sw.plot_size * 1e4)

            j = T.index(obs_year)
            pops['prd']['sw'].append(sw.population(nsw[j]))
            pops['prd']['aw'].append(aw.population(naw[j]))


    stats = { 'sw': {}, 'aw': {} }

    for species in [ 'sw', 'aw' ]:

        from utils.stats import dent_blackie, theil, MSEP
        from scipy.stats import ttest_ind, mannwhitneyu, pearsonr

        obs = pops['obs'][species]
        prd = pops['prd'][species]

        F, pval, (a, b) = dent_blackie(obs, prd)
        stats[species]['F'] = (F, pval)

        stats[species]['theil'] = theil(obs, prd)
        stats[species]['t'] = ttest_ind(prd, obs)
        stats[species]['mannwhitneyu'] = mannwhitneyu(prd, obs)
        stats[species]['pearsonr'] = pearsonr(prd, obs)

        mc, sc, rc, msep = MSEP(obs, prd)
        stats[species]['msep_msep'] = msep
        stats[species]['msep_mc']   = mc
        stats[species]['msep_sc']   = sc
        stats[species]['msep_rc']   = rc

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
        plt.savefig('plots/obsprd_%s_%s.png' % (flavour, species))
        plt.close()
        # plt.title({'sw': 'Spruce', 'aw': 'Aspen'}[species])

    return stats


###############################################################################
# load results

if __name__ == '__main__':

    with open('out/abb_sw_mort_model.pkl', 'r') as f:
        plots = pickle.load(f)

    stats = obsprd_plot(plots, 'sw_mort_model')
    pprint(stats)

    plt.show()





