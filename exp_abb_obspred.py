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
        stats[species]['msep'] = msep
        stats[species]['mc']   = mc
        stats[species]['sc']   = sc
        stats[species]['rc']   = rc

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





def mktable(stats, label, caption):

    from string import Template

    class LaTexTemplate(Template):
        delimiter = '@'

    table = LaTexTemplate(r"""
\begin{table}
  \begin{center}
    \begin{tabular}{lrrcrr} \toprule
                  & \multicolumn{2}{c}{Spruce} & \hspace{2em} & \multicolumn{2}{c}{Aspen} \\
                  & statistic  & p-value        &              & statistic & p-value       \\ \midrule
      $F$         & @sw_F     & @sw_F_p         &              & @aw_F     & @aw_F_p       \\
      $t$         & @sw_t     & @sw_t_p         &              & @aw_t     & @aw_t_p       \\
      $Theil$     & @sw_theil & @sw_theil_p     &              & @aw_theil & @aw_theil_p   \\
      $U$         & @sw_mannwhitneyu & @sw_mannwhitneyu_p &  & @aw_mannwhitneyu & @aw_mannwhitneyu_p   \\
      $Pearson r$ & @sw_pearsonr & @sw_pearsonr_p &            & @aw_pearsonr & @aw_pearsonr_p   \\
      $MSEP$      & @sw_msep & @sw_msep_p       &              & @aw_msep & @aw_msep_p   \\
      $MC$      & @sw_mc & @sw_mc_p       &              & @aw_mc & @aw_mc_p   \\
      $SC$      & @sw_sc & @sw_sc_p       &              & @aw_sc & @aw_sc_p   \\
      $RC$      & @sw_rc & @sw_rc_p       &              & @aw_rc & @aw_rc_p   \\
                          \bottomrule
    \end{tabular}
  \end{center}
  \label{@label}
  \caption{@caption}
\end{table}
""")

    statvals = {}

    for species in stats:
        for stat in stats[species]:
            if isinstance(stats[species][stat], tuple):
                val, pval = stats[species][stat]
                val  = "%.3g" % val
                pval = "%.3g" % pval
            else:
                val  = stats[species][stat]
                val  = "%.3g" % val
                pval = '--'

            statvals[species + '_' + stat]        = val
            statvals[species + '_' + stat + '_p'] = pval

    return table.substitute(label=label, caption=caption, **statvals)

###############################################################################
# load results

if __name__ == '__main__':

    with open('out/abb_with_comp.pkl', 'r') as f:
        plots = pickle.load(f)

    stats = obsprd_plot(plots, 'with_comp')
    pprint(stats)

    print mktable(stats, 'tab:obsprdstats', 'Observed vs predicted statistics')


    with open('out/abb_with_nocomp.pkl', 'r') as f:
        plots = pickle.load(f)

    stats = obsprd_plot(plots, 'with_nocomp')
    pprint(stats)

    print mktable(stats, 'tab:obsprdstatsnc', 'Observed vs predicted statistics (no comp)')

    # plt.show()
