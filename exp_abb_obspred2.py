"""Alberta boreal SW and AW competition experiment.

Generate obs vs prd plots for a particular flavour.
"""

import matplotlib
matplotlib.rcParams['font.size'] = 14

import numpy as np
import cPickle as pickle
import pylab as plt

from pprint import pprint

from exp_abb_kernels import abb_init_kernels, ignore


###############################################################################
# plotting routine

def obsprd_plot(plots, flavour):

    pops = {
        'obs': { 'sw': [], 'aw': [] },
        'prd': { 'sw': [], 'aw': [] },
    }

    for plotname in plots:

        if plotname in ignore:
            continue

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

        tmp = np.asarray(pops['obs']['sw'])
        if any(tmp > 15000):
            print plotname


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
        plt.plot(x,   yex, '-r', label='y=x')
        plt.plot(x,   yhat, '-b', label='best fit')
        plt.ylabel('observed')
        plt.xlabel('predicted')
        plt.legend(loc='best')
        plt.savefig('plots/obsprd_%s_%s.png' % (flavour, species))
        plt.close()
        # plt.title({'sw': 'Spruce', 'aw': 'Aspen'}[species])

    return stats





def mktable(comp, ncomp, aw, label, caption):

    from string import Template

    class LaTexTemplate(Template):
        delimiter = '@'

    table = LaTexTemplate(r"""
\begin{table}
  \begin{center}
    \begin{tabular}{lrrcrrcrr} \toprule
 & \multicolumn{2}{c}{Spruce}
 & \hspace{2em}
 & \multicolumn{2}{c}{Spruce}
 & \hspace{2em}
 & \multicolumn{2}{c}{Aspen} \\
 & \multicolumn{2}{c}{(with comp)}
 & \hspace{2em}
 & \multicolumn{2}{c}{without comp}
 & \hspace{2em}
 & \multicolumn{2}{c}{(with and without comp)} \\
 & statistic & p-value &
 & statistic & p-value &
 & statistic & p-value \\ \midrule
 $F$         & @F           & @F_p
   & & @F_nc            & @F_nc_p
   & & @F_aw            & @F_aw_p \\
 $t$         & @t           & @t_p
   & & @t_nc            & @t_nc_p
   & & @t_aw            & @t_aw_p \\
 Theil       & @theil   & @theil_p
   & & @theil_nc        & @theil_nc_p
   & & @theil_aw        & @theil_aw_p \\
 $U$         & @mannwhitneyu & @mannwhitneyu_p
   & & @mannwhitneyu_nc & @mannwhitneyu_nc_p
   & & @mannwhitneyu_aw & @mannwhitneyu_aw_p \\
 Pearson $r$ & @pearsonr     & @pearsonr_p
   & & @pearsonr_nc     & @pearsonr_nc_p
   & & @pearsonr_aw     & @pearsonr_aw_p \\
 MSEP        & @msep         & @msep_p
   & & @msep_nc         & @msep_nc_p
   & & @msep_aw         & @msep_aw_p  \\
 MC          & @mc           & @mc_p
   & & @mc_nc           & @mc_nc_p
   & & @mc_aw           & @mc_aw_p \\
 SC          & @sc           & @sc_p
   & & @sc_nc           & @sc_nc_p
   & & @sc_aw           & @sc_aw_p \\
 RC          & @rc           & @rc_p
   & & @rc_nc           & @rc_nc_p
   & & @rc_aw           & @rc_aw_p \\ \bottomrule
    \end{tabular}
  \end{center}
  \label{@label}
  \caption{@caption}
\end{table}
""")

    statvals = {}

    for stat in comp:
        if isinstance(comp[stat], tuple):
            val, pval = comp[stat]
            val  = "%.3g" % abs(val)
            pval = "%.3g" % pval

            valnc, pvalnc = ncomp[stat]
            valnc  = "%.3g" % abs(valnc)
            pvalnc = "%.3g" % pvalnc

            valaw, pvalaw = aw[stat]
            valaw  = "%.3g" % abs(valaw)
            pvalaw = "%.3g" % pvalaw

        else:
            val  = comp[stat]
            val  = "%.3g" % val
            pval = '--'

            valnc  = ncomp[stat]
            valnc  = "%.3g" % valnc
            pvalnc = '--'

            valaw  = aw[stat]
            valaw  = "%.3g" % valaw
            pvalaw = '--'


        statvals[stat]        = val
        statvals[stat + '_p'] = pval

        statvals[stat + '_nc']   = valnc
        statvals[stat + '_nc_p'] = pvalnc

        statvals[stat + '_aw']   = valaw
        statvals[stat + '_aw_p'] = pvalaw


    return table.substitute(label=label, caption=caption, **statvals)

###############################################################################
# load results

if __name__ == '__main__':

    with open('out/abb_with_comp.pkl', 'r') as f:
        plots = pickle.load(f)

    stats1 = obsprd_plot(plots, 'with_comp')
    pprint(stats1)


    with open('out/abb_with_nocomp.pkl', 'r') as f:
        plots = pickle.load(f)

    stats2 = obsprd_plot(plots, 'with_nocomp')
    pprint(stats2)

    print mktable(stats1['sw'], stats2['sw'], stats1['aw'], 'tab:obsprdstats', 'Observed vs predicted statistics')

    # plt.show()
