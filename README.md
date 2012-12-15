PyIPM
=====

Python IPM framework.

This repository also contains numerical experiments and plotting
routines for several projects:

1. The 'A guide to efficient integration schemes for integral
   projection models (IPMs)' paper analysing the accuracy and
   efficiency of various integration schemes by Andria Dawson and
   Matthew Emmett.

   To run all of the numerical experiments for this paper:

   $ python exp_eff_run.py

   To generate the plots:

   $ python exp_eff_plot.py


2. The 'ARTR' project by the Clark Lab at Duke University.

   To run the ARTR experiment:

   $ python exp_artr.py

3. The 'ABB' project consists of several steps:

   To run all of the experiments:

   $ python exp_abb_run.py

   To generate various plots and statistics:

   $ python exp_abb_eviction.py   # compute eviction rates
   $ python exp_abb_comp.py       # compare with/without competition in growth model
   $ python exp_abb_ks.py         # generate ks plots for with/without competition models
   $ python exp_abb_mort.py       # compare modeled vs const mortality model
   $ python exp_abb_obspred.py    # generate obs/pred plots and statistics

