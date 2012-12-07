"""Grab a few columns from pspksclean.rdata and save in a Python pickle."""


import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
R = ro.r

import numpy as np

psp = {}

R("load('pspksclean.rdata')")
for column in [ 'plot', 'year', 'tree', 'dbh', 'spec', 'plotsize' ]:
    R("tmp = psp$" + column.upper())
    psp[column] = np.asarray(R["tmp"])

np.savez_compressed('psp.npz', **psp)