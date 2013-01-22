"""Alberta boreal PSP statistics."""


import numpy as np
import pylab as plt
import cPickle as pickle
from exp_abb_kernels import ignore
import xlrd


book   = xlrd.open_workbook("kernels/abb/PSPLatLong.xls")
sdsloc = book.sheet_by_name("SDSLOC")

plotmap = {}
for row in range(sdsloc.nrows):
    v = sdsloc.cell_value(row, 1)
    if isinstance(v, float):
        plotmap[int(v)] = row

with open('out/abb_with_comp.pkl', 'r') as f:
    plots = pickle.load(f)


print "pspplot,lat,lon,nsr"

for plotname in plots:
    if plotname in ignore:
        continue

    row = plotmap[int(plotname)]

    lat, lon = sdsloc.cell_value(row, 21), sdsloc.cell_value(row, 20)
    nsr = int(sdsloc.cell_value(row, 5))

    if plotname == '603':
        lat, lon = 53.2958205, 117.8706427

    print "%s,%s,-%s,%s" % (plotname, lat, lon, nsr)



