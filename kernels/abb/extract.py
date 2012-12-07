
import numpy as np
import pandas

psp = np.load('psp.npz')
psp = { k: pandas.Series(v) for k, v in psp.iteritems() }
psp = pandas.DataFrame(psp)
psp = psp.dropna()

names = set(psp["plot"])

SW = 18                  # integer level corresponding to R 'SW' level
AW = 3                   # integer level corresponding to R 'AW' level

plotsizes = []

for name in names:
    p = psp[ psp["plot"] == name ]
    years = sorted(set(p["year"]))

    if len(years) < 2:
        continue

    # cache plotsize
    psizes = list(set(p["plotsize"]))
    if len(psizes) > 1:
        continue
    plotsizes.append({ 'plot': int(name), 'plotsize': psizes[0] })

    # check number of sw and aw for first year
    fy = min(years)
    swfy = p[ (p["year"] == fy) & (p["spec"] == SW) ]
    awfy = p[ (p["year"] == fy) & (p["spec"] == AW) ]

    nsw = len(set(swfy["tree"]))
    naw = len(set(awfy["tree"]))

    if naw > nsw or nsw < 30 or naw < 10:
        continue

    # check mean size of sw at first year
    mdbh = swfy["dbh"].mean()
    if mdbh < 160.0:
        continue

    # grad tree ids at first year and filter
    swtrees = list(set(swfy["tree"]))
    awtrees = list(set(awfy["tree"]))

    sw = p[ (p["spec"] == SW) & (p["tree"].isin(swtrees)) ]
    aw = p[ (p["spec"] == AW) & (p["tree"].isin(awtrees)) ]


    # save!
    sw["spec"] = 'SW'
    aw["spec"] = 'AW'

    plot = pandas.concat([ sw, aw ])
    plot["year"] = np.asarray(plot["year"], np.int16)
    del plot['plotsize']
    del plot['tree']
    del plot['plot']
    plot.to_csv(str(int(name)) + '.csv', index=False)

plotsizes = pandas.DataFrame(plotsizes)
plotsizes.to_csv("plotsizes.csv", index=False)

