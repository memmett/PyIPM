
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
ignore = []
keep = []


for name in names:

    p = psp[ psp["plot"] == name ]

    # cache plotsize
    psizes = list(set(p["plotsize"]))
    if len(psizes) > 1:
        continue
    plotsizes.append({ 'plot': int(name), 'plotsize': psizes[0] })


    # check measurement years
    years = sorted(set(p["year"]))
    if len(years) < 2:
        continue

    # check species composition
    fy  = min(years)
    pfy = p[ p["year"] == fy ]
    nfy = len(pfy)

    species = set(pfy["spec"])

    # skewed = False
    # for spec in species:
    #     if spec in [ SW, AW ]:
    #         continue
    #     nspec = len(pfy[ pfy["spec"] == spec ])
    #     if float(nspec) / nfy > 0.3:
    #         # print spec, nspec, nfy
    #         skewed = True

    # if skewed:
    #     print "skipping %s: too many trees of other species present" % name
    #     ignore.append(name)
    #     continue

    swfy = pfy[ pfy["spec"] == SW ]
    awfy = pfy[ pfy["spec"] == AW ]

    nsw = len(swfy)
    naw = len(awfy)

    if float(nsw)/nfy < 0.2 or float(naw)/nfy < 0.2:
        print "skipping %s: not enough spruce and/or aspen present" % name
        ignore.append(str(int(name)))
        continue

    print 'kept: %s' % name
    keep.append(str(int(name)))

    # # check mean size of sw at first year
    # mdbh = swfy["dbh"].mean()
    # if mdbh < 160.0:
    #     continue

    # grab tree ids at first year and filter
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

print len(ignore), len(keep)

plotsizes = pandas.DataFrame(plotsizes)
plotsizes.to_csv("plotsizes.csv", index=False)

