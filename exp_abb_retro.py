"""Alberta boreal retro-spective data statistics."""

import pylab as plt
import numpy as np

import pandas

sw = pandas.read_csv('spruce_growth_data.csv').dropna()
aw = pandas.read_csv('aspen_growth_data.csv').dropna()

density = pandas.read_csv('plot.density')

idx = sw.live == 1

print 'nsw', len(set(sw.tree[idx]))
print 'naw', len(set(aw.tree))

td = pandas.read_csv('tree_data.csv').dropna()
idx = td.live == 0

print len(set(td.tree[idx]))
print len(set(td.tree))

print len(set(td.tree)), len(set([x for x in td.tree if x[2] == 'L'])), len(set([x for x in td.tree if x[2] == 'D']))

print len(td.tree)


#raise SystemExit



# dbh histograms
fig, ax = plt.subplots(1, 2, sharey=True)

ax[0].hist(sw.dbh, 50, color='0.5')
ax[0].set_title('Spruce')
ax[0].set_ylabel('frequency')
ax[0].set_xlabel('dbh (mm)')

ax[1].hist(aw.dbh, 50, color='0.5')
ax[1].set_title('Aspen')
ax[1].set_xlabel('dbh (mm)')

fig.savefig('plots/retrodbhhist.png')


# sph plots

sph = {}
sph['sw'] = density.spruce
sph['aw'] = density.hw

plt.figure()
idx = np.argsort(sph['sw'])
x   = np.arange(len(idx))

width = 0.35

plt.bar(x, sph['sw'][idx], width, color='black', label='spruce')
plt.bar(x+width, sph['aw'][idx], width, color='white', hatch='/', label='aspen')

plt.xlabel('psp')
plt.ylabel('sph (#/ha)')
plt.legend(loc='best')

plt.savefig('plots/retrosph.png')


plt.show()
