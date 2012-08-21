#
# Generate various error plots.
#

import cPickle as pickle
import matplotlib as mpl
mpl.rc('legend', fontsize='small')

import matplotlib.pylab as plt


from utils.plot import plot_errors
from mesh_sizes import mesh_sizes


with open('measurements.pkl') as f:
    measurements = pickle.load(f)

plot_errors(measurements, mesh_sizes)

plt.show()
