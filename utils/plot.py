"""Plotting utilities."""

import numpy as np
import matplotlib.pylab as plt

def unpack(errs):
    pop = []
    growth = []
    for err in errs:
        pop.append(err[0])
        growth.append(err[1])

    return np.asarray(pop), np.asarray(growth)

linestyles = {
    'MidPoint':            '-ok',
    'Zuidema':             '-sk',
    'ClenshawCurtis(9)':   '-.db',
    'GaussLegendre(9)':    '-.^r',
    'AdjGaussLegendre(9)': '-.vm',
    'GENClarkGL(9)':       '--<g',
    'INTClarkGL(9)':       '-->c',

    # 'ClenshawCurtis':      ':ob',
    # 'ClenshawCurtis(3)':   '-.^b',
    # 'ClenshawCurtis(5)':   '--vb',
    # 'GaussLegendre':       ':or',
    # 'GaussLegendre(3)':    '-.^r',
    # 'GaussLegendre(5)':    '--vr',
    # 'AdjGaussLegendre(3)': '-.^m',
    # 'AdjGaussLegendre(5)': '--vm',
    # 'GENClarkQP':          '-og',
    # 'GENClarkGL(3)':       '-.^g',
    # 'GENClarkGL(5)':       '--vg',
    # 'INTClarkQP':          '-oc',
    # 'INTClarkGL(3)':       '-.^c',
    # 'INTClarkGL(5)':       '--vc',
    }


def plot_projections(kernel, method, projections, populations):
    """Plot *projections* and *populations* vs time."""

    N = projections.shape[1]
    T = kernel.T

    # plot projections
    plt.figure()
    plt.title('Projections: %s: %s (%d)'
              % (kernel.name, method.name, N))
    for j, t in enumerate(T):
        plt.plot(kernel.x, projections[j], '-r')
    plt.xlabel('x')
    plt.ylabel('n')

    # plot population vs time
    plt.figure()
    plt.title('Population: %s: %s (%d)'
              % (kernel.name, method.name, N))
    plt.plot(T, populations, 'o-k')
    plt.xlabel('time')
    plt.ylabel('population')



def plot_errors(errors, mesh_sizes, linestyles=linestyles):
    """Plot erros and costs for all kernels and methods in *errors*."""

    for kernel in errors:

        # error vs computational cost
        fig, ax = plt.subplots(1, 2, sharey=True, figsize=(7, 6))
        fig.subplots_adjust(top=0.75)

        for method in sorted(errors[kernel]):
            errs, cnts, tims = errors[kernel][method]
            pop, growth = unpack(errs)
            print kernel, method
            print '  ', pop
            print '  ', growth
            ax[0].loglog(cnts, abs(pop),
                         linestyles.get(method, '-k'), label=method)
            ax[1].loglog(cnts, abs(growth),
                         linestyles.get(method, '-k'), label=method)

        ax[0].set_title('population error')
        ax[1].set_title('growth rate error')
        ax[0].set_ylabel('absolute error')
        for a in ax:
          a.set_xlabel('no. of kernel evals.')

        plt.ylim(1e-16, 1e1)
        plt.legend(bbox_to_anchor=(0.5, 0.8), loc=8, ncol=3, borderaxespad=0.0,
                   bbox_transform=fig.transFigure)

        plt.savefig('plots/ErrorCost_%s.png' % kernel)

        # error vs run time
        fig, ax = plt.subplots(1, 2, sharey=True, figsize=(7, 6))
        fig.subplots_adjust(top=0.75)

        for method in sorted(errors[kernel]):
            errs, cnts, tims = errors[kernel][method]
            pop, growth = unpack(errs)
            tims = np.asarray(tims)
            tims[ tims <= 0.0 ] = 1e-4
            # if not np.any(tims <= 0.0):
            ax[0].loglog(np.asarray(tims), abs(pop),
                         linestyles.get(method, '-k'), label=method)
            ax[1].loglog(np.asarray(tims), abs(growth),
                         linestyles.get(method, '-k'), label=method)

        ax[0].set_title('population error')
        ax[1].set_title('growth rate error')
        ax[0].set_ylabel('absolute error')
        for a in ax:
          a.set_xlabel('run time (s)')

        plt.ylim(1e-16, 1e1)
        plt.legend(bbox_to_anchor=(0.5, 0.8), loc=8, ncol=3, borderaxespad=0.0,
                   bbox_transform=fig.transFigure)

        plt.savefig('plots/ErrorTime_%s.png' % kernel)
