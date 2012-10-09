"""Plotting utilities."""

import numpy as np
import matplotlib.pylab as plt

def unpack(errs):
    pop = []
    growth = []
    for err in errs:
        pop.append(err[0])
        growth.append(err[1])

    return pop, growth

linestyles = {
    'MidPoint':            '-ok',
    'Zuidema':             '-sk',
    'ClenshawCurtis':      ':ob',
    'ClenshawCurtis(3)':   '-.^b',
    'ClenshawCurtis(5)':   '--vb',
    'ClenshawCurtis(9)':   '-sb',
    'GaussLegendre':       ':or',
    'GaussLegendre(3)':    '-.^r',
    'GaussLegendre(5)':    '--vr',
    'GaussLegendre(9)':    '-sr',
    'AdjGaussLegendre(3)': '-.^m',
    'AdjGaussLegendre(5)': '--vm',
    'AdjGaussLegendre(9)': '-sm',
    'GENClarkQP':          '-og',
    'GENClarkGL(3)':       '-.^g',
    'GENClarkGL(5)':       '--vg',
    'GENClarkGL(9)':       '-sg',
    'INTClarkQP':          '-oc',
    'INTClarkGL(3)':       '-.^c',
    'INTClarkGL(5)':       '--vc',
    'INTClarkGL(9)':       '-sc',
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

        # error vs mesh size
        plt.figure(figsize=(7, 6))
        plt.axes([0.125, 0.1, 0.8, 0.55])
        # plt.title('Error vs mesh size: ' + str(kernel))
        for method in sorted(errors[kernel]):
            if method.find('ClarkCC') > 0:
                continue

            errs, cnts, tims = errors[kernel][method]
            pop, gwroth = unpack(errs)
            plt.loglog(mesh_sizes[kernel, method], 
                       pop, linestyles.get(method, '-k'), label=method)
        plt.xlabel('N - mesh size')
        plt.ylabel('absolute error')
        plt.legend(bbox_to_anchor=(0.5, 1.1), loc=8,
                   ncol=3, borderaxespad=0.)
        plt.savefig('plots/ErrorMesh_%s.png' % kernel)

        # error vs computational cost
        plt.figure(figsize=(7, 6))
        plt.axes([0.125, 0.1, 0.8, 0.55])
        # plt.title('Error vs computational cost: ' + str(kernel))
        for method in sorted(errors[kernel]):
            if method.find('ClarkCC') > 0:
                continue

            errs, cnts, tims = errors[kernel][method]
            pop, gwroth = unpack(errs)
            plt.loglog(cnts, pop, linestyles.get(method, '-k'), label=method)
        plt.xlabel('cost (number of kernel evaluations)')
        plt.ylabel('absolute error')
        plt.legend(bbox_to_anchor=(0.5, 1.1), loc=8,
                   ncol=3, borderaxespad=0.)
        plt.savefig('plots/ErrorCost_%s.png' % kernel)


        # error vs run time
        plt.figure(figsize=(7, 6))
        plt.axes([0.125, 0.1, 0.8, 0.55])

        for method in sorted(errors[kernel]):
            if method.find('ClarkCC') > 0:
                continue

            errs, cnts, tims = errors[kernel][method]
            pop, gwroth = unpack(errs)
            tims = np.asarray(tims)
            if not np.any(tims <= 0.0):
                plt.loglog(np.asarray(tims), pop, linestyles.get(method, '-k'), label=method)

        plt.xlabel('run time (s)')
        plt.ylabel('absolute error')
        plt.legend(bbox_to_anchor=(0.5, 1.1), loc=8,
                   ncol=3, borderaxespad=0.)
        plt.savefig('plots/ErrorTime_%s.png' % kernel)
