import h5py as h5  # library to manipulate HDF5 files
import numpy as np  # library for scientific computing
import convention_read as cr
import matplotlib.pyplot as plt

f = h5.File('OWEZ.h5', 'r')
fd = f['devices']

N = len(f['time/timestamp'])

# %% Range check

for d in fd:  # each d corresponds to a device
    for v in fd[d]:  # each variable v measured by that device
        if 'range' in fd[d][v].attrs:  # is a range defined for v?
            lb, ub = fd[d][v].attrs['range']  # range bounds
            for l in fd[d][v]:  # locations l where such a device is present
                # go over statistics s of interest available for v
                for s in {'Mean', 'Min', 'Max'}:
                    x = fd[d][v][l][s][:]  # array of stat. s values of var. v
                    xlb = np.nanmin(x)
                    xub = np.nanmax(x)
                    # if some range bound is violated, provide this as output
                    print(d, v, l, s + ':',
                          '{!s} < {!s}'.format(xlb, lb) if xlb < lb else '',
                          ';',
                          '{!s} > {!s}'.format(xub, ub) if xub > ub else '')

# %% Check whether min <= avg <= max

for g in fd:  # each d corresponds to a device
    for v in fd[d]:  # each variable v measured by that device
        for l in fd[d][v]:  # each location l where such a device is present
            xa = fd[d][v][l]['Mean'][:]
            xi = fd[d][v][l]['Min'][:]
            a = ~np.isnan(xa) & ~np.isnan(xi)
            b = (xa[a] + 2 ** cr.accuracy(xa[a]) <
                 xi[a] - 2 ** cr.accuracy(xi[a]))
            if any(b):
                print(d, v, l, 'Min', np.max(xi[a][b] - xa[a][b]))
            xm = fd[d][v][l]['Max'][:]
            a = ~np.isnan(xa) & ~np.isnan(xm)
            b = (xa[a] - 2 ** cr.accuracy(xa[a]) >
                 xm[a] + 2 ** cr.accuracy(xm[a]))
            if any(b):
                print(d, v, l, 'Max', np.max(xa[a][b] - xm[a][b]))

# %% Check for outliers

for d in fd:  # each d corresponds to a device
    for v in fd[d]:  # each variable v measured by that device
        for l in fd[d][v]:  # each location l where such a device is present
            for s in ['Min', 'Max', 'Mean', 'StdDev']:
                x = x[~np.isnan(x)]
                x = np.float32(fd[d][v][l][s][:])
                if 'range' in fd[d][v].attrs:  # is a range defined for v?
                    lb, ub = fd[d][v].attrs['range']  # range bounds
                    x = x[(x >= lb) & (x <= ub)]
                μ = np.mean(x)
                σ = np.std(x)
                y = np.abs(x - μ) / σ
                print(d, v, l, s, μ, σ, np.mean(y > 3))

# %% Outlier plots

sw = 2.5
sh = 2.5

plt.ioff()

plt.rc('xtick', labelsize=8)  # fontsize of the tick labels
plt.rc('ytick', labelsize=8)  # fontsize of the tick labels

for d in fd:  # each d corresponds to a device
    for v in fd[d]:  # each variable v measured by that device
        n = len(fd[d][v])
        fig, axarr = plt.subplots(n, 4,
                                  squeeze=False, sharex='col', sharey='col',
                                  figsize=((4 + 0.5) * sw,
                                           (n + 0.5) * sh))
        fig.suptitle("OWEZ {}/{} ".format(d, v)
                     + '(' + fd[d][v].attrs['description']
                     + ((", units: " + fd[d][v].attrs['unit'])
                        if "unit" in fd[d][v].attrs
                        else '')
                     + ')', y=1.1)
        for i, l in enumerate(fd[d][v]):
            axarr[n - i - 1, 0].set_ylabel(l)
        for j, s in enumerate(['Min', 'Max', 'Mean', 'StdDev']):
            axarr[0, j].set_title(s)
            X = np.nan * np.ones((N, n), 'f4')
            for i, l in enumerate(fd[d][v]):
                X[:, i] = fd[d][v][l][s][:]
            Z = np.zeros(X.shape)
            dX = np.abs(X[1:, :] - X[:-1, :])
            if v in {'wd', 'cd7', 'cd11', 'vd'}:
                dX = np.where(dX > 180, 360 - dX, dX)
            Z[:-1] += dX / 2
            Z[1:] += dX / 2
            Zm = np.nanmedian(Z)
            Zm = Zm if Zm > 0 else np.nanstd(Z)
            Xl, Xu = np.nanpercentile(X, [.1, 99.9])
            Xn, Xx = np.nanmin(X), np.nanmax(X)
            if v in {'wd', 'cd7', 'cd11', 'vd'}:
                Xl, Xu = 0, 360
                if s == 'StdDev':
                    Xu = np.inf
            elif 'range' in fd[d][v].attrs:
                lb, ub = fd[d][v].attrs['range']  # range bounds
                if s == 'StdDev':
                    lb, ub = 0, (ub - lb) / 2  # bound on std
                Xl, Xu = np.maximum(Xl, lb), np.minimum(Xu, ub)
            Zu = np.nanpercentile(Z, 99)
            print(d, v, l, s, Xl, Xu, Zu)
            for i, l in enumerate(fd[d][v]):
                ax = axarr[n - i - 1, j]
                if 'range' in fd[d][v].attrs:
                    if v not in {'wd', 'cd7', 'cd11', 'vd'}:
                        ax.axvline(lb, color='r', linewidth=3)
                        ax.axvline(ub, color='r', linewidth=3)
                x = X[:, i]
                if np.all(np.isnan(x)):
                    continue
                z = Z[:, i]
                pvs = [1.5625, 3.125, 6.25, 12.5, 25,
                       75, 87.5, 93.75, 96.875, 98.4375]
                if v not in {'wd', 'cd7', 'cd11', 'vd'}:
                    ax.axvline(np.nanmedian(x), color='b',
                               linewidth=2, linestyle='dashed')
                    for p in np.nanpercentile(x, pvs):
                        ax.axvline(p, color='b',
                                   linewidth=1, linestyle='dotted')
                ax.axhline(np.nanmedian(z), color='b',
                           linewidth=2, linestyle='dashed')
                for p in np.nanpercentile(z, pvs):
                    ax.axhline(p, color='b',
                               linewidth=1, linestyle='dotted')
                a = (x < Xl) | (x > Xu) | (z > Zu)
                ax.plot(np.where(a, x, np.nan), np.where(a, z, np.nan),
                        color='k', alpha=.25, linestyle='-', marker='.')
                ax.axis('tight')
                ax.set_xlim(Xn, Xx)
                lZu = np.log10(Zu)
                ax.set_yscale('symlog', basey=10,
                              linthreshy=Zu, linscaley=lZu if lZu > 0 else 1)
        fig.savefig("OWEZ-{}_{}.pdf".format(d, v),
                    bbox_inches='tight', pad_inches=.5)
        plt.close(fig)
