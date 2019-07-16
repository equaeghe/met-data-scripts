import netCDF4 as nc  # library to manipulate netCDF4 files
import numpy as np  # library for scientific computing
import convention_read as cr
import matplotlib.pyplot as plt

f = nc.Dataset('FINO1.nc', 'r')  # load our netCDF4 version of FINO1 data-set

# %% Range check

for g in f.groups:  # each group g corresponds to a device
    for v in f[g].groups:  # each variable v measured by that device
        if 'valid_range' in f[g][v].ncattrs():  # is a range defined for q?
            l, u = f[g][v].getncattr('valid_range')  # range bounds
            # go over statistics s of interest available for v
            for s in set(f[g][v].variables) & {'avg', 'min', 'max'}:
                x = f[g][v][s][:]  # array of statistic s values of variable v
                xl = np.nanmin(x)
                xu = np.nanmax(x)
                # if some range bound is violated, provide this as output
                print(g, v,
                      s + ':', '{!s} < {!s}'.format(xl, l) if xl < l else '',
                          ';', '{!s} > {!s}'.format(xu, u) if xu > u else '')

# %% Check whether min <= avg <= max  and 2 * std <= max - min

for g in f.groups:  # each group g corresponds to a device
    for v in f[g].groups:  # each variable v measured by that device
        if 'avg' in f[g][v].variables:
            xa = f[g][v]['avg'][:]
            if 'min' in f[g][v].variables:
                xi = f[g][v]['min'][:]
                a = ~np.isnan(xa) & ~np.isnan(xi)
                b = (xa[a] + 2 ** cr.accuracy(xa[a]) <
                     xi[a] - 2 ** cr.accuracy(xi[a]))
                if any(b):
                    print(g, v, 'min', xa[a][b], xi[a][b])
            if 'max' in f[g][v].variables:
                xm = f[g][v]['max'][:]
                a = ~np.isnan(xa) & ~np.isnan(xm)
                b = (xa[a] - 2 ** cr.accuracy(xa[a]) >
                     xm[a] + 2 ** cr.accuracy(xm[a]))
                if any(b):
                    print(g, v, 'max', xa[a][b], xm[a][b])
            if {'std', 'min', 'max'} in set(f[g][v].variables):
                n = 600 * f[g][v].getncattr('sampling_frequency')
#                c = np.sqrt(n / (n - 1)) / 2
                c = .5  # more conservative than above
                xs = f[g][v]['std'][:]
                a = ~np.isnan(xs) & ~np.isnan(xi) & ~np.isnan(xm)
                b = c * (np.abs(xm[a] - xi[a])
                         + 2 ** cr.accuracy(xm[a]) + 2 ** cr.accuracy(xi[a])
                         < (xs[a] - 2 ** cr.accuracy(xm[a])))
                if any(b):
                    print(g, v, 'std',
                          c * np.abs(xm[a][b] - xi[a][b]), xs[a][b])

# %% Check for outliers

for g in f.groups:  # each group g corresponds to a device
    for v in f[g].groups:  # each variable v measured by that device
        for s in set(f[g][v].variables) - {'flag'}:
            x = f[g][v][s][:]
            if 'valid_range' in f[g][v].ncattrs():  # is a range defined for v?
                l, u = f[g][v].getncattr('valid_range')  # range bounds
                x = x[(x >= l) & (x <= u)]
            μ = np.nanmean(x)
            σ = np.nanstd(x)
            a = ~np.isnan(x)
            y = np.abs(x[a] - μ) / σ
            print(g, v, s, μ, σ, np.mean(y > 3))

# %% Outlier plots

sw = 2.5
sh = 2.5

plt.ioff()

plt.rc('xtick', labelsize=8)  # fontsize of the tick labels
plt.rc('ytick', labelsize=8)  # fontsize of the tick labels

for g in f.groups:  # each group g corresponds to a device
    heights = f[g]['height_' + g][:]
    for v in f[g].groups:  # each variable v measured by that device
        rows = len(heights)
        cols = len(f[g][v].variables) - 1
        fig, axarr = plt.subplots(rows, cols,
                                  squeeze=False, sharex='col', sharey='col',
                                  figsize=((cols + 0.5) * sw,
                                           (rows + 0.5) * sh))
        fig.suptitle("FINO1 {}/{} ".format(g, v)
                     + '(' + f[g][v].long_name
                     + ((", units: " + f[g][v].units)
                        if "units" in f[g][v].ncattrs()
                        else '')
                     + ')', y=1.1)
        nh = len(heights)
        for i, h in enumerate(heights):
            axarr[nh - i - 1, 0].set_ylabel(str(h) + " m")
        for j, s in enumerate(set(f[g][v].variables) - {'flag'}):
            axarr[0, j].set_title(s)
            X = f[g][v][s][:].data
            Z = np.zeros(X.shape)
            dX = np.abs(X[1:, :] - X[:-1, :])
            if v == 'Wr':
                dX = np.where(dX > 180, 360 - dX, dX)
            Z[:-1] += dX / 2
            Z[1:] += dX / 2
            Zm = np.nanmedian(Z)
            Zm = Zm if Zm > 0 else np.nanstd(Z)
            Xl, Xu = np.nanpercentile(X, [.1, 99.9])
            Xn, Xx = np.nanmin(X), np.nanmax(X)
            if v == 'Wr':
                Xl, Xu = 0, 360
                if s == 'std':
                    Xu = np.inf
            elif 'valid_range' in f[g][v].ncattrs():
                l, u = f[g][v].getncattr('valid_range')  # range bounds
                if s == 'std':
                    l, u = 0, (u - l) / 2  # bound on std
                Xl, Xu = np.maximum(Xl, l), np.minimum(Xu, u)
            Zu = np.nanpercentile(Z, 99)
            print(g, v, s, Xl, Xu, Zm, Zu)
            for i, h in enumerate(heights):
                ax = axarr[nh - i - 1, j]
                if 'valid_range' in f[g][v].ncattrs():
                    if not v == 'Wr':
                        ax.axvline(l, color='r', linewidth=3)
                        ax.axvline(u, color='r', linewidth=3)
                x = X[:, i]
                if np.all(np.isnan(x)):
                    continue
                z = Z[:, i]
                pvs = [1.5625, 3.125, 6.25, 12.5, 25,
                       75, 87.5, 93.75, 96.875, 98.4375]
                if v != 'Wr':
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
#                ax.axis('tight')
                ax.set_xlim(Xn, Xx)
                lZu = np.log10(Zu)
                ax.set_yscale('symlog', basey=10,
                              linthreshy=Zu, linscaley=lZu if lZu > 0 else 1)
        fig.savefig("FINO1-{}_{}.pdf".format(g, v),
                    bbox_inches='tight', pad_inches=.5)
        plt.close(fig)
