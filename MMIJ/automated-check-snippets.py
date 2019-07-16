import netCDF4 as nc  # library to manipulate netCDF4 files
import numpy as np  # library for scientific computing
import matplotlib.pyplot as plt


# load our netCDF4 version of MMIJ data-set
f = nc.Dataset('MMIJ-compound-sem.nc', 'r')

# %% Range check

for g in f.groups:  # each group g corresponds to a field (air, wind, prec)
    for d in f[g].groups:  # each device d
        for v in f[g][d].variables:  # each variable v measured by that device
            if 'valid_range' in f[g][d][v].ncattrs():
                l, u = f[g][d][v].getncattr('valid_range')  # range bounds
                # go over statistics s of interest available for v
                for s in {'avg', 'min', 'max'}:
                    x = f[g][d][v][:][s]
                    xl = np.nanmin(x)
                    xu = np.nanmax(x)
                    # if some range bound is violated, provide this as output
                    print(g, d, v, s + ':',
                          '{!s} < {!s}'.format(xl, l) if xl < l else '', ';',
                          '{!s} > {!s}'.format(xu, u) if xu > u else '')

# %% Check whether min <= avg <= max and 2 * std <= max - min

for g in f.groups:  # each group g corresponds to a field (air, wind, prec)
    for d in f[g].groups:  # each device d
        for v in f[g][d].variables:  # each variable v measured by that device
            if not hasattr(f[g][d][v].dtype, 'names'):
                continue
            if f[g][d][v].dtype.names is None:
                continue
            xa = f[g][d][v][:]['avg']
            xi = f[g][d][v][:]['min']
            a = ~np.isnan(xa) & ~np.isnan(xi)
            b = xa[a] < xi[a]
            if any(b):
                print(g, d, v, 'min', np.max(xi[a][b] - xa[a][b]))
            xm = f[g][d][v][:]['max']
            a = ~np.isnan(xa) & ~np.isnan(xm)
            b = xa[a] > xm[a]
            if any(b):
                print(g, d, v, 'max', np.max(xa[a][b] - xm[a][b]))
            xs = f[g][d][v][:]['std']
            a = ~np.isnan(xs) & ~np.isnan(xi) & ~np.isnan(xm)
            b = np.abs(xm[a] - xi[a]) < 2 * xs[a]
            if any(b):
                print(g, d, v, 'std',
                      np.max(2 * xs[a] - np.abs(xm[a] - xi[a])))

# %% Check for outliers

for g in f.groups:  # each group g corresponds to a device
    for d in f[g].groups:  # each device d
        for v in f[g][d].variables:  # each variable v measured by that device
            if not hasattr(f[g][d][v].dtype, 'names'):
                continue
            if f[g][d][v].dtype.names is None:
                continue
            for s in {'min', 'max', 'avg', 'std'}:
                x = f[g][d][v][:][s]
                if 'valid_range' in f[g][d][v].ncattrs():
                    l, u = f[g][d][v].getncattr('valid_range')  # range bounds
                    x = x[(x >= l) & (x <= u)]
                a = np.isnan(x)
                if a.all():
                    continue
                μ = np.nanmean(x)
                σ = np.nanstd(x)
                y = np.abs(x[~a] - μ) / σ
                print(g, d, v, s, μ, σ, np.mean(y > 3))

# %% Outlier plots

sw = 2.5
sh = 2.5

plt.ioff()

plt.rc('xtick', labelsize=8)  # fontsize of the tick labels
plt.rc('ytick', labelsize=8)  # fontsize of the tick labels

for g in f.groups:
    if g == 'air':
        places = f[g]['height_' + g][:]
    for d in f[g].groups:  # each device d
        if g == 'wind':
            heights = f[g][d]['height_' + d][:]
        elif g == 'prec':
            places = f[g][d]['location_' + d][:]
        for v in f[g][d].variables:  # each variable v measured by that device
            if not hasattr(f[g][d][v].dtype, 'names'):
                continue
            if f[g][d][v].dtype.names is None:
                continue
            if g == 'wind':
                if v.startswith('True'):
                    places = heights
                else:
                    Θ = (f[g][d]['direction_' + d][:] if d == 'CA' else
                         f['direction'][:])
                    places = [((ih, iθ), "{} m, {}°".format(h, θ))
                              for iθ, θ in enumerate(Θ)
                              for ih, h in enumerate(heights)
                              if (h, θ) not in {(27., 226.5), (27., 346.5),
                                                (58.5, 226.5), (58.5, 346.5),
                                                (92., 46.5), (92., 166.5),
                                                (92., 286.5)}]
            rows = len(places)
            fig, axarr = plt.subplots(rows, 4, squeeze=False,
                                      sharex='col', sharey='col',
                                      figsize=((4 + 0.5) * sw,
                                               (rows + 0.5) * sh))
            fig.suptitle("MMIJ {}/{} ".format(d, v)
                         + '(' + f[g][d][v].long_name
                         + ((", units: " + f[g][d][v].units)
                            if "units" in f[g][d][v].ncattrs()
                            else '')
                         + ')', y=1.1)
            nh = len(places)
            for i, h in enumerate(places):
                if (g != 'wind') or v.startswith('True'):
                    hlabel = str(h) + ' m'
                else:
                    hlabel = h[1]
                axarr[nh - i - 1, 0].set_ylabel(hlabel)
            for j, s in enumerate(['min', 'max', 'avg', 'std']):
                axarr[0, j].set_title(s)
                X = f[g][d][v][:][s]
                Z = np.zeros(X.shape)
                dX = np.abs(X[1:, :] - X[:-1, :])
                if v.endswith('Wd') and s == 'avg':
                    dX = np.where(dX > 180, 360 - dX, dX)
                Z[:-1] += dX / 2
                Z[1:] += dX / 2
                Zm = np.nanmedian(Z)
                Zm = Zm if Zm > 0 else np.nanstd(Z)
                Xl, Xu = np.nanpercentile(X, [.1, 99.9])
                Xn, Xx = np.nanmin(X), np.nanmax(X)
                if v.endswith('Wd'):
                    Xl, Xu = 0, 360
                    if s == 'std':
                        Xu = np.inf
                elif 'valid_range' in f[g][d][v].ncattrs():
                    l, u = f[g][d][v].getncattr('valid_range')  # range bounds
                    if s == 'std':
                        l, u = 0, (u - l) / 2  # bound on std
                    Xl, Xu = np.maximum(Xl, l), np.minimum(Xu, u)
                Zu = np.nanpercentile(Z, 99)
                Zu = Zu if Zu > 0 else 3 * Zm
                print(g, d, v, s, Xl, Xu, Zm, Zu)
                for i, h in enumerate(places):
                    ax = axarr[nh - i - 1, j]
                    if 'valid_range' in f[g][d][v].ncattrs():
                        if not v.endswith('Wd'):
                            ax.axvline(l, color='r', linewidth=3)
                            ax.axvline(u, color='r', linewidth=3)
                    if (g == 'wind') and not v.startswith('True'):
                        x = X[:, h[0][0], h[0][1]]
                        z = Z[:, h[0][0], h[0][1]]
                    elif d == 'PM':
                        x = X[:, 0, i]
                        z = Z[:, 0, i]
                    else:
                        x = X[:, i]
                        z = Z[:, i]
                    if np.all(np.isnan(x)):
                        continue
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
                    ax.axis('tight')
                    ax.set_xlim(Xn, Xx)
                    lZu = np.log10(Zu)
                    ax.set_yscale('symlog', basey=10,
                                  linthreshy=Zu,
                                  linscaley=lZu if lZu > 0 else 1)
            fig.savefig("MMIJ-{}_{}.pdf".format(d, v),
                        bbox_inches='tight', pad_inches=.5)
            plt.close(fig)
