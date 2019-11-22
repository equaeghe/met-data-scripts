import netCDF4 as nc
import numpy as np


f = nc.Dataset('MMIJ.nc', 'r')

for c in f.groups:
    for i in f[c].groups:
        for v in f[c][i].groups:
            ε_a = ε_r = 0
            if 'uncertainty_abs' in f[c][i][v].ncattrs():
                ε_a = f[c][i][v].getncattr('uncertainty_abs')
            if 'uncertainty_rel' in f[c][i][v].ncattrs():
                ε_r = f[c][i][v].getncattr('uncertainty_rel')
            if (ε_a == 0) and (ε_r == 0):
                continue
            print(c, i, v)
            x_min = f[c][i][v]['min'][:]
            x_max = f[c][i][v]['max'][:]
            x_avg = f[c][i][v]['avg'][:]
            x_std = f[c][i][v]['std'][:]
            n = 2400
            δ = (x_max - x_min) / 2
            τ2 = (δ/n)**2
            σ_min2 = ε_a**2 + ε_r**2 * x_min**2
            σ_max2 = ε_a**2 + ε_r**2 * x_max**2
            σ_avg2 = (ε_a**2 + ε_r**2 * (x_avg**2 + x_std**2)) / n
            σ_std2 = (ε_a**2 + ε_r**2 * (x_avg**2 + 3 * x_std**2)) / n
            ε_min = np.sqrt(σ_min2 + τ2)
            ε_max = np.sqrt(σ_max2 + τ2)
            ε_avg = np.sqrt(σ_avg2 + τ2)
            ε_std = np.sqrt(σ_std2 + τ2)
            x_std_unbiased = np.sqrt(
                np.maximum(x_std**2 - (ε_a**2 + ε_r**2 * x_avg**2), 0))
            print(
                np.round(100 * np.nanmean(ε_min / np.abs(x_min)), 3),
                np.round(100 * np.nanmean(ε_max / np.abs(x_max)), 4),
                np.round(100 * np.nanmean(ε_avg / np.abs(x_avg)), 5),
                np.round(100 * np.nanmean(ε_std / x_std), 2),
                np.round(100 * (1 - np.nanmean(x_std_unbiased / x_std)), 1)
            )
