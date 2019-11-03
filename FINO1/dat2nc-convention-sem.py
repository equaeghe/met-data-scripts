import zipfile as zf
import pandas as pd
import netCDF4 as nc
import numpy as np
from ruamel.yaml import YAML as yaml


def yaml_load(f):
    return yaml(typ='safe').load(f)


# %%

# The input file must be a zip file with dat files that can be obtained via
# http://fino.bsh.de/
input_file = "FINO1_20040101_20161231-with_fixes.zip"

# The output file is a netCDF4-file
output_file = "FINO1-compound.nc"


dsd = {}
with zf.ZipFile(input_file) as finoz:
    path, *files = finoz.namelist()
    for f in [f for f in files if f.endswith('.dat')]:
        key = f.split('/')[1].split('.')[0] # foo/key.bar
        dsd[key] = pd.read_csv(finoz.open(f), sep='\t', header=4, skiprows=[5],
                               index_col=0, parse_dates=True,
                               na_values=['-999.99', '-999.9', '-999'])

times = pd.date_range('2004-01-01 00:10:00', '2017-01-01 00:00:00',
                      freq='10min')

# %%

f = nc.Dataset(output_file, 'w')
with open('metadata-global.yaml') as g:
    for attribute, value in yaml_load(g).items():
        f.setncattr_string(attribute, value)
# We override the comment from 'metadata-global.yaml', as it describes also how
# the datasets have been structured within the created file. This differs for
# this file from the finally chosen solution. Also further below we make less
# use of the metadata files, but provide custom metadata.
f.comment = ("The statistics datasets have as values a compound data structure "
             "with as possible components the minimum ('min'), "
             "maximum ('max'), average ('avg'), and standard deviation ('std') "
             "of the samples within the given time or a value ('val') "
             "characterizing the measurement for the time interval. The "
             "estimated accuracy (error) for the statistics is\n"
             "\t* the instrument's provided error for 'min' and 'max',\n"
             "\t* the maximum of the standard or propagated error for 'avg' "
             "and 'std'.\n\n"
             "The 'min', 'max', 'avg', and 'std' statistics are recorded as a "
             "floating point numbers for which the accuracy is encoded using "
             "the following convention: if x is the number and e its accuracy, "
             "then we store y = (2 * k ± 1) * d, where\n"
             "\t* d = 2 ** floor(log2(e)) / 2 and\n"
             "\t* k = round(((x / d) ∓ 1) / 2).\n"
             "So y is the odd multiple of d nearest to x, implying that the "
             "possible values for y are spaced 2 * d apart. This step size "
             "2 * d is the largest power of two not larger than e, so that "
             "|x - y| <= e / 2. Because y is an odd multiple of d, the "
             "denominator of y seen as an irreducable fraction is 1 / d, so "
             "that d can effectively be determined from y. To recover d from "
             "y, let n = np.abs(np.float32(y)).view(np.int32) be its bit-level "
             "integer representation, then b = (n >> 23) - 127 is the binary "
             "exponent and m = n & np.int32(2**23 - 1) the 23-bit significand. "
             "So then d = (m & -m) * 2**-23 * 2**b.\n\n" # https://stackoverflow.com/questions/46093123, https://stackoverflow.com/questions/18806481
             "There is also always a 'flag' component in the compound data "
             "structure, providing information about quality and missingness. "
             "The flag values are encoded using the following "
             "enumeration:\n"
             "\t* 'missing': -1,\n"
             "\t* 'raw': 0,\n"
             "\t* 'doubtful quality': 1,\n"
             "\t* 'quality controlled': 2."
             "\n\n"
             "Custom dataset attributes:\n"
             "\t* 'uncertainty_abs' describes the absolute accuracy of the "
             "sampled values.\n"
             "\t* 'uncertainty_rel' describes the relative accuracy of the "
             "sampled values as a number between 0 and 1.\n"
             "\t* 'sampling_frequency' describes the frequency in Hz with "
             "which the quantity has been sampled, which determines the number "
             "n of samples.\n"
             "\t* 'standard_error_avg' describes the standard error of the avg "
             "as a fraction 1/sqrt(n) of the std values.\n"
             "\t* 'standard_error_std' describes the standard error of the std "
             "as a fraction 1/sqrt(2·(n-1)) of the std values.\n"
             "\t* 'uncertainty_propagated_avg' describes the accuracy of the "
             "avg values as a fraction 1/sqrt(n) of the sampled values' "
             "(absolute) accuracy, as derived from 'uncertainty_abs' and "
             "'uncertainty_rel'.\n"
             "\t* 'uncertainty_propagated_std' describes the accuracy of the "
             "std values as a fraction 1/sqrt(n-1) of the sampled values' "
             "(absolute) accuracy, as derived from 'uncertainty_abs' and 'uncertainty_rel'."
             "\n\n"
             "Custom instrument attributes:\n"
             "\t* 'manufacturer' describes the manufacturer of the instrument, "
             "typically by giving its name.\n"
             "\t* 'part_number' describes the manufacturer's part number for "
             "the instrument.")

# %%

flags = {"missing": -1,
         "raw": 0, "doubtful quality": 1, "quality controlled": 2}
flag = f.createEnumType('i4', "flag", flags) # i4 instead of i1 because it will
                                             # be included in compound type and
                                             # padded up to 4 bytes anyway; the
                                             # padding seems to be random
                                             # instead of constant, resulting
                                             # in inefficient compression
allstats_np = np.dtype({'names': ["min", "max", "avg", "std", "flag"],
                        'formats': ['f4', 'f4', 'f4', 'f4', flag]},
                       align=False)
allstats = f.createCompoundType(allstats_np, "all_statistics")
maxavgstd_np = np.dtype({'names': ["max", "avg", "std", "flag"],
                         'formats': ['f4', 'f4', 'f4', flag]},
                        align=False)
maxavgstd = f.createCompoundType(maxavgstd_np, "maxavgstd_statistics")
avgstd_np = np.dtype({'names': ["avg", "std", "flag"],
                      'formats': ['f4', 'f4', flag]},
                     align=False)
avgstd = f.createCompoundType(avgstd_np, "avgstd_statistics")
avg_np = np.dtype({'names': ["avg", "flag"], 'formats': ['f4', flag]},
                  align=False)
avg = f.createCompoundType(avg_np, "avg_statistic")
value_np = np.dtype({'names': ["val", "flag"], 'formats': ['f4', flag]},
                    align=False)
value = f.createCompoundType(value_np, "value")

# %%

## GROUP HIERARCHY ##

# INSTRUMENT GROUPS
instruments_locations = {}
with open('metadata-instruments.yaml') as g:
    for instrument, metadata in yaml_load(g).items():
        instruments_locations[instrument] = metadata['locations']
        del metadata['locations']
        f.createGroup(instrument)
        f[instrument].setncatts(metadata)

# %%

## DIMENSIONS ##

# GLOBAL
# auxiliary
f.createDimension("bounds_edges", 2)
# time
f.createDimension("time", len(times))
f.createVariable('time', 'u4', ("time",), zlib=True, complevel=9,
                 fill_value=False)
f["time"][:] = nc.date2num(times.to_pydatetime(), "seconds since 1970-01-01")
f["time"].units = "seconds since 1970-01-01"
f["time"].delta_t = "0000-00-00 00:10:00"
f["time"].actual_range = np.array([times[0].timestamp(),
                                   times[-1].timestamp()], dtype='u4')
f["time"].axis = 'T'
f["time"].bounds = 'time_bounds'
# time bounds (of statistics, denote open-closed intervals(?))
f.createVariable('time_bounds', 'u4', ("time", "bounds_edges"),
                 zlib=True, complevel=9, fill_value=False)
f["time_bounds"][:, 0] = f["time"][:] - 600
f["time_bounds"][:, 1] = f["time"][:]
f["time_bounds"].units = "seconds since 1970-01-01"
f["time_bounds"].long_name = ("bounds of (closed-open) time intervals "
                              "over which the statistics have been calculated")
# height/level
f.createDimension("level", 9)
f.createVariable("level", 'u1', ("level",), fill_value=False)
f["level"][:] = np.array(range(9), dtype='u1')
f.createVariable("height", 'f4', ("level",), fill_value=False)
f["height"][:] = np.array([21, 33.5, 41, 51, 61, 71, 81, 91, 101], dtype='f4')
f["height"].setncatts(
    {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
     'long_name': "boom height above LAT (Lowest Astronomical Tide)"})

# %%

# PER-INSTRUMENT
for instrument, locations in instruments_locations.items():
    # levels
    levels = "level_" + instrument
    f[instrument].createDimension(levels, len(locations['levels']))
    f[instrument].createVariable(levels, 'u1', (levels), fill_value=False)
    f[instrument][levels].compress = "level"
    f[instrument][levels][:] = np.array(locations['levels'], 'u1')
    # measurement height
    height = "height_" + instrument
    f[instrument].createVariable(height, 'u1', (levels), fill_value=False)
    f[instrument][height][:] = np.array(locations['measurement_height'], 'f4')
    f[instrument][height].setncatts(
        {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
         'long_name': "measurement height above LAT (Lowest Astronomical Tide)"})
    # position
    position = "position_" + instrument
    f[instrument].createVariable(position, np.str, (levels), fill_value=False)
    f[instrument][position][:] = np.array(locations['position'], np.str)
    # orientation
    if locations['orientation'] is not None:
        orientation = "orientation_" + instrument
        f[instrument].createVariable(orientation, 'f4', (levels), fill_value=False)
        f[instrument][orientation][:] = np.array(locations['orientation'], 'f4')
        f[instrument][orientation].setncatts(
            {'units': "degree", 'long_name': "instrument orientation on mast"})
    # mast distance
    if locations['mast_distance'] is not None:
        distance = "distance_" + instrument
        f[instrument].createVariable(distance, 'f4', (levels), fill_value=False)
        f[instrument][distance][:] = np.array(locations['mast_distance'], 'f4')
        f[instrument][orientation].setncatts(
            {'units': "m",
             'long_name': "horizontal distance to outer edge of mast"})

# %%

## SIGNALS ##
level2fnheight = {i: str(h) + "m"
                  for i, h in enumerate([20, 33, 40, 50, 60, 70, 80, 90, 100])}
with open('metadata-quantities.yaml') as g:
    quantities_metadata = yaml_load(g)
extra_uncertainty_abs = {'PS': {'Ni'}, 'PYR': {'Gs'}}

# %%
dtypes = {}
for instrument, quantities in quantities_metadata.items():
    dtypes[instrument] = {}
    for name, properties in quantities.items():
        stats = properties['statistics']
        del properties['statistics']
        if stats == ['min', 'max', 'avg', 'std']:
            dtype = allstats
        elif stats == ['max', 'avg', 'std']:
            dtype = maxavgstd
        elif stats == ['avg', 'std']:
            dtype = avgstd
        elif stats == ['avg']:
            dtype = avg
        elif stats == ['val']:
            dtype = value
        else:
            raise ValueError("No dtype known for ‘{}’.".format(stats))
        dtypes[instrument][name] = dtype
        s = f[instrument].createVariable(name, dtype,
                                         ("time", "level_" + instrument),
                                         zlib=True, complevel=9,
                                         fill_value=False)
        s.setncatts(properties)
        if instrument in extra_uncertainty_abs:
            if name in extra_uncertainty_abs[instrument]:
                s.setncattr('uncertainty_abs', np.float32(0.005))
                # based on resolution in CSV
        # Fill cell_methods property as needed
        if dtype is not value:
            T = np.float32(600.) # number of seconds in 10 minutes
            ν = s.sampling_frequency
            n = int(ν * T)
            cell_method_details = (
                "(interval: 10 minutes comment: "
                "sampled at {} Hz, so {} samples)").format(ν, n)
            s.cell_methods = ""
            if dtype is allstats:
                s.cell_methods += (
                    "['min'] time: minimum " + cell_method_details + '\n')
            if dtype in {allstats, maxavgstd}:
                s.cell_methods += (
                    "['max'] time: maximum " + cell_method_details + '\n')
            if dtype in {allstats, maxavgstd, avgstd, avg}:
                s.cell_methods += (
                    "['avg'] time: mean " + cell_method_details)
                s.uncertainty_propagated_avg = np.float32(1/np.sqrt(n))
            if dtype in {allstats, maxavgstd, avgstd}:
                s.cell_methods += ('\n' +
                    "['std'] time: standard_deviation " + cell_method_details)
                s.standard_error_avg = np.float32(1/np.sqrt(n))
                s.standard_error_std = np.float32(1/np.sqrt(2 * (n-1)))
                s.uncertainty_propagated_std = np.float32(1/np.sqrt(n-1))
        # Fill variables with default values
        if dtype is value:
            # pre-fill values with NaN and flags with 'missing'
            s[:] = (np.nan, flags["missing"])
            continue
        if dtype is avg:
            # pre-fill values with NaN, accuracies with largest int8,
            # and flags with 'missing'
            s[:] = (np.nan, flags["missing"])
            continue
        if dtype is avgstd:
            # pre-fill values with NaN, accuracies with largest int8,
            # and flags with 'missing'
            s[:] = (np.nan, np.nan, flags["missing"])
            continue
        if dtype is maxavgstd:
            # pre-fill values with NaN, accuracies with largest int8,
            # and flags with 'missing'
            s[:] = (np.nan, np.nan, np.nan, flags["missing"])
            continue
        if dtype is allstats:
            # pre-fill values with NaN, accuracies with largest int8,
            # and flags with 'missing'
            s[:] = (np.nan, np.nan, np.nan, np.nan, flags["missing"])
            continue

# %%

def accround(x, ε):
    y = x
    b = ε != 0.0
    δ = 2 ** np.floor(np.log2(ε[b])) / 2
    k = np.round(((x[b] / δ) - 1) / 2)
    # TODO: δ and k must always be feasible in the sense that the convention holds
    y[b] = (2 * k + 1) * δ
    return y # TODO: case ε == 0.0: if x != 0.0 set last significand bit to 1
             # or is this irrelevant due to limited resolution for all quantities?

def error_abs(x, a, r):
    return np.fmax(a, np.abs(r * x))

for instrument, quantities in quantities_metadata.items():
    ls = instruments_locations[instrument]['levels'] # levels where instrument can be found
    for name, properties in quantities.items():
        dtype = dtypes[instrument][name]
        s = f[instrument][name]
        a = s.uncertainty_abs if "uncertainty_abs" in s.ncattrs() else np.nan
        r = s.uncertainty_rel if "uncertainty_rel" in s.ncattrs() else np.nan
        for j, l in enumerate(ls):
            fullnamel = s.bsh_name + '_' + level2fnheight[l]
            ds = dsd[fullnamel]
            ks = times.get_indexer(ds.index)
            if dtype is value:
                s[ks, j] = np.array(
                    list(zip(ds['Value'].values, ds['Quality'].values)), dtype)
                continue
            if dtype is avg:
                dsavg = ds['Value'].values
                ε_avg = s.uncertainty_propagated_avg * error_abs(dsavg, a, r)
                s[ks, j] = np.array(
                    list(zip(accround(dsavg, ε_avg),
                             ds['Quality'].values)), dtype)
                continue
            # in all other cases both std and avg are present
            dsstd = ds['Deviation'].values
            ε_std = np.fmax(s.uncertainty_propagated_std * a,
                            s.standard_error_std * dsstd)
            dsavg = ds['Value'].values
            ε_avg = np.fmax(s.uncertainty_propagated_avg * error_abs(dsavg, a, r),
                            s.standard_error_avg * dsstd)
            if dtype is avgstd:
                s[ks, j] = np.array(
                    list(zip(accround(dsavg, ε_avg), accround(dsstd, ε_std),
                             ds['Quality'].values)), dtype)
                continue
            # in all other cases max is present
            dsmax = ds['Maximum'].values
            ε_max = error_abs(dsmax, a, r)
            if dtype is maxavgstd:
                s[ks, j] = np.array(
                    list(zip(accround(dsmax, ε_max),
                             accround(dsavg, ε_avg), accround(dsstd, ε_std),
                             ds['Quality'].values)), dtype)
                continue
            # in the remaining case (allstat) min is also present
            dsmin = ds['Minimum'].values
            ε_min = error_abs(dsmin, a, r)
            if dtype is allstats:
                s[ks, j] = np.array(
                    list(zip(accround(dsmin, ε_min), accround(dsmax, ε_max),
                             accround(dsavg, ε_avg), accround(dsstd, ε_std),
                             ds['Quality'].values)), dtype)
                continue

# %%

f.close()
