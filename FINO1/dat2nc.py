import zipfile as zf
import pandas as pd
import netCDF4 as nc
import numpy as np
from ruamel.yaml import YAML as yaml
import signifl as sf  # https://github.com/equaeghe/signifl
from collections import OrderedDict


def yaml_load(f):
    return yaml(typ='safe').load(f)


# %%

# The input file must be a zip file with dat files that can be obtained via
# http://fino.bsh.de/
input_file = "FINO1_20040101_20161231-with_fixes.zip"

# The output file is a netCDF4-file
output_file = "FINO1.nc"


dsd = {}
with zf.ZipFile(input_file) as finoz:
    path, *files = finoz.namelist()
    for f in [f for f in files if f.endswith('.dat')]:
        key = f.split('/')[1].split('.')[0]  # foo/key.bar
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

# %%

flags = OrderedDict([("missing", -1), ("raw", 0),
                     ("doubtful_quality", 1), ("quality_controlled", 2)])
flag = f.createEnumType('i1', "flag", flags)

#%%

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
f['time'].standard_name = 'time'
f['time'].long_name = 'time'
f["time"][:] = nc.date2num(times.to_pydatetime(), "seconds since 1970-01-01")
f["time"].units = "seconds since 1970-01-01"
f["time"].calendar = "standard"
f["time"].delta_t = "0000-00-00 00:10:00"
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
f.createVariable("height", 'f4', ("level",), fill_value=False,
                 least_significant_digit=0)
f["height"][:] = np.float32([21, 33.5, 41, 51, 61, 71, 81, 91, 101])
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
    f[instrument][levels][:] = np.uint8(locations['levels'])
    # measurement height
    height = "height_" + instrument
    f[instrument].createVariable(height, 'f4', (levels), fill_value=False,
                        least_significant_digit=1)
    f[instrument][height][:] = np.float32(locations['measurement_height'])
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
        f[instrument].createVariable(orientation, 'f4', (levels),
                                     fill_value=False,
                                     least_significant_digit=0)
        f[instrument][orientation][:] = np.float32(locations['orientation'])
        f[instrument][orientation].setncatts(
            {'units': "degree", 'long_name': "instrument orientation on mast"})
    # mast distance
    if locations['mast_distance'] is not None:
        distance = "distance_" + instrument
        f[instrument].createVariable(distance, 'f4', (levels), fill_value=False,
                            least_significant_digit=1)
        f[instrument][distance][:] = np.float32(locations['mast_distance'])
        f[instrument][orientation].setncatts(
            {'units': "m",
             'long_name': "horizontal distance to outer edge of mast"})

# %%

## SIGNALS ##
level2fnheight = {i: str(h) + "m"
                  for i, h in enumerate([20, 33, 40, 50, 60, 70, 80, 90, 100])}
instrument_coordinates = {
    'CA': "height_CA position_CA orientation_CA distance_CA",
    'HTT': "height_HTT position_HTT",
    'PM': "height_PM position_PM orientation_PM distance_PM",
    'PS': "height_PS position_PS orientation_PS distance_PS",
    'PTB': "height_PTB position_PTB",
    'PYR': "height_PYR position_PYR orientation_PYR distance_PYR",
    'UA': "height_UA position_UA orientation_UA distance_UA",
    'WV': "height_WV position_WV orientation_WV distance_WV"
}
with open('metadata-quantities.yaml') as g:
    quantities_metadata = yaml_load(g)

# %%

stat2cms = {'min': 'minimum', 'max': 'maximum',
            'avg': 'mean', 'std': 'standard_deviation'}

quantities_statistics = {}
for instrument, quantities in quantities_metadata.items():
    quantities_statistics[instrument] = {}
    for name, properties in quantities.items():
        quantities_statistics[instrument][name] = set(properties['statistics'])
        del properties['statistics']
        g = f[instrument].createGroup(name)
        g.setncattr_string('coordinates', instrument_coordinates[instrument])
        for attribute, value in properties.items():
            if attribute == 'flag_values':
                g.setncattr(attribute, np.uint8(value))
            elif attribute in {'valid_range', 'sampling_frequency',
                               'uncertainty_abs', 'uncertainty_rel'}:
                g.setncattr(attribute, np.float32(value))
            else:
                g.setncattr_string(attribute, value)
        s = g.createVariable('flag', 'i1', ("time", "level_" + instrument),
                             zlib=True, complevel=9, fletcher32=True,
                             fill_value=flags["missing"])
        s.flag_values = np.int8(list(flags.values()))
        s.flag_meanings = ' '.join(flags.keys())
        for stat in quantities_statistics[instrument][name]:
            s = g.createVariable(stat, 'f4', ("time", "level_" + instrument),
                                 zlib=True, complevel=9, fletcher32=True,
                                 fill_value=-np.nan)
            # Fill cell_methods property as needed
            if stat != 'val':
                T = np.float32(600.)  # number of seconds in 10 minutes
                ν = g.sampling_frequency
                n = int(ν * T)
                g.setncattr('number_of_observations', n)
                s.cell_methods = ("time: {} (interval: {} s)".format(
                                                        stat2cms[stat], 1 / ν))

# %%

resolution = np.float32(1e-2)

for instrument, quantities in quantities_metadata.items():
    ls = instruments_locations[instrument]['levels']
    for name in quantities:
        stats = quantities_statistics[instrument][name]
        g = f[instrument][name]
        a = g.uncertainty_abs if 'uncertainty_abs' in g.ncattrs() else 0
        r = g.uncertainty_rel if 'uncertainty_rel' in g.ncattrs() else 0
        for j, l in enumerate(ls):
            fullnamel = g.bsh_name + '_' + level2fnheight[l]
            ds = dsd[fullnamel]
            ks = times.get_indexer(ds.index)
            # flags
            g['flag'][ks, j] = ds['Quality'].values
            # statistics
            # first only value
            if stats == {'val'}:
                g['val'][ks, j] = ds['Value'].values
                continue
            # in all other cases, avg is present
            samples = g.sampling_frequency * 600  # number of samples in 10 minutes
            dsavg = ds['Value'].values
            if stats == {'avg'}:
                ε_avg = np.fmax(np.sqrt((a ** 2 + (r * dsavg) ** 2) / samples),
                                resolution / 2)
                g['avg'][ks, j] = sf.encode(dsavg, ε_avg)
                del ε_avg
                continue
            # in all other cases also std is present
            dsstd = ds['Deviation'].values
            ε_avg = np.fmax(np.sqrt((a ** 2 + (r * dsavg) ** 2 + 3 * (r ** dsstd) ** 2)
                                    / samples),
                            resolution / 2)
            g['avg'][ks, j] = sf.encode(dsavg, ε_avg)
            ε_std = np.fmax(np.sqrt((a ** 2 + (r * dsavg) ** 2 + (r ** dsstd) ** 2)
                                    / samples),
                            resolution / 2)
            g['std'][ks, j] = sf.encode(dsstd, ε_std)
            del ε_std
            for stat in (stats - {'avg', 'std'}):
                statname = 'Maximum' if stat == 'max' else 'Minimum'
                dsext = ds[statname].values
                ε_ext = np.fmax(np.sqrt(a ** 2 + (r * dsext) ** 2),
                                resolution / 2)
                g[stat][ks, j] = sf.encode(dsext, ε_ext)
                del ε_ext

# %%

f.close()
