import netCDF4 as nc
import numpy as np
import scipy.stats as sps
import pandas as pd
import collections as cl
from ruamel.yaml import YAML as yaml

import signifl as sf  # https://github.com/equaeghe/signifl


def yaml_load(f):
    return yaml(typ='safe').load(f)


# %%

# The input file must be a csv file (possibly gzipped) that can be obtained via
# https://www.windopzee.net/meteomast-ijmuiden-mmij/data/
input_file = "MMIJ-ALLSTA.csv.gz"

# The output file is a netCDF4-file
output_file = "MMIJ.nc"


pdf = pd.read_csv(input_file, encoding='iso-8859-15',
                  sep=';', header=[1,2], na_values=['NaN'],
                  parse_dates=True, index_col=0)

# %%

f = nc.Dataset(output_file, 'w')
with open('metadata-global.yaml') as g:
    for attribute, value in yaml_load(g).items():
        f.setncattr_string(attribute, value)
f.comment = yaml_load(""">-
    Custom dataset attributes:
    \t* 'quality_indicators' is a global attribute that describes the values of
    the (custom) 'quality' attribute.\n
    \t* 'uncertainty_abs' describes the absolute uncertainty of the sampled
    values.\n
    \t* 'uncertainty_rel' describes the relative uncertainty of the sampled
    values as a number between 0 and 1.\n

    Custom instrument attributes:\n
    \t* 'manufacturer' describes the manufacturer of the instrument, typically
    by giving its name.\n
    \t* 'manufacturer_name' gives the manufacturer's (marketing) name for the
    instrument.\n
    \t* 'part_number' describes the manufacturer's part number for the
    instrument.\n

    The statistics' values are recorded as floating point numbers for which the
    uncertainty is encoded using the following convention: if x is the number
    and e its uncertainty, then we store y = (2 * k ± 1) * d/2, where\n
    \t* d = 2 ** floor(log2(e)) and\n
    \t* k = floor(x / d).\n
    So y is the odd multiple of d/2 nearest to x. Because y is an odd multiple
    of a power of two, the denominator of y seen as an irreducable fraction is
    2/d, so that d can effectively be determined from y.
""")

# %%

## GROUP HIERARCHY ##
with open('structure.yaml') as g:
    structure = yaml_load(g)

with open('metadata-instruments.yaml') as g:
    instrument_metadata = yaml_load(g)

for c, c_data in structure['groups'].items():  # MEASUREMENT GROUPS
    f.createGroup(c)
    for d in c_data['groups']:  # INSTRUMENT GROUPS
        f[c].createGroup(d)
        if d in instrument_metadata:
            f[c][d].setncatts(instrument_metadata[d])

## DIMENSIONS & COORDINATES ##

def create_coordinate(group, coordinate, data):
    v = group.createVariable(coordinate, data['dtype'],
                             tuple(data['dimensions']), **data['options'])
    del data['dtype']
    del data['dimensions']
    del data['options']
    if 'values' in data:
        values = np.array(data['values'])
        if values.dtype == 'f8':
            values = np.array(values, 'f4')
        v[:] = values
        del data['values']
    for attribute, value in data.items():
        v.setncattr_string(attribute, value)

dimension2coordinate = {}

# global
for dimension, length in structure['dimensions'].items():
    f.createDimension(dimension, length)
    dimension2coordinate[dimension] = dimension
for coordinate, data in structure['coordinates'].items():
    create_coordinate(f, coordinate, data)

# per-group or per-instrument
for c, c_data in structure['groups'].items():  # MEASUREMENT GROUPS
    # per-group coordinates
    if 'dimensions' in c_data:
        for dimension, length in c_data['dimensions'].items():
            f[c].createDimension(dimension, length)
            dimension2coordinate[dimension] = '/'.join([c, dimension])
    if 'coordinates' in c_data:
        for coordinate, data in c_data['coordinates'].items():  # GROUP COORDINATES
            create_coordinate(f[c], coordinate, data)
    # per-instrument coordinates
    for d, d_data in c_data['groups'].items():  # INSTRUMENT GROUPS
        if 'dimensions' in d_data:
            for dimension, length in d_data['dimensions'].items():
                f[c][d].createDimension(dimension, length)
                dimension2coordinate[dimension] = '/'.join([c, d, dimension])
        if 'coordinates' in d_data:
            for coordinate, data in d_data['coordinates'].items():  # INSTRUMENT COORDINATES
                create_coordinate(f[c][d], coordinate, data)

# ‘manual’ data entry
# time
f["time"][:] = nc.date2num(pdf.index.to_pydatetime(),
                           "seconds since 1970-01-01")
f["time_bounds"][:, 0] = f["time"][:]
f["time_bounds"][:, 1] = f["time"][:] + 600

# %%

## SIGNALS ##
# processing functions
def process_data(pdf, dss, g, h, l):
    resolution = np.float32(10e-5)
    ε_rel = g.getncattr('uncertainty_rel') if 'uncertainty_rel' in g.ncattrs() else 0.
    ε_abs = g.getncattr('uncertainty_abs') if 'uncertainty_abs' in g.ncattrs() else 0.
    data_avg = pdf[dss + 'avg'].values[:, 0]
    sq_ε_avg = ε_abs ** 2 + (ε_rel * data_avg) ** 2
    del data_avg
    data_std = pdf[dss + 'std'].values[:, 0]
    sq_ε_std = (ε_rel * data_std) ** 2
    data_min = pdf[dss + 'min'].values[:, 0]
    data_max = pdf[dss + 'max'].values[:, 0]
    z = sps.norm.ppf(1 - 1/2400)  # quantile corresponding to 1/2400 exceedence
    sq_τ = (np.fmin(data_max - data_min, z * data_std) / 2 / 2400) ** 2
    del data_std
    del data_min
    del data_max
    for statistic in {'min', 'max', 'avg', 'std'}:
        data = pdf[dss + statistic].values[:, 0]
        if statistic in {'min', 'max'}:
            scale = np.fmax(np.sqrt(sq_τ + ε_abs ** 2 + (data * ε_rel) ** 2),
                            resolution / 2)
        if statistic in {'avg', 'std'}:
            if statistic == 'avg':
                sq_accuracy = sq_ε_avg + sq_ε_std
            elif statistic == 'std':
                sq_accuracy = sq_ε_avg + 3 * sq_ε_std
            scale = np.fmax(np.sqrt(sq_τ + sq_accuracy / 2400), resolution / 2)
        if l is None:
            g[statistic][:, h] = sf.encode(data, scale)
        else:
            g[statistic][:, h, l] = sf.encode(data, scale)


def create_signal(group, instrument, name, dimensions, attributes):
    statistics = {'min': 'minimum', 'max': 'maximum',
                  'avg': 'mean', 'std': 'standard_deviation'}
    g = f[group][instrument].createGroup(name)
    g.setncatts(attributes)
    chunksizes = [1] * len(dimensions)
    chunksizes[0] = len(f['time'])
    for statistic in statistics.keys():
        s = g.createVariable(statistic, 'f4', tuple(dimensions.keys()),
                             chunksizes=chunksizes,
                             zlib=True, complevel=9, fill_value=-np.nan)
        # we can derive the cell method from the statistic name
        s.cell_methods = "time: {} (interval: 0.25 s)".format(
                                                        statistics[statistic])
    if name.startswith('True'):
        name = name[4:]
    boom = None
    location = None
    for dim in dimensions:
        if dim.startswith('level'):
            level = dim
        if dim.startswith('boom'):
            boom = dim
        if dim.startswith('location'):
            location = dim
    for lev in dimensions[level]:
        h = np.where(dimensions[level][:] == lev)[0][0]
        if boom is not None:
            for bm in dimensions[boom]:
                l = np.where(dimensions[boom][:] == bm)[0][0]
                if name == 'Ws':
                    if (((lev != '92') and (bm in {'180', '300'})) or
                        ((bm not in {'180', '300'}) and (lev == '92'))):
                        continue
                dss = ('MMIJ_' + 'H' + lev + 'B' + bm + '_' + name + '_' +
                        attributes['quality'] + '_')
                if dss.startswith('MMIJ_H85B000_WsMag'):
                    dss = dss.replace('WsMag', 'WSMag')
                print(dss, h, l)
                process_data(pdf, dss, g, h, l)
        elif location is not None:
            for loc in dimensions[location]:
                l = np.where(dimensions[location][:] == loc)[0][0]
                dss = ('MMIJ_' + 'H' + lev + '_' + name
                       + '_' + loc.item().decode()
                       + '_' + attributes['quality'] + '_')
                print(dss, h, l)
                process_data(pdf, dss, g, h, l)
        else:
            dss = ('MMIJ_' + 'H' + lev + '_' + name + '_' +
                    attributes['quality'] + '_')
            if dss.startswith('MMIJ_H85_WsHor'):
                dss = dss.replace('WsHor', 'Ws')
            print(dss, h)
            process_data(pdf, dss, g, h, None)


# load and process the data and metadata
with open('metadata-quantities.yaml') as g:
    quantities_metadata = yaml_load(g)

for c, c_data in structure['groups'].items():  # MEASUREMENT GROUPS
    for d, d_data in c_data['groups'].items():  # INSTRUMENT GROUPS
        for v, v_data in d_data['variables'].items():  # VARIABLES
            dimensions = cl.OrderedDict(
                ((dimension, f[dimension2coordinate[dimension]])
                 for dimension in v_data['dimensions'])
            )
            del v_data['dimensions']
            attributes = quantities_metadata[d][v]
            attributes.update(v_data)
            for k, val in attributes.items():
                if k in {'valid_range', 'uncertainty_abs', 'uncertainty_rel'}:
                    attributes[k] = np.float32(val)
            create_signal(c, d, v, dimensions, attributes)

# %%

f.close()
