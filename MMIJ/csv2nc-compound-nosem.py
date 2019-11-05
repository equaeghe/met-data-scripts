import netCDF4 as nc
import numpy as np
import pandas as pd
import collections as cl
from ruamel.yaml import YAML as yaml


def yaml_load(f):
    return yaml(typ='safe').load(f)


# %%

# The input file must be a csv file (possibly gzipped) that can be obtained via
# https://www.windopzee.net/meteomast-ijmuiden-mmij/data/
input_file = "MMIJ-ALLSTA.csv.gz"

# The output file is a netCDF4-file
output_file = "MMIJ-compound-nosem.nc"


pdf = pd.read_csv(input_file, encoding='iso-8859-15',
                  sep=';', header=[1,2], na_values=['NaN'],
                  parse_dates=True, index_col=0)

# %%

f = nc.Dataset(output_file, 'w')
with open('metadata-global.yaml') as g:
    for attribute, value in yaml_load(g).items():
        f.setncattr_string(attribute, value)
f.comment = yaml_load(""">-
    The statistics datasets have as values a compound data structure with the
    minimum ('min'), maximum ('max'), average ('avg'), and standard deviation
    ('std') of the samples within the given time interval.\n

    Custom dataset attributes:\n
    \t* 'quality_indicators' is a global attribute that describes the values of
    the (custom) 'quality' attribute.\n
    \t* 'uncertainty_abs' describes the absolute uncertainty of the sampled
    values.\n
    \t* 'uncertainty_rel' describes the relative uncertainty of the sampled
    values as a number between 0 and 1.\n
    \t* 'uncertainty_propagated' describes the uncertainty of the avg and std
    values as a fraction of the sampled values' (absolute) uncertainty, as
    derived from 'uncertainty_abs' and 'uncertainty_rel'.\n

    Custom instrument attributes:\n
    \t* 'manufacturer' describes the manufacturer of the instrument, typically
    by giving its name.\n
    \t* 'manufacturer_name' gives the manufacturer's (marketing) name for the
    instrument.\n
    \t* 'part_number' describes the manufacturer's part number for the
    instrument.
""")


np_stats = np.dtype([("min", 'f4'), ("max", 'f4'),
                     ("avg", 'f4'), ("std", 'f4')])
stats = f.createCompoundType(np_stats, "statistics")

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
def process_data(pdf, dss, data_out, attributes):
    print(dss)
    quantize = lambda r: 2 ** np.floor(np.log2(r) - 1)
    resolution = 10 ** -5
    uncertainty_rel = attributes.get('uncertainty_rel', 0.)
    uncertainty_abs = np.fmax(attributes.get('uncertainty_abs', 0.), resolution)
    data_in = {}
    for statistic in {'min', 'max', 'avg', 'std'}:
        data_in[statistic] = pdf[dss + statistic].values[:,0]
    for statistic in {'min', 'max'}:
        data = data_in[statistic]
        scale = quantize(np.fmax(data * uncertainty_rel, uncertainty_abs))
        data_out[statistic] = np.around(data / scale) * scale
    # propagate error to avg and std
    scale = quantize(np.fmax(uncertainty_abs/50, resolution))
    for statistic in {'avg', 'std'}:
        data_out[statistic] = np.around(data_in[statistic] / scale) * scale
    return data_out

def create_signal(group, instrument, name, dimensions, attributes):
    s = f[group][instrument].createVariable(name,
            stats, tuple(dimensions.keys()),
            zlib=True, complevel=9, fill_value=False)
    s.setncatts(attributes)
    # we can derive the cell method from the statistic name
    cell_method_details = (
        "(interval: 10 minutes comment: sampled at 4 Hz, so 2400 samples)")
    s.cell_methods = (
        "['min'] time: minimum " + cell_method_details + '\n'
        "['max'] time: maximum " + cell_method_details + '\n'
        "['avg'] time: mean " + cell_method_details + '\n'
        "['std'] time: standard_deviation " + cell_method_details)
    # we can derive the propagated error from the statistics name
    if 'uncertainty_abs' in s.ncattrs():
        s.uncertainty_propagated = (
            "['avg'] uncertainty/50 (50 approx. equal to sqrt(2400))" + '\n'
            "['std'] uncertainty/50 (50 approx. equal to sqrt(2400-1))")
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
                        s[:, h, l] = (np.nan, np.nan, np.nan, np.nan)
                        continue
                dss = ('MMIJ_' + 'H' + lev + 'B' + bm + '_' + name + '_' +
                       attributes['quality'] + '_')
                if dss.startswith('MMIJ_H85B000_WsMag'):
                    dss = dss.replace('WsMag', 'WSMag')
                s[:, h, l] = process_data(pdf, dss, s[:, h, l], attributes)
        elif location is not None:
            for loc in dimensions[location]:
                l = np.where(dimensions[location][:] == loc)[0][0]
                dss = ('MMIJ_' + 'H' + lev + '_' + name
                       + '_' + loc.item().decode()
                       + '_' + attributes['quality'] + '_')
                s[:, h, l] = process_data(pdf, dss, s[:, h, l], attributes)
        else:
            dss = ('MMIJ_' + 'H' + lev + '_' + name + '_' +
                   attributes['quality'] + '_')
            if dss.startswith('MMIJ_H85_WsHor'):
                dss = dss.replace('WsHor', 'Ws')
            s[:, h] = process_data(pdf, dss, s[:, h], attributes)


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
