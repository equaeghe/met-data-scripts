import netCDF4 as nc
import numpy as np
import scipy.stats as sps
import pandas as pd
import collections as cl
import signifl as sf  # https://github.com/equaeghe/signifl


# The input file must be a csv file (possibly gzipped) that can be obtained via
# https://www.windopzee.net/meteomast-ijmuiden-mmij/data/
input_file = "MMIJ-ALLSTA.csv.gz"

# The output file is a netCDF4-file
output_file = "MMIJ.nc"


pdf = pd.read_csv(input_file, encoding='iso-8859-15',
                  sep=';', header=[1,2], na_values=['NaN'],
                  parse_dates=True, index_col=0)

f = nc.Dataset(output_file, 'w')
f.Conventions = "CF-1.7" # http://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html
f.title = "Statistics from the offshore meteorological mast IJmuiden"
f.history = ("Created from MMIJ-ALLSTA.csv of 2016-12-21 on 2018-05-29 "
             "by Erik Quaeghebeur using a custom Python import script "
             "with netCDF4, Pandas, and numpy modules.")
f.references = "http://www.meteomastijmuiden.nl/data/instrumentatierapport/"
f.time_coverage_start = "2011-11-01T00:00Z"
f.time_coverage_end = "2016-03-12T00:00Z" # excluded
f.comment = ("The statistics datasets have as values a compound data structure"
             "with the minimum ('min'), maximum ('max'), average ('avg'), and "
             "standard deviation ('std') of the samples within the given time "
             "interval."
             "\n\n"
             "Custom dataset attributes:\n"
             "\t* 'quality_indicators' is a global attribute that describes "
             "the values of the (custom) 'quality' attribute.\n"
             "\t* 'uncertainty_abs' describes the absolute uncertainty of the "
             "sampled values.\n"
             "\t* 'uncertainty_rel' describes the relative uncertainty of the "
             "sampled values as a number between 0 and 1."
             "\n\n"
             "Custom instrument attributes:\n"
             "\t* 'manufacturer' describes the manufacturer of the instrument, "
             "typically by giving its name.\n"
             "\t* 'manufacturer_name' gives the manufacturer's (marketing) "
             "name for the instrument.\n"
             "\t* 'part_number' describes the manufacturer's part number for "
             "the instrument."
             "\n\n"
             "The statistics' values are recorded as floating point numbers "
             "for which the uncertainty is encoded using the following "
             "convention: if x is the number and e its uncertainty, then we "
             "store y = (2 * k ± 1) * d/2, where\n"
             "	* d = 2 ** floor(log2(e)) and\n"
             "	* k = floor(x / d).\n"
             "So y is the odd multiple of d/2 nearest to x. Because y is an "
             "odd multiple of a power of two, the denominator of y seen as an "
             "irreducable fraction is 2/d, so that d can effectively be "
             "determined from y.")
f.quality_indicators = ( # this and the quality attribute used are non-standard
    "Q1 - ISO 17025 approved signal quality, in accordance with IEC61400-12\n"
    "Q2 - Signal measured under QA of another MEASNET member\n"
    "Q3 - Signal measured under external QA, checked by ECN "
                                                "or other MEASNET member.\n"
    "Q4 - Signal measured under external QA, but not checked.\n"
    "Q5 - Signal not calibrated or calibration not checked")

## GROUP HIERARCHY ##

# MEASUREMENT GROUPS
for c in {"air", "wind", "prec"}: f.createGroup(c)

# INSTRUMENT GROUPS
for d in {'HMP', 'PTB', 'pseudo'}: f["air"].createGroup(d)
for d in {'CA', 'WV', 'USA'}: f["wind"].createGroup(d)
for d in {'PM', 'PD'}: f["prec"].createGroup(d)


## DIMENSIONS ##

# GLOBAL
# auxiliary
f.createDimension("bounds_edges", 2)
# time
f.createDimension("time", len(pdf))
f.createVariable('time', 'u4', ("time",), zlib=True, complevel=9,
                 fill_value=False)
f["time"][:] = nc.date2num(pdf.index.to_pydatetime(),
                           "seconds since 1970-01-01")
f['time'].standard_name = 'time'
f['time'].long_name = 'time'
f["time"].units = "seconds since 1970-01-01"
f["time"].calendar = "standard"
f["time"].delta_t = "0000-00-00 00:10:00"
f["time"].axis = 'T'
f["time"].bounds = 'time_bounds'
# time bounds (of statistics, denote closed-open intervals)
f.createVariable('time_bounds', 'u4', ("time", "bounds_edges"),
                 zlib=True, complevel=9, fill_value=False)
f["time_bounds"][:, 0] = f["time"][:]
f["time_bounds"][:, 1] = f["time"][:] + 600
f["time_bounds"].units = "seconds since 1970-01-01"
f["time_bounds"].long_name = ("bounds of (closed-open) time intervals "
                              "over which the statistics have been calculated")
# boom/direction
f.createDimension("boom", 3)
f.createVariable("boom", 'S3', ("boom",), fill_value=False)
f["boom"][:] = np.array(['000', '120', '240'], dtype='S3')
f["boom"].long_name = "the label of the boom on which the instrument is placed"
f.createVariable("direction", 'f4', ("boom",), fill_value=False)
f["direction"][:] = np.float32([46.5, 166.5, 286.5])
f["direction"].units = "degree"
f["direction"].long_name = ("the direction of the boom on which the "
                            "instrument is placed")


# AIR DIMENSIONS
# level/height
f["air"].createDimension("level_air", 2)
f["air"].createVariable("level_air", 'S2', ("level_air",), fill_value=False)
f["air/level_air"][:] = np.array(['21', '90'], dtype=str)
f["air/level_air"].axis = 'Z'
f["air"].createVariable("height_air", 'f4', ("level_air",), fill_value=False)
f["air/height_air"].setncatts(
    {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
     'long_name': "height above LAT (Lowest Astronomical Tide)"})
f["air/height_air"][:] = np.float32([21, 90])

# WIND/CA DIMENSIONS
# level/height
f["wind/CA"].createDimension("level_CA", 3)
f["wind/CA"].createVariable("level_CA", 'S2', ("level_CA",), fill_value=False)
f["wind/CA/level_CA"][:] = np.array(['27', '58', '92'], dtype='S2')
f["wind/CA/level_CA"].axis = 'Z'
f["wind/CA"].createVariable("height_CA", 'f4', ("level_CA",), fill_value=False)
f["wind/CA/height_CA"].setncatts(
    {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
     'long_name': "height above LAT (Lowest Astronomical Tide)"})
f["wind/CA/height_CA"][:] = np.float32([27, 58.5, 92])
# boom/direction
f["wind/CA"].createDimension("boom_CA", 5)
f["wind/CA"].createVariable("boom_CA", 'S3', ("boom_CA",), fill_value=False)
f["wind/CA/boom_CA"][:] = np.array(['000', '120', '180', '240', '300'],
                                   dtype='S3')
f["wind/CA"].createVariable("direction_CA", 'f4', ("boom_CA",),
                            fill_value=False)
f["wind/CA/direction_CA"].units = "degree"
f["wind/CA/direction_CA"][:] = np.float32([46.5, 166.5, 226.5, 286.5, 346.5])
## (level, boom)-mask
#f["wind/CA"].createDimension("present_CA", 8)
#level_boom_CA = np.dtype([("level", 'S2'), ("boom", 'S3')])
#level_boom_CA_t = f["wind/CA"].createCompoundType(level_boom_CA,
                                                  #"level_boom_CA")
#f["wind/CA"].createVariable("present_CA", (level_boom_CA_t,), ("present_CA"),
                            #fill_value=False)
#f["wind/CA/present_CA"][:] = np.array([('27', '000'), ('27', '120'),
                                       #('27', '240'), ('58', '000'),
                                       #('58', '120'), ('58', '240'),
                                       #('92', '180'), ('92', '300')],
                                      #dtype=level_boom_CA)
#f["wind/CA/present_CA"].long_name = ("Instruments were present and so data is "
                                     #"stored for the level/boom-combinations "
                                     #"listed in this variable")

# WIND/WV DIMENSIONS
# level/height
f["wind/WV"].createDimension("level_WV", 3)
f["wind/WV"].createVariable("level_WV", 'S2', ("level_WV",), fill_value=False)
f["wind/WV/level_WV"][:] = np.array(['27', '58', '87'], dtype='S2')
f["wind/WV/level_WV"].axis = 'Z'
f["wind/WV"].createVariable("height_WV", 'f4', ("level_WV",), fill_value=False)
f["wind/WV/height_WV"].setncatts(
    {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
     'long_name': "height above LAT (Lowest Astronomical Tide)"})
f["wind/WV/height_WV"][:] = np.float32([27, 58.5, 87])

# WIND/USA DIMENSIONS
# level/height
f["wind/USA"].createDimension("level_USA", 1)
f["wind/USA"].createVariable("level_USA", 'S2', ("level_USA",),
                             fill_value=False)
f["wind/USA/level_USA"][:] = np.array(['85'], dtype='S2')
f["wind/USA/level_USA"].axis = 'Z'
f["wind/USA"].createVariable("height_USA", 'f4', ("level_USA",),
                             fill_value=False)
f["wind/USA/height_USA"].setncatts(
    {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
     'long_name': "height above LAT (Lowest Astronomical Tide)"})
f["wind/USA/height_USA"][:] = np.float32([85])

# PREC/PD DIMENSIONS
# level/height
f["prec/PD"].createDimension("level_PD", 1)
f["prec/PD"].createVariable("level_PD", 'S2', ("level_PD",), fill_value=False)
f["prec/PD/level_PD"][:] = np.array(['27'], dtype='S2')
f["prec/PD/level_PD"].axis = 'Z'
f["prec/PD"].createVariable("height_PD", 'f4', ("level_PD",), fill_value=False)
f["prec/PD/height_PD"].setncatts(
    {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
     'long_name': "height above LAT (Lowest Astronomical Tide)"})
f["prec/PD/height_PD"][:] = np.float32([27])
# location
f["prec/PD"].createDimension("location_PD", 1)
f["prec/PD"].createVariable("location_PD", 'S1', ("location_PD",),
                            fill_value=False)
f["prec/PD/location_PD"][:] = np.array(['U'], dtype='S1')

# PREC/PM DIMENSIONS
# level/height
f["prec/PM"].createDimension("level_PM", 1)
f["prec/PM"].createVariable("level_PM", 'S2', ("level_PM",), fill_value=False)
f["prec/PM/level_PM"][:] = np.array(['21'], dtype=str)
f["prec/PM/level_PM"].axis = 'Z'
f["prec/PM"].createVariable("height_PM", 'f4', ("level_PM",), fill_value=False)
f["prec/PM/height_PM"].setncatts(
    {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
     'long_name': "height above LAT (Lowest Astronomical Tide)"})
f["prec/PM/height_PM"][:] = np.float32([21])
# location
f["prec/PM"].createDimension("location_PM", 2)
f["prec/PM"].createVariable("location_PM", 'S1', ("location_PM",),
                            fill_value=False)
f["prec/PM/location_PM"][:] = np.array(['l', 'r'], dtype='S1')

## SIGNALS ##
def process_data(pdf, dss, g, h, l):
    resolution = np.float32(10e-5)
    ε_rel = g.getncattr('uncertainty_rel') if 'uncertainty_rel' in g.ncattrs() else 0.
    ε_abs = g.getncattr('uncertainty_abs') if 'uncertainty_abs' in g.ncattrs() else 0.
    data_avg = pdf[dss + 'avg'].values[:,0]
    sq_ε_avg = ε_abs ** 2 + (ε_rel * data_avg) ** 2
    del data_avg
    data_std = pdf[dss + 'std'].values[:,0]
    sq_ε_std = (ε_rel * data_std) ** 2
    data_min = pdf[dss + 'min'].values[:,0]
    data_max = pdf[dss + 'max'].values[:,0]
    z = sps.norm.ppf(1 - 1/2400)  # quantile corresponding to 1/2400 exceedence
    sq_τ = (np.fmin(data_max - data_min, z * data_std) / 2 / 2400) ** 2
    del data_std
    del data_min
    del data_max
    for statistic in {'min', 'max', 'avg', 'std'}:
        data = pdf[dss + statistic].values[:,0]
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
                if name is 'Ws':
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
                dss = ('MMIJ_' + 'H' + lev + '_' + name + '_' +
                        loc.decode() + '_' + attributes['quality'] + '_')
                print(dss, h, l)
                process_data(pdf, dss, g, h, l)
        else:
            dss = ('MMIJ_' + 'H' + lev + '_' + name + '_' +
                    attributes['quality'] + '_')
            if dss.startswith('MMIJ_H85_WsHor'):
                dss = dss.replace('WsHor', 'Ws')
            print(dss, h)
            process_data(pdf, dss, g, h, None)

# AIR
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_air", f["air/level_air"])])
## HMP
f["air/HMP"].setncatts({
    'manufacturer': 'Vaisala',
    'manufacturer_name': 'HUMICAP',
    'part_number': 'HMP155',
    'references': ('http://www.vaisala.com/Vaisala%20Documents/'
        'Brochures%20and%20Datasheets/HMP155-Datasheet-B210752EN-E-LoRes.pdf'),
    'long_name': 'humidity and temperature probe'})
### Rh
#HMP Rh {'calibration_uncertainty': '(20 degree Celsius) 0.6 pp between 0% and 40%, 1 pp between 40% and 97%.', 'accuracy': '15 … 25 degree Celsius: 1 pp between 0% and 90%, 1.7 pp between 90% and 100%; -20 … 40 degree Celsius: (1 + 0.008·reading) pp [pp = percentage point]'}
attributes = {'long_name': "relative humidity",
              'standard_name': "relative_humidity",
              'units': "%", 'coordinates': "height_air",
              'quality': 'Q1',
              'valid_range': np.float32([0.0, 100.0]),
              'uncertainty_abs': np.float32(1.)} # there are some values > 100
create_signal("air", "HMP", "Rh", dimensions, attributes)
### Tair
#HMP Tair {'accuracy': '-80 … 20 degree Celsius: (0.176 - 0.0028·reading) degree Celsius; 20 … 60 degree Celsius: (0.07 + 0.0025·reading) degree Celsius', 'sensor': 'Pt100 RTD Class F0.1 IEC 60751'}
attributes = {'long_name': "air temperature",
              'standard_name': "air_temperature",
              'units': "degree Celsius",
              'coordinates': "height_air boom direction",
              'quality': 'Q1',
              'valid_range': np.float32([-80.0, 60.0]),
              'uncertainty_abs': np.float32(0.12)}
create_signal("air", "HMP", "Tair", dimensions, attributes)
## PTB
f["air/PTB"].setncatts({
    'manufacturer': 'Vaisala',
    'part_number': 'PTB210',
    'references': ('http://www.vaisala.com/Vaisala%20Documents/'
        'Brochures%20and%20Datasheets/PTB210-Datasheet-B210942EN-B-LoRes.pdf'),
    'long_name': 'digital barometer'})
### Pair
attributes = {'long_name': "air pressure", 'standard_name': "air_pressure",
              'units': "hPa",
              'coordinates': "height_air boom direction",
              'quality': 'Q1',
              'valid_range': np.float32([500.0, 1100.0]),
              'uncertainty_abs': np.float32(0.1)}
create_signal("air", "PTB", "Pair", dimensions, attributes)
## pseudo
### AirDensity
air_density = lambda P, T: P * 100 / 287.05 / (T + 273.15)
attributes = {'long_name': "calculated air density",
              'standard_name': "air_density", 'units': "kg/m3",
              'coordinates': "height_air boom direction",
              'quality': 'Q1',
              'valid_range': np.float32([air_density(500, 60),
                                       air_density(1100, -80)]),
              'uncertainty_abs': np.float32(air_density(0.1, 60))}
create_signal("air", "pseudo", "AirDensity", dimensions, attributes)

# WIND
## CA
f["wind/CA"].setncatts({
    'manufacturer': 'Adolf Thies GmbH & Co. KG',
    'manufacturer_name': 'Wind Transmitter “First Class” Advanced',
    'part_number': '4.3351.00.140',
    'references': ('https://www.thiesclima.com/en/Products/'
                                                'Wind-First-Class/?art=219'),
    'long_name': 'cup anemometer'})
### Ws
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_CA", f["wind/CA/level_CA"]),
                             ("boom_CA", f["wind/CA/boom_CA"])])
attributes = {'long_name': "wind speed", 'standard_name': "wind_speed",
              'units': "m/s",
              'coordinates': "height_CA direction_CA",
              'quality': 'Q1', 'valid_range': np.float32([0.3, 75.]),
              'uncertainty_abs': np.float32(0.2),
              'uncertainty_rel': np.float32(0.01)}
create_signal("wind", "CA", "Ws", dimensions, attributes)
### TrueWs
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_CA", f["wind/CA/level_CA"])])
attributes = {'long_name': "tower shadow-corrected wind speed",
              'standard_name': "wind_speed", 'units': "m/s",
              'coordinates': "height_CA",
              'quality': 'Q1',
              'valid_range': np.float32([0.3, 75.]),
              'uncertainty_abs': np.float32(0.2 / np.sqrt(2))}
create_signal("wind", "CA", "TrueWs", dimensions, attributes)
## WV
f["wind/WV"].setncatts({
    'manufacturer': 'Adolf Thies GmbH & Co. KG',
    'manufacturer_name': 'Thies Wind Vane First Class',
    'part_number': 'P 6200H',
    'references': 'http://www.wind-pgc.com/files/en_ds_windvane_thies_first.pdf',
    'long_name': 'wind vane'})
### Wd
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_WV", f["wind/WV/level_WV"]),
                             ("boom", f["boom"])])
attributes = {'long_name': "wind direction",
              'standard_name': "wind_from_direction",
              'units': "degree",
              'coordinates': "height_WV direction",
              'quality': 'Q1',
              'uncertainty_abs': np.float32(1.),
              'comment': (
                  "The avg values lie in [0,360], but the min and max values "
                  "may lie outside this interval, to make sure that "
                  "min < avg < max.")}
create_signal("wind", "WV", "Wd", dimensions, attributes)
### TrueWd
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_WV", f["wind/WV/level_WV"])])
attributes = {'long_name': "tower shadow-corrected wind direction",
              'standard_name': "wind_from_direction", 'units': "degree",
              'coordinates': "height_WV", 'quality': 'Q1',
              'uncertainty_abs': np.float32(1. / np.sqrt(2)),
              'comment': (
                  "The avg values lie in [0,360], but the min and max values "
                  "may lie outside this interval, to make sure that "
                  "min < avg < max.")}
create_signal("wind", "WV", "TrueWd", dimensions, attributes)
## USA
f["wind/USA"].setncatts({
    'manufacturer': 'METEK GmbH',
    'manufacturer_name': 'previously: USA-1; now: uSonic-3 Scientific)',
    'references': ('http://metek.de/wp-content/uploads/2014/05/'
        'Metek-Ultrasonic-Anemometer-uSonic-3-Scientific-USA-1-Datasheet.pdf'),
    'long_name': 'ultrasonic wind sensor'})
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_USA", f["wind/USA/level_USA"]),
                             ("boom", f["boom"])])
## SSon
attributes = {'long_name': "sonic status",
              'coordinates': "height_USA direction",
              'quality': 'Q1',
              'flag_values': np.float32([-1000., 0.]),
              'flag_meanings': "Not_OK OK"}
create_signal("wind", "USA", "SSon", dimensions, attributes)
### WsXSon
attributes = {'long_name': "wind speed X", 'standard_name': "x_wind",
              'units': "m/s",
              'coordinates': "height_USA direction",
              'quality': 'Q1',
              'ancillary_variables': "SSon",
              'valid_range': np.float32([-60.,  60.]),
              'uncertainty_abs': np.float32(0.1),
              'uncertainty_rel': np.float32(0.02)}
create_signal("wind", "USA", "WsXSon", dimensions, attributes)
### WsYSon
attributes = {'long_name': "wind speed Y", 'standard_name': "y_wind",
              'units': "m/s",
              'coordinates': "height_USA direction",
              'quality': 'Q1',
              'ancillary_variables': "SSon",
              'valid_range': np.float32([-60.,  60.]),
              'uncertainty_abs': np.float32(0.1),
              'uncertainty_rel': np.float32(0.02)}
create_signal("wind", "USA", "WsYSon", dimensions, attributes)
### WsZSon
attributes = {'long_name': "wind speed Z",
              'standard_name': "upward_air_velocity",
              'units': "m/s",
              'coordinates': "height_USA direction",
              'quality': 'Q1',
              'ancillary_variables': "SSon",
              'valid_range': np.float32([-60.,  60.]),
              'uncertainty_abs': np.float32(0.1),
              'uncertainty_rel': np.float32(0.02)}
create_signal("wind", "USA", "WsZSon", dimensions, attributes)
### WsMag
attributes = {'long_name': "calculated wind speed magnitude", 'units': "m/s",
              'coordinates': "height_USA direction",
              'quality': 'Q1',
              'valid_range': np.float32([0.,  np.sqrt(3) * 60.]),
              'uncertainty_abs': np.float32(0.1 * np.sqrt(3 - 8/np.pi))} # Maxwell-Boltzmann assumed
create_signal("wind", "USA", "WsMag", dimensions, attributes)
### WsHor
attributes = {'long_name': "calculated wind speed",
              'standard_name': "wind_speed", 'units': "m/s",
              'coordinates': "height_USA direction",
              'quality': 'Q1',
              'valid_range': np.float32([0.,  np.sqrt(2) * 60.]),
              'uncertainty_abs': np.float32(0.1 * np.sqrt(2 - np.pi/2))} # Rayleigh assumed
create_signal("wind", "USA", "WsHor", dimensions, attributes)
### TrueWsHor
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_USA", f["wind/USA/level_USA"])])
attributes = {'long_name': "tower shadow-corrected wind speed",
              'standard_name': "wind_speed", 'units': "m/s",
              'coordinates': "height_USA",
              'quality': 'Q1',
              'valid_range': np.float32([0.,  np.sqrt(2) * 60.]),
              'uncertainty_abs': np.float32(0.1 * np.sqrt(1 - np.pi/4))}
create_signal("wind", "USA", "TrueWsHor", dimensions, attributes)

# PREC
## PD
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_PD", f["prec/PD/level_PD"]),
                             ("location_PD", f["prec/PD/location_PD"])])
### Prec
attributes = {'long_name': "precipitation presence",
              'coordinates': "height_PD", 'quality': 'Q1',
              'flag_values': np.float32([0., 100.]),
              'flag_meanings': "no yes"}
create_signal("prec", "PD", "Prec", dimensions, attributes)
## PM
f["prec/PM"].setncatts({
    'manufacturer': 'Adolf Thies GmbH & Co. KG',
    'manufacturer_name': 'Thies Disdro Laser precipitation monitor',
    'part_number': '5.4110.00.xxx',
    'references': ('https://www.thiesclima.com/pdf/en/Products/'
                                'Precipitation-Electrical-devices/?art=774'),
    'long_name': 'laser precipitation monitor'})
dimensions = cl.OrderedDict([("time", f["time"]),
                             ("level_PM", f["prec/PM/level_PM"]),
                             ("location_PM", f["prec/PM/location_PM"])])
### qual
#PM qual {'resolution': '1 pp', 'description': 'momentaneous measurement quality', 'range': array([   0.,  100.])}
attributes = {'long_name': "momentaneous measurement quality", 'units': "%",
              'coordinates': "height_PM", 'quality': 'Q5'}
create_signal("prec", "PM", "qual", dimensions, attributes)
### Synop
#PM Synop {'resolution': '1', 'values': '51, 52, 53: Drizzle (also freezing); 58, 59: Drizzle with rain; 61, 63, 65: Rain (also freezing); 68, 69: Rain and/or Drizzle with snow; 71, 73, 75: Snow; 87, 88: Ice pellets, Soft hail; 77: Snow grains (also ice prisms), Ice crystals/needles; 89, 90: Hail', 'description': 'synoptic code (ww 4677), specifies the kind of precipitation, if any'}
attributes = {'long_name': ("synoptic code (ww 4677), "
                            "specifies the kind of precipitation, if any"),
              'coordinates': "height_PM", 'quality': 'Q5',
              'ancillary_variables': "OK qual"}
create_signal("prec", "PM", "Synop", dimensions, attributes)
### intens
#PM intens {'resolution_tentative': '0.016599655151399162 mm/min or about 1 mm/h', 'range': array([  5.00000000e-03,   2.50000000e+02])}
attributes = {'long_name': "intensity of precipitation",
              'standard_name': "lwe_precipitation_rate", 'units': "mm/min",
              'coordinates': "height_PM", 'quality': 'Q5',
              'ancillary_variables': "OK qual",
              'uncertainty_rel': np.float32(0.15)}
create_signal("prec", "PM", "intens", dimensions, attributes)
### Prec
attributes = {'long_name': "precipitation presence",
              'coordinates': "height_PM", 'quality': 'Q5',
              'flag_values': np.float32([0., 100.]),
              'flag_meanings': "no yes",
              'ancillary_variables': "OK qual"}
create_signal("prec", "PM", "Prec", dimensions, attributes)
### amount
#PM amount {'resolution': '1 mm'}
attributes = {'long_name': "amount of precipitation since the last sensor reset",
              'units': "mm", 'coordinates': "height_PM",
              'quality': 'Q5', 'valid_min': np.float32(0.),
              'ancillary_variables': "OK qual"}
create_signal("prec", "PM", "amount", dimensions, attributes)
### visib
#PM visib {'resolution': '10 m'}
attributes = {'long_name': "visibility distance",
              'standard_name': "visibility_in_air", 'units': "m",
              'coordinates': "height_PM", 'quality': 'Q5',
              'valid_range': np.float32([0., 100000.]),
              'ancillary_variables': "OK qual"}
create_signal("prec", "PM", "visib", dimensions, attributes)
### OK
attributes = {'long_name': "sensor status", 'units': "%",
              'coordinates': "height_PM", 'quality': 'Q1',
              'flag_values': np.float32([0., 100.]),
              'flag_meanings': "OK not_OK"}
create_signal("prec", "PM", "OK", dimensions, attributes)

f.close()
