import zipfile as zf
import pandas as pd
import netCDF4 as nc
import numpy as np
import signifl as sf  # https://github.com/equaeghe/signifl
from collections import OrderedDict

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
f.date_metadata_modified = "2019-10-15"
f.Conventions = "CF-1.6"  # http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html
f.title = "Statistics from the offshore meteorological mast FINO1"
f.institution = "Bundesamt für Seeschifffahrt und Hydrographie (BSH)"
f.source = "meteorological mast"
f.history = ("Created on 2018-02-07 from FINO1 CSV files downloaded on "
             "2017-03-22 by Erik Quaeghebeur using a custom Python import "
             "script with zipfile, netCDF4, Pandas, and numpy modules.")
f.references = (
    "https://www.bsh.de/"
          "DE/THEMEN/Beobachtungssysteme/Messnetz-MARNET/FINO/fino_node.html\n"
    "FINO1_Metadaten_for_dissemination.pdf")
f.time_coverage_start = "2004-01-01T00:00Z"  # excluded
f.time_coverage_end = "2017-01-01T00:00Z"
f.comment = ("The statistics datasets have as values a compound data "
             "structure with as possible components the minimum ('min'), "
             "maximum ('max'), average ('avg'), and standard deviation "
             "('std') of the samples within the given time or a value ('val') "
             "characterizing the measurement for the time interval. The "
             "estimated uncertainty for the statistics is\n"
             "\t* the instrument's provided error for 'min' and 'max',\n"
             "\t* the propagated error for 'avg' and 'std'.\n\n"
             "The statistics' values are recorded as floating point numbers "
             "for which the uncertainty is encoded using the following "
             "convention: if x is the number and e its uncertainty, then we "
             "store y = (2 * k ± 1) * d/2, where\n"
             "	* d = 2 ** floor(log2(e)) and\n"
             "	* k = floor(x / d).\n"
             "So y is the odd multiple of d/2 nearest to x. Because y is an "
             "odd multiple of a power of two, the denominator of y seen as an "
             "irreducable fraction is 2/d, so that d can effectively be "
             "determined from y.\n\n"
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
             "\t* 'uncertainty_abs' describes the absolute uncertainty of the "
             "sampled values.\n"
             "\t* 'uncertainty_rel' describes the relative uncertainty of the "
             "sampled values as a number between 0 and 1.\n"
             "\t* 'sampling_frequency' describes the frequency in Hz with "
             "which the quantity has been sampled, which determines the number "
             "n of samples."
             "\n\n"
             "Custom instrument attributes:\n"
             "\t* 'manufacturer' describes the manufacturer of the "
             "instrument, typically by giving its name.\n"
             "\t* 'part_number' describes the manufacturer's part number for "
             "the instrument.")
f.creator_name = "Energieonderzoek Centrum Nederland (ECN)"
f.creator_type = "institution"
f.publisher = "Erik Quaeghebeur"
f.publisher_email = "E.R.G.Quaeghebeur@tudelft.nl"
f.publisher_type = "person"
f.publisher_institution = "Delft University of Technology"

# %%

flags = OrderedDict([("missing", -1), ("raw", 0),
                     ("doubtful_quality", 1), ("quality_controlled", 2)])
flag = f.createEnumType('i1', "flag", flags)

#%%

## GROUP HIERARCHY ##

# INSTRUMENT GROUPS
instruments = {
    'CA': {'long_name': "cup anemometer",
           'manufacturer': "Vector Instruments Windspeed Ltd.",
           'part_number': "A100LK/PC3/WR",
           'references': ("http://www.windspeed.co.uk/ws/index.php?"
                                    "option=displaypage&op=page&Itemid=52")},
    'HTT': {'long_name': "hygro-thermo transmitter",
            'manufacturer': "Thies", 'part_number': "1.1005.50.512",
            'references': ("https://www.thiesclima.com"
                    "/pdf/en/Products/Humidity-Electrical-devices/?art=696")},
    'PM': {'long_name': "precipitation monitor",
           'manufacturer': "Thies", 'part_number': "5.4103.10.000",
           'references': ("https://www.thiesclima.com"
               "/pdf/en/Products/Precipitation-Electrical-devices/?art=791")},
    'PS': {'long_name': "precipitation sensor",
           'manufacturer': "Thies", 'part_number': "5.4103.20.xxx",
           'references': ("https://www.thiesclima.com"
               "/pdf/en/Products/Precipitation-Electrical-devices/?art=794")},
    'PTB': {'long_name': "barometer",
            'manufacturer': "Vaisala", 'part_number': "PTB 100 A",
            'references': ("http://www.eol.ucar.edu/isf/facilities/isff"
                                        "/sensors/vaisala/ptb100/PTB100.pdf")},
    'PYR': {'long_name': "pyranometer",
            'manufacturer': "Kipp & Zonen", 'part_number': "CM11",
            'references': ("http://www.kippzonen.com/Download"
                                "/72/Manual-Pyranometers-CMP-series-English")},
    'UA': {'long_name': "ultrasonic anemometer",
           'manufacturer': "Gill Instruments", 'part_number': "R3-50",
           'references': ("http://gillinstruments.com"
                                                "/data/datasheets/r3-50.pdf")},
    'WV': {'long_name': "wind vane", 'manufacturer': "Thies",
           'part_number': "4.3120.22.012",
           'references': ("https://www.thiesclima.com"
                                    "/pdf/en/Products/Wind-Classic/?art=275")}
}
for i, a in instruments.items():
    f.createGroup(i)
    f[i].setncatts(a)

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

# PER-INSTRUMENT
# level
instruments_locations = {
    # instrument: (levels, measurement_height, position, orientation, mast distance)
    'CA': (range(1, 9), [34.1, 41.6, 51.6, 61.6, 71.6, 81.6, 91.6, 102.5],
           ["boom"] * 7 + ["vertical strut"],
           [143, 142, 140, 142, 143, 139, 135, np.nan],
           [6, 6, 5.5, 5.5, 4, 3, 3, np.nan]),
    'HTT': ([1, 2, 3, 5, 8], [34.9, 42.4, 52.4, 72.4, 101],
            ["mast internal"] * 5, None, None),
    'PM': ([0, 8], [23.7, 101.2], ["container top", "mast external"], [0, 180],
           [1.54, 0.07]),
    'PS': ([0], [23.7], ["container top"], [0], [1.54]),
    'PTB': ([0, 7], [21, 93],
            ["inside container", "mast internal"], None, None),
    'PYR': ([1, 7], [34.8, 93], ["individual boom"] * 2, [203, 180], [1, 0.4]),
    'UA': ([2, 4, 6], [42.1, 62.1, 82.1], ["boom"] * 3, [308, 308, 311],
           [1.08] * 3),
    'WV': ([1, 3, 5, 7], [34.1, 51.6, 71.6, 91.6],
           ["boom"] * 4, [307, 310, 307, 315], [0.65] * 4)}
for i, locations in instruments_locations.items():
    # levels
    levels = "level_" + i
    f[i].createDimension(levels, len(locations[0]))
    f[i].createVariable(levels, 'u1', (levels), fill_value=False)
    f[i][levels].compress = "level"
    f[i][levels][:] = np.array(locations[0], 'u1')
    # measurement height
    height = "height_" + i
    f[i].createVariable(height, 'f4', (levels), fill_value=False,
                        least_significant_digit=1)
    f[i][height][:] = np.array(locations[1], 'f4')
    f[i][height].setncatts(
        {'units': "m", 'positive': 'up', 'axis': 'Z', 'standard_name': 'height',
         'long_name': "measurement height above LAT (Lowest Astronomical Tide)"})
    # position
    if locations[2] is not None:
        position = "position_" + i
        f[i].createVariable(position, np.str, (levels), fill_value=False)
        f[i][position][:] = np.array(locations[2], np.str)
    # orientation
    if locations[3] is not None:
        orientation = "orientation_" + i
        f[i].createVariable(orientation, 'f4', (levels), fill_value=False,
                            least_significant_digit=0)
        f[i][orientation][:] = np.array(locations[3], 'f4')
        f[i][orientation].setncatts(
            {'units': "degree", 'long_name': "instrument orientation on mast"})
    # mast distance
    if locations[4] is not None:
        distance = "distance_" + i
        f[i].createVariable(distance, 'f4', (levels), fill_value=False,
                            least_significant_digit=1)
        f[i][distance][:] = np.array(locations[4], 'f4')
        f[i][orientation].setncatts(
            {'units': "m",
             'long_name': "horizontal distance to outer edge of mast"})

## SIGNALS ##
level2fnheight = {i: str(h) + "m"
                  for i, h in enumerate([20, 33, 40, 50, 60, 70, 80, 90, 100])}
instruments_quantities = {
    'CA': {
        ("Windgeschwindigkeit", "Wg"): ({'min', 'max', 'avg', 'std'}, {
            'long_name': "wind speed", 'standard_name': "wind_speed",
            'coordinates': "height_CA position_CA orientation_CA distance_CA",
            'units': "m/s", 'valid_range': np.float32([0.1, 75.]),
            'uncertainty_abs': np.float32(0.1),
            'uncertainty_rel': np.float32(0.01),
            'sampling_frequency': np.float32(1.)})},
    'HTT': {
        ("Luftfeuchte", "Lf"): ({'avg'}, {
            'long_name': "relative humidity",
            'standard_name': "relative_humidity",
            'coordinates': "height_HTT position_HTT",
            'units': "%", 'valid_range': np.float32([10., 100.]),
            'uncertainty_abs': np.float32(3.),
            'sampling_frequency': np.float32(1.)}),
        ("Lufttemperatur", "Lt"): ({'avg'}, {
            'long_name': "air temperature", 'standard_name': "air_temperature",
            'coordinates': "height_HTT position_HTT",
            'units': "degree Celsius", 'uncertainty_abs': np.float32(0.1),
            'sampling_frequency': np.float32(1.)})},
    'PM': {
        ("Niederschlag", "Ns"): ({'val'}, {
            'long_name': "precipitation presence",
            'coordinates': "height_PM position_PM orientation_PM distance_PM",
            'flag_values': np.uint8([0, 1]),
            'flag_meanings': "no yes"})},
    'PS': {
        ("Niederschlagsintensitaet", "Ni"): ({'avg'}, {
            'long_name': "intensity of precipitation",
            'standard_name': "lwe_precipitation_rate",
            'units': "mA",  # stands for mm/min via formula in manual
            'coordinates': "height_PS position_PS orientation_PS distance_PS",
            'valid_range': np.float32([4., 20.]),  # [0,10] mm/min
            'sampling_frequency': np.float32(1.)})},
    'PTB': {
        ("Luftdruck", "Ld"): ({'avg'}, {
            'long_name': "air pressure", 'standard_name': "air_pressure",
            'coordinates': "height_PTB position_PTB",
            'units': "hPa",
            'valid_range': np.float32([800., 1060.]),
            'uncertainty_abs': np.float32(0.3),
            'sampling_frequency': np.float32(1.)})},
    'PYR': {
        ("Globalstrahlung", "Gs"): ({'avg'}, {
            'long_name': "global radiation", 'units': "W/m2",
            'coordinates': ("height_PYR position_PYR "
                            "orientation_PYR distance_PYR"),
            'valid_range': np.float32([0., 4000.]),
            'uncertainty_rel': np.float32(0.03),  # from metadata sheets; not in specsheet
            'sampling_frequency': np.float32(1.)})},
    'UA': {
        ("Windgeschwindigkeit_U_Anemometer", "Wg"): (
          {'min', 'max', 'avg', 'std'}, {
            'long_name': "wind speed", 'standard_name': "wind_speed",
            'coordinates': "height_UA position_UA orientation_UA distance_UA",
            'units': "m/s", 'valid_range': np.float32([0., 45.]),
            'uncertainty_abs': np.float32(0.01),
            'uncertainty_rel': np.float32(0.01),
            'sampling_frequency': np.float32(50.)}),
        ("Windrichtung", "Wr"): ({'avg', 'std'}, {
            'long_name': "wind direction",
            'standard_name': "wind_from_direction",
            'coordinates': "height_UA position_UA orientation_UA distance_UA",
            'units': "degree", 'valid_range': np.float32([0, 359]),
            'uncertainty_abs': np.float32(1.),
            'sampling_frequency': np.float32(50.)})},
    'WV': {
        ("Windrichtung", "Wr"): ({'max', 'avg', 'std'}, {
            'long_name': "wind direction",
            'standard_name': "wind_from_direction",
            'coordinates': "height_WV position_WV orientation_WV distance_WV",
            'units': "degree", 'valid_range': np.float32([0, 360]),
            'uncertainty_abs': np.float32(2.),
            'sampling_frequency': np.float32(1.)})}
}

# %%

stat2cms = {'min': 'minimum', 'max': 'maximum',
            'avg': 'mean', 'std': 'standard_deviation'}

for i, q in instruments_quantities.items():
    for name, properties in q.items():
        stats = properties[0]
        g = f[i].createGroup(name[1])
        g.setncatts(properties[1])
        s = g.createVariable('flag', 'i1', ("time", "level_" + i),
                             zlib=True, complevel=9, fletcher32=True,
                             fill_value=flags["missing"])
        s.flag_values = np.int8(list(flags.values()))
        s.flag_meanings = ' '.join(flags.keys())
        for stat in stats:
            s = g.createVariable(stat, 'f4', ("time", "level_" + i),
                                 zlib=True, complevel=9, fletcher32=True,
                                 fill_value=-np.nan)
            # Fill cell_methods property as needed
            if stat is not 'val':
                T = np.float32(600.)  # number of seconds in 10 minutes
                ν = g.sampling_frequency
                n = int(ν * T)
                g.setncattr('number_of_observations', n)
                s.cell_methods = ("time: {} (interval: {} s)".format(
                                                        stat2cms[stat], 1 / ν))

# %%

resolution = np.float32(1e-2)

for i, q in instruments_quantities.items():
    ls = instruments_locations[i][0]  # levels where instrument can be found
    for name, properties in q.items():
        stats = properties[0]
        g = f[i][name[1]]
        a = g.uncertainty_abs if 'uncertainty_abs' in g.ncattrs() else 0
        r = g.uncertainty_rel if 'uncertainty_rel' in g.ncattrs() else 0
        for j, l in enumerate(ls):
            fullnamel = name[0] + '_' + level2fnheight[l]
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
