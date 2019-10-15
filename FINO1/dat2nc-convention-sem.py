import zipfile as zf
import pandas as pd
import netCDF4 as nc
import numpy as np

#%%

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

#%%

f = nc.Dataset(output_file, 'w')
f.date_metadata_modified = "2019-10-15"
f.Conventions = "CF-1.6" # http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html
f.title = "Statistics from the offshore meteorological mast Fino 1"
f.institution = "Bundesamt für Seeschifffahrt und Hydrographie (BSH)"
f.source = "meteorological mast"
f.history = ("Created on 2017-09-06 from Fino 1 CSV files downloaded on "
             "2017-03-22 by Erik Quaeghebeur using a custom Python import "
             "script with zipfile, netCDF4, Pandas, and numpy modules.")
f.references = (
    "http://www.bsh.de/de/Meeresdaten/Beobachtungen/MARNET-Messnetz"
                                                        "/FINO_1/index.jsp\n"
    "FINO1_Metadaten_for_dissemination.pdf")
f.time_coverage_start = "2004-01-01T00:00Z" # excluded
f.time_coverage_end = "2017-01-01T00:00Z"
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
             "\t* 'accuracy_abs' describes the absolute accuracy of the "
             "sampled values.\n"
             "\t* 'accuracy_rel' describes the relative accuracy of the "
             "sampled values as a number between 0 and 1.\n"
             "\t* 'sampling_frequency' describes the frequency in Hz with "
             "which the quantity has been sampled, which determines the number "
             "n of samples.\n"
             "\t* 'standard_error_avg' describes the standard error of the avg "
             "as a fraction 1/sqrt(n) of the std values.\n"
             "\t* 'standard_error_std' describes the standard error of the std "
             "as a fraction 1/sqrt(2·(n-1)) of the std values.\n"
             "\t* 'accuracy_propagated_avg' describes the accuracy of the avg "
             "values as a fraction 1/sqrt(n) of the sampled values' (absolute) "
             "accuracy, as derived from 'accuracy_abs' and 'accuracy_rel'.\n"
             "\t* 'accuracy_propagated_std' describes the accuracy of the std "
             "values as a fraction 1/sqrt(n-1) of the sampled values' (absolute) "
             "accuracy, as derived from 'accuracy_abs' and 'accuracy_rel'."
             "\n\n"
             "Custom instrument attributes:\n"
             "\t* 'manufacturer' describes the manufacturer of the instrument, "
             "typically by giving its name.\n"
             "\t* 'part_number' describes the manufacturer's part number for "
             "the instrument.")
f.creator_name = "Energieonderzoek Centrum Nederland (ECN)"
f.creator_type = "institution"
f.publisher = "Erik Quaeghebeur"
f.publisher_email = "E.R.G.Quaeghebeur@tudelft.nl"
f.publisher_type = "person"
f.publisher_institution = "Delft University of Technology"

#%%

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
    f[i].createVariable(height, 'u1', (levels), fill_value=False)
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
        f[i].createVariable(orientation, 'f4', (levels), fill_value=False)
        f[i][orientation][:] = np.array(locations[3], 'f4')
        f[i][orientation].setncatts(
            {'units': "degree", 'long_name': "instrument orientation on mast"})
    # mast distance
    if locations[4] is not None:
        distance = "distance_" + i
        f[i].createVariable(distance, 'f4', (levels), fill_value=False)
        f[i][distance][:] = np.array(locations[4], 'f4')
        f[i][orientation].setncatts(
            {'units': "m",
             'long_name': "horizontal distance to outer edge of mast"})

## SIGNALS ##
level2fnheight = {i: str(h) + "m"
                  for i, h in enumerate([20, 33, 40, 50, 60, 70, 80, 90, 100])}
instruments_quantities = {
    'CA': {
        ("Windgeschwindigkeit", "Wg"): (allstats, {
            'long_name': "wind speed", 'standard_name': "wind_speed",
            'units': "m/s", 'valid_range': np.array([0.1, 75.], dtype='f4'),
            'accuracy_abs': np.float32(0.1), 'accuracy_rel': np.float32(0.01),
            'sampling_frequency': np.float32(1.)})},
    'HTT': {
        ("Luftfeuchte", "Lf"): (avg, {
            'long_name': "relative humidity",
            'standard_name': "relative_humidity",
            'units': "%", 'valid_range': np.array([10., 100.], dtype='f4'),
            'accuracy_abs': np.float32(3.),
            'sampling_frequency': np.float32(1.)}),
        ("Lufttemperatur", "Lt"): (avg, {
            'long_name': "air temperature", 'standard_name': "air_temperature",
            'units': "degree Celsius", 'accuracy_abs': np.float32(0.1),
            'sampling_frequency': np.float32(1.)})},
    'PM': {
        ("Niederschlag", "Ns"): (value, {
            'long_name': "precipitation presence",
            'flag_values': np.array([0, 1], dtype='u1'),
            'flag_meanings': "no yes"})},
    'PS': {
        ("Niederschlagsintensitaet", "Ni"): (avg, {
            'long_name': "intensity of precipitation",
            'standard_name': "lwe_precipitation_rate",
            'units': "mA", # stands for mm/min via formula in manual
            'valid_range': np.array([4., 20.], dtype='u1'), # [0,10] mm/min
            'accuracy_abs': np.float32(0.005), # based on resolution in CSV
            'sampling_frequency': np.float32(1.)})},
    'PTB': {
        ("Luftdruck", "Ld"): (avg, {
            'long_name': "air pressure", 'standard_name': "air_pressure",
            'units': "hPa",
            'valid_range': np.array([800., 1060.], dtype='f4'),
            'accuracy_abs': np.float32(0.3),
            'sampling_frequency': np.float32(1.)})},
    'PYR': {
        ("Globalstrahlung", "Gs"): (avg, {
            'long_name': "global radiation", 'units': "W/m2",
            'valid_range': np.array([0., 4000.], dtype='f4'),
            'accuracy_rel': np.float32(0.03),  # TODO: unclear where 3% comes from
            'accuracy_abs': np.float32(0.005), # based on resolution in CSV
            'sampling_frequency': np.float32(1.)})},
    'UA': {
        ("Windgeschwindigkeit_U_Anemometer", "Wg"): (allstats, {
            'long_name': "wind speed", 'standard_name': "wind_speed",
            'units': "m/s", 'valid_range': np.array([0., 45.], dtype='f4'),
            'accuracy_abs': np.float32(0.01),
            'accuracy_rel': np.float32(0.01),
            'sampling_frequency': np.float32(50.)}),
        ("Windrichtung", "Wr"): (avgstd, {
            'long_name': "wind direction",
            'standard_name': "wind_from_direction",
            'units': "degree", 'valid_range': np.array([0, 359], dtype='f4'),
            'accuracy_abs': np.float32(1.),
            'sampling_frequency': np.float32(50.)})},
    'WV': {
        ("Windrichtung", "Wr"): (maxavgstd, {
            'long_name': "wind direction",
            'standard_name': "wind_from_direction",
            'units': "degree", 'valid_range': np.array([0, 360], dtype='f4'),
            'accuracy_abs': np.float32(2.),
            'sampling_frequency': np.float32(1.)})}
}

#%%

for i, q in instruments_quantities.items():
    for name, properties in q.items():
        dtype = properties[0]
        s = f[i].createVariable(name[1], dtype, ("time", "level_" + i),
                                zlib=True, complevel=9, fill_value=False)
        s.setncatts(properties[1])
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
                s.accuracy_propagated_avg = np.float32(1/np.sqrt(n))
            if dtype in {allstats, maxavgstd, avgstd}:
                s.cell_methods += ('\n' +
                    "['std'] time: standard_deviation " + cell_method_details)
                s.standard_error_avg = np.float32(1/np.sqrt(n))
                s.standard_error_std = np.float32(1/np.sqrt(2 * (n-1)))
                s.accuracy_propagated_std = np.float32(1/np.sqrt(n-1))
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

#%%

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

for i, q in instruments_quantities.items():
    ls = instruments_locations[i][0] # levels where instrument can be found
    for name, properties in q.items():
        dtype = properties[0]
        s = f[i][name[1]]
        a = s.accuracy_abs if "accuracy_abs" in s.ncattrs() else np.nan
        r = s.accuracy_rel if "accuracy_rel" in s.ncattrs() else np.nan
        for j, l in enumerate(ls):
            fullnamel = name[0] + '_' + level2fnheight[l]
            ds = dsd[fullnamel]
            ks = times.get_indexer(ds.index)
            if dtype is value:
                s[ks, j] = np.array(
                    list(zip(ds['Value'].values, ds['Quality'].values)), dtype)
                continue
            if dtype is avg:
                dsavg = ds['Value'].values
                ε_avg = s.accuracy_propagated_avg * error_abs(dsavg, a, r)
                s[ks, j] = np.array(
                    list(zip(accround(dsavg, ε_avg),
                             ds['Quality'].values)), dtype)
                continue
            # in all other cases both std and avg are present
            dsstd = ds['Deviation'].values
            ε_std = np.fmax(s.accuracy_propagated_std * a,
                            s.standard_error_std * dsstd)
            dsavg = ds['Value'].values
            ε_avg = np.fmax(s.accuracy_propagated_avg * error_abs(dsavg, a, r),
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

#%%

f.close()
