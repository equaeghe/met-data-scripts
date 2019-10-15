import zipfile as zf
import pandas as pd
import numpy as np
import h5py as h5
import datetime as dt
import pytz
import signifl as sf  # https://github.com/equaeghe/signifl

#%% READ DATA FROM EXCEL FILES INTO PANDAS DATAFRAME

# The input file must be a zip file containing all the xls files available at
# https://www.noordzeewind.nl/kennis/meteogegevens.html
input_file = "OWEZ_M_181.xls-files.zip"

# The output file is an HDF5-file
output_file = "OWEZ.h5"

dsd = {}
with zf.ZipFile(input_file) as owezz:
    path, *files = owezz.namelist()
    for f in files:
        key = f.split("/")[1].split(".")[0] # foo/key.bar
        dsd[key] = pd.read_excel(owezz.open(f), skiprows=[0],
                                 index_col=[0,1,2,3,4,5], na_values=-99999)
    del key

df = pd.concat([dsd[key] for key in sorted(dsd.keys())])
del dsd

channels = list(range(0, 39)) + list(range(50, 59))
df.columns = [col + "." + str(channel) for channel in channels
              for col in ["Channel", "Max", "Min", "Mean", "StdDev"]]
df.index.names = ["year", "month", "day", "hour", "minute", "seconds"]


#%% CREATE HDF5 FILE AND ADD GENERAL METADATA

f = h5.File(output_file, "w")
f.attrs["date_metadata_modified"] = "2019-10-15"
f.attrs["HDF5-conversion"] = (
    "The original Excel data files (66, total 938 MB) were converted to this "
    "single HDF5 file using the Python packages Pandas (Excel data extraction),"
    " NumPy (data processing), and h5py (HDF5 manipulation) and the h5repack "
    "utility. The result is a much smaller size (74 MB), convenient data "
    "access, and the addition of descriptive annotations.")
f.attrs["annotations"] = (
    "The annotations attached to the groups and datasets in this file contain "
    "not only a general description of those groups or datasets, but also "
    "information about, for example, device identifiers, units, resolution, "
    "precision, sampling rate, error propagation, data storage details, etc.")
f.attrs["institution"] = "NoordzeeWind"
f.attrs["data_processing"] = (
    "The data as recorded in the source Excel files has been modified. "
    "Namely, the storage format, standard 8-byte 'float64' format in Excel, "
    "was reduced for storage and transmission efficiency purposes to 4-byte "
    "'float32' and individual values were stored in such a way as to encode "
    "their uncertainty given the information available about resolution and "
    "precision. This is a 'lossy' but compression-improving transformation. "
    "(N.B.: Some relevant information has not been found and in that case "
    "conservative bounds were assumed, so further storage size reduction is "
    "possible in case those data become available.)")
f.attrs["uncertainty_encoding"] = (
    "The statistics' values are recorded as floating point numbers for which "
    "the uncertainty is encoded using the following convention: if x is the "
    "number and e its uncertainty, then we store y = (2 * k ± 1) * d/2, where\n"
    "	* d = 2 ** floor(log2(e)) and\n"
    "	* k = floor(x / d).\n"
    "So y is the odd multiple of d/2 nearest to x. Because y is an odd "
    "multiple of a power of two, the denominator of y seen as an irreducable "
    "fraction is 2/d, so that d can effectively be determined from y.")
f.attrs["description"] = (
    "This file mainly contains the datasets of the met mast at the Egmond "
    "aan Zee Offshore Wind Farm (OWEZ), which was the first wind farm in the "
    "Netherlands to be built offshore. Measurements are available from July "
    "2005, before construction of the wind farm, untill December 2011. "
    "Construction started in April 2006.")
f.attrs["references"] = (
    "Data files "
    "<https://www.noordzeewind.nl/nl_nl/kennis/meteogegevens.html>\n"
    "User manual data files meteorological mast NoordzeeWind "
    "(version 2, 2007-10-01) by H. J. Kouwenhoven "
    "<https://www.noordzeewind.nl/nl_nl/kennis/meteogegevens/_jcr_content/par/"
    "iconlist/iconlistsection/link.stream/1554383874387/"
    "2f65120cc6dbeb967cfa34e08f95a95304e22e39/"
    "r03-manual-data-files-meteo-mast-noordzeewind.pdf>\n"
    "Surrounding obstacles influencing the OWEZ meteo mast measurements "
    "(version 2, 2007-08) by A. Curvers "
    "<https://www.noordzeewind.nl/nl_nl/kennis/meteogegevens/_jcr_content/par/"
    "iconlist/iconlistsection_1245724476/link.stream/1554383869782/"
    "bca17955510aad4a2be63909cf037499351150de/"
    "owez-r-181-t0-undisturbed-wind.pdf>\n"
    "Current Profiles at the Offshore Wind Farm Egmond aan Zee (2010-11-11) "
    "by J. W. Wagenaar & P. J. Eecen "
    "<https://www.ecn.nl/publicaties/ECN-E--10-076>\n"
    "3D Turbulence at the Offshore Wind Farm Egmond aan Zee (2010-10-08) "
    "by J. W. Wagenaar & P. J. Eecen "
    "<https://www.ecn.nl/publicaties/ECN-E--10-075>\n"
    "The OWEZ Meteorological Mast: Analysis of mast-top displacements (2008) "
    "by P.J. Eecen & E. Branlard "
    "<https://www.ecn.nl/publicaties/ECN-E--08-067>\n"
    "Personal Communication (2016-11/12) with Marc de Hoop "
    "of Mierij Meteo Systems B. V.\n"
    "Personal Communication (2016-12) with Peter J. Eecen of ECN\n"
    "Personal Communication (2016-12) with Sicco Kamminga of Nortek B. V.")
f.attrs["source"] = (
    "The original data are (at the time of writing) available from the OWEZ "
    "website <http://www.noordzeewind.nl/>. The data themselves are provided "
    "as per-month Excel files. Next to the data files themselves, there are "
    "multiple relevant reports available there as well.")
f.attrs["creator_name"] = "NoordzeeWind"
f.attrs["creator_type"] = "institution"
f.attrs["version"] = "2018-02-27"
f.attrs["publisher"] = "Erik Quaeghebeur"
f.attrs["publisher_email"] = "E.R.G.Quaeghebeur@tudelft.nl"
f.attrs["publisher_type"] = "person"
f.attrs["publisher_institution"] = "Delft University of Technology"


#%% ADD TIME INFORMATION TO GROUP IN HDF5 FILE

t = f.create_group("time")
t.attrs["description"] = (
    "This group contains the time instances of the datasets contained under "
    "the 'device' and 'channel' groups. 'timestamp' contains "
    "the full-precision information (up to second resolution), "
    "while the others allow for convenient time-based filtering of the data.")
dts = [dt.datetime(*ymdhms, tzinfo=pytz.UTC) for ymdhms in df.index.values]
t.create_dataset("timestamp", dtype="i4",
                 data=[dto.timestamp() for dto in dts],
                 shuffle=True, fletcher32=True,
                 compression="gzip", compression_opts=9)
t["timestamp"].attrs["description"] = (
    "The time instances encoded in unix time, the number of seconds "
    "that have elapsed since 1970-01-01T00:00:00 (UTC).")
for timeres in ["year", "month", "day", "hour", "minute"]:
    t.create_dataset(timeres, dtype="i4",
                     data=[getattr(dto, timeres) for dto in dts],
                     shuffle=True, fletcher32=True,
                     compression="gzip", compression_opts=9)
f.flush()

#%% PREPARE METADATA FOR DATASETS

channel_metadata = {
    0: {"channel": "00", "device": "WM",
        "location": "NW21", "quantity": "wd"},
    1: {"channel": "01", "device": "WM",
        "location": "NW21", "quantity": "hws"},
    2: {"channel": "02", "device": "WM",
        "location": "NW21", "quantity": "vws"},
    3: {"channel": "03", "device": "WM",
        "location": "NW116", "quantity": "wd"},
    4: {"channel": "04", "device": "WM",
        "location": "NW116", "quantity": "hws"},
    5: {"channel": "05", "device": "WM",
        "location": "NW116", "quantity": "vws"},
    6: {"channel": "06", "device": "WM",
        "location": "NW70", "quantity": "wd"},
    7: {"channel": "07", "device": "WM",
        "location": "NW70", "quantity": "hws"},
    8: {"channel": "08", "device": "WM",
        "location": "NW70", "quantity": "vws"},
    9: {"channel": "09", "device": "CA",
        "location": "NW116", "quantity": "ws"},
    10: {"channel": "10", "device": "CA",
         "location": "NE116", "quantity": "ws"},
    11: {"channel": "11", "device": "CA",
         "location": "S116", "quantity": "ws"},
    12: {"channel": "12", "device": "CA",
         "location": "NW70", "quantity": "ws"},
    13: {"channel": "13", "device": "CA",
         "location": "NE70", "quantity": "ws"},
    14: {"channel": "14", "device": "CA",
         "location": "S70", "quantity": "ws"},
    15: {"channel": "15", "device": "CA",
         "location": "NW21", "quantity": "ws"},
    16: {"channel": "16", "device": "CA",
         "location": "NE21", "quantity": "ws"},
    17: {"channel": "17", "device": "CA",
         "location": "S21", "quantity": "ws"},
    18: {"channel": "18", "device": "HMP",
         "location": "S116", "quantity": "at"},
    19: {"channel": "19", "device": "HMP",
         "location": "S70", "quantity": "at"},
    20: {"channel": "20", "device": "HMP",
         "location": "S21", "quantity": "at"},
    21: {"channel": "21", "device": "HMP",
         "location": "S116", "quantity": "rh"},
    22: {"channel": "22", "device": "HMP",
         "location": "S70", "quantity": "rh"},
    23: {"channel": "23", "device": "HMP",
         "location": "S21", "quantity": "rh"},
    24: {"channel": "24", "device": "RPT",
         "location": "T20", "quantity": "ap"},
    25: {"channel": "25", "device": "PD",
         "location": "NW70", "quantity": "pl"},
    26: {"channel": "26", "device": "PD",
         "location": "NE70", "quantity": "pl"},
    27: {"channel": "27", "device": "ST",
         "location": "U3_8", "quantity": "wt"},
    28: {"channel": "28", "device": "AC",
         "location": "T116", "quantity": "ma"},
    29: {"channel": "29", "device": "AC",
         "location": "T116", "quantity": "pa"},
    30: {"channel": "30", "device": "WV",
         "location": "NW116", "quantity": "wd"},
    31: {"channel": "31", "device": "WV",
         "location": "NE116", "quantity": "wd"},
    32: {"channel": "32", "device": "WV",
         "location": "S116", "quantity": "wd"},
    33: {"channel": "33", "device": "WV",
         "location": "NW70", "quantity": "wd"},
    34: {"channel": "34", "device": "WV",
         "location": "NE70", "quantity": "wd"},
    35: {"channel": "35", "device": "WV",
         "location": "S70", "quantity": "wd"},
    36: {"channel": "36", "device": "WV",
         "location": "NW21", "quantity": "wd"},
    37: {"channel": "37", "device": "WV",
         "location": "NE21", "quantity": "wd"},
    38: {"channel": "38", "device": "WV",
         "location": "S21", "quantity": "wd"},
    50: {"channel": "50", "device": "AWAC",
         "location": "U17", "quantity": "wl"},
    51: {"channel": "51", "device": "AWAC",
         "location": "U17", "quantity": "wt"},
    52: {"channel": "52", "device": "AWAC",
         "location": "U17", "quantity": "vh"},
    53: {"channel": "53", "device": "AWAC",
         "location": "U17", "quantity": "vp"},
    54: {"channel": "54", "device": "AWAC",
         "location": "U17", "quantity": "vd"},
    55: {"channel": "55", "device": "AWAC",
         "location": "U17", "quantity": "cv7"},
    56: {"channel": "56", "device": "AWAC",
         "location": "U17", "quantity": "cv11"},
    57: {"channel": "57", "device": "AWAC",
         "location": "U17", "quantity": "cd7"},
    58: {"channel": "58", "device": "AWAC",
         "location": "U17", "quantity": "cd11"}
}
device_metadata = {
    "AC": {"description": "acceleration sensor", "make": "Seika",
           "type": "XB2i", "owez_code": "AC SB2i",
           "reference": (
               "B1, B2, B3: Accelerometers of high overload resistance with "
               "integrated electronics for measurement of acceleration in "
               "the frequency range 0 to several 100 Hz (data sheet) from "
               "SEIKA Mikrosystemtechnik GmbH "
               "<http://www.seika.de/english/pdf_e/B_e.pdf>"),
           "note": (
               "The XB2i is an SB2i sensor box with two B accelerometers.")},
    "AWAC": {"description": "acoustic wave and current sensor",
             "make": "Nortek AS", "type": "AWAC", "owez_code": "ADCP",
             "reference": (
                 "AWAC: Acoustic Wave And Current Profiler (data sheet) from "
                 "Nortek AS "
                 "<http://www.nortek-as.com/lib/brochures/datasheet-awac>"),
             "note": (
                 "The version is unknown. However, from the generic October "
                 "2014 information brochure from Nortek about this line of "
                 "devices and the given depth of 17 m, we assume this device "
                 "uses a 1 MHz acoustic frequency. Based on that, we assume to "
                 "know from the brochure the sampling frequencies for the "
                 "different quantities; these are conservative in the sense "
                 "that they are the highest described as possible and "
                 "consequently the most precise Mean and StdDev estimates.")},
    "CA": {"description": "cup anemometer", "make": "Mierij Meteo",
           "type": "018", "owez_code": "WS 018",
           "reference": (
               "MW 11 Anemometer (data sheet) ""from Mierij Meteo B. V. "
               "<http://mierijmeteo.nl/"
               "wp-content/uploads/2014/07/MW-11-final.pdf>"),
           "note": "This corresponds to the current (2016) model MW 11."},
    "HMP": {"description": "relative humidity and ambient temperature sensors",
            "make": "Vaisala Oyj", "type": "HMP 233", "owez_code": "RHTT 261",
            "reference": (
                "HMP230 Series Transmitters User's Guide (2002-05) from "
                "Vaisala Oyj "
                "<http://www.vaisala.com/Vaisala%20Documents/"
                "User%20Guides%20and%20Quick%20Ref%20Guides/"
                                                      "HMP230_UserGuide.pdf>"),
                "note": (
                    "This is a transmitter with two separate sensors, "
                    "a Pt 100 sensor for temperature and a 'HUMICAP' sensor "
                    "for relative humidity.")},
    "PD": {"description": "precipitation sensor", "make": "Mierij Meteo",
           "type": "PD 205", "owez_code": "PD 205"},
    "RPT": {"description": "barometric pressure sensor", "make": "Druck",
            "type": "RPT 410V", "owez_code": "DP 910",
            "reference": (
                "RPT 410 Barometric Pressure Sensor (data sheet) "
                "from Druck Inc. "
                "<http://veronics.com/products/"
                "pressure_transducers-sensors/rpt410.pdf>")},
    "ST": {"description": "seawater temperature sensor",
           "make": "Mierij Meteo", "owez_code": "ST 808",
           "reference": (
               "SEM 1503/P & SEM1504/P DIN Rail Mounted ""Pt100 Transmitter "
               "(data sheet) from Status Instruments "
               "<http://www.status.co.uk/files/Products/411.pdf>"),
           "note": (
               "This is a Pt 100 sensor connected to "
               "a SEM1504/P transmitter.")},
    "WM": {"description": "3-axis ultrasonic anemometer", "make": "Gill",
           "type": "Windmaster 1086M", "owez_code": "3D WM4"},
    "WV": {"description": "wind vane", "make": "Mierij Meteo", "type": "524",
           "owez_code": "WD 524",
           "note": (
               "This likely corresponds to the current (2016) model MW 12.")}
}
quantity_metadata = {
    "AC": {
        "ma": {
            "description": "North-South acceleration",
            "note": ("The North-South or 'X' direction corresponds to movement "
                     "along a meridian."),
            "resolution": 0.01, # 0.001g (approx. 0.01 m/s2 with g = 9.81 m/s2);
                                # smallest difference between ordered quantized
                                # values is 0.007324218 m/s2
            "uncertainty_abs": 0.01, # resolution
            "sampling_frequency": 33,
            "samples": 19800,
            "unit": "m/s2"},
        "pa": {
            "description": "West-East acceleration",
            "note": ("The West-east or 'Y' direction corresponds to movement "
                     "along a parallel."),
            "resolution": 0.01, # 0.001g (approx. 0.01 m/s2 with g = 9.81 m/s2);
                                # smallest difference between ordered quantized
                                # values is 0.007324218 m/s2
            "uncertainty_abs": 0.01, # resolution
            "sampling_frequency": 33,
            "samples": 19800,
            "unit": "m/s2"}
    },
    "AWAC": {
        "cd11": {
            "description": "current direction 11 m",
            "resolution": 0.01,
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "note": (
                "From correspondence with the manufacturer, it seems likely "
                "that the 'raw' values (such as Max and Min) are based on "
                "acoustic measurements over 1 minute, so that the effective "
                "sampling frequency may well be lower than 1 Hz."),
            "sampling_frequency": 1,
            "samples": 600,
            "unit": "degree"},
        "cd7": {
            "description": "current direction 7 m",
            "resolution": 0.01,
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "note": (
                "From correspondence with the manufacturer, it seems likely "
                "that the 'raw' values (such as Max and Min) are based on "
                "acoustic measurements over 1 minute, so that the effective "
                "sampling frequency may well be lower than 1 Hz."),
            "sampling_frequency": 1,
            "samples": 600,
            "unit": "degree"},
        "cv11": {
            "uncertainty_rel": 0.01,
            "uncertainty_abs": 0.005,
            "description": "current velocity 11 m",
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "note": (
                "From correspondence with the manufacturer, it seems likely "
                "that the 'raw' values (such as Max and Min) are based on "
                "acoustic measurements over 1 minute, so that the effective "
                "sampling frequency may well be lower than 1 Hz."),
            "resolution": 0.01,
            "sampling_frequency": 1,
            "samples": 600,
            "unit": "m/s"},
        "cv7": {
            "uncertainty_rel": 0.01,
            "uncertainty_abs": 0.005,
            "description": "current velocity 7 m",
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "note": (
                "From correspondence with the manufacturer, it seems likely "
                "that the 'raw' values (such as Max and Min) are based on "
                "acoustic measurements over 1 minute, so that the effective "
                "sampling frequency may well be lower than 1 Hz."),
            "resolution": 0.01,
            "sampling_frequency": 1,
            "samples": 600,
            "unit": "m/s"},
        "vd": {
            "uncertainty_abs": 2, # 2–4°
            "description": "wave direction",
            "note": (
                "From correspondence with the manufacturer, we know that "
                "acoustic measurements over 512 s (8 minutes 32 s) were the "
                "basis for the 'raw' values (such as Max and Min), which means "
                "that the effective sampling frequency may well be lower than "
                "2 Hz. The uncertainty was also an assessment expressed in "
                "correspondence from the manufacturer."),
            "resolution": 0.01,
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "sampling_frequency": 2,
            "samples": 1200,
            "unit": "degree"},
        "vh": {
            "uncertainty_rel": 0.01,
            "description": "wave height",
            "note": (
                "From correspondence with the manufacturer, it seems very "
                "likely that this is the H_m0 wave height and not the "
                "classical definition of significant wave height; they also "
                "mentioned a 5 cm uncertainty for this quantity. Also, acoustic "
                "measurements over 512 s (8 minutes 32 s) were the basis for "
                "the 'raw' values (such as Max and Min), which means that the "
                "effective sampling frequency may well be lower than 4 Hz."),
            "resolution": 0.01,
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "sampling_frequency": 4,
            "samples": 2400,
            "unit": "m"},
        "vp": {
            "description": "wave period",
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "resolution": 0.01,
            "sampling_frequency": 2,
            "samples": 1200,
            "unit": "s"},
        "wl": {
            "description": "water level",
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "resolution": 0.01,
            "sampling_frequency": 4,
            "samples": 2400,
            "unit": "m"},
        "wt": {
            "uncertainty_abs": 0.1,
            "description": "water temperature",
            "resolution_stats": 0.01, # The recorded Mean & StdDev values have
                                      # the same resolution as Max & Min values.
            "range": [ -4.,  40.],
            "resolution": 0.01,
            "sampling_frequency": 1,
            "samples": 600,
            "unit": "degree Celsius"}
    },
    "CA": {
        "ws": {
            "uncertainty_abs": 0.5,
            "description": "wind speed",
            "range": [  0.,  50.],
            "sampling_frequency": 4,
            "samples": 2400,
            "starting_speed": 0.5,
            "unit": "m/s"}
    },
    "HMP": {
        "at": {
            "uncertainty_abs": 0.1, # at 20°C
            "description": "ambient temperature",
            "range": [-40.,  80.],
            "temperature_dependence": "0.005 degree Celsius/ degree Celsius",
            "unit": "degree Celsius"},
        "rh": {
            "uncertainty_abs": 1, # maximum: 1% (0–90%), 2% (90–100%);
                               # salt solutions: 2% (0–90%), 3% (90–100%);
                               # these are absolute values (percentage points),
                               # not relative values
            "description": "relative humidity",
            "range": [   0.,  100.],
            "unit": "%"}
    },
    "PD": {
        "pl": {
            "description": "precipitation level",
            "note": (
                "This is described as a 'yes'/'no' sensor. From correspondence "
                "with the manufacturer, I learned it is meant to be a "
                "'no'/'light'/'medium'/'heavy' sensor, but that in pratice it "
                "is difficult to decide between the latter three levels. The "
                "log-histogram of the data itself contains a sharp peak around "
                "1 that contains a large majority of the data points and a "
                "crudely uniform distribution of data between 1 and 4. "
                "Supposedly, 1 corresponds to 'no' and higher values to "
                "increasing intensity of precipitation; there are a few data "
                "points between 0 and 1, but these are likely faulty "
                "measurements."),
            "range": [ 0.,  5.],
            "resolution": 7.9e-5} # smallest difference between ordered
                                  # quantized values
    },
    "RPT": {
        "ap": {
            "uncertainty_abs": 0.5, # 0.5 mbar (around 20°C);
                                 # 1 mbar (between -10°C and 50°C)',
            "description": "ambient air pressure",
            "range": [600., 1100.],
            "resolution": 0.02,
            "unit": "mbar"}
        },
    "ST": {
        "wt": {
            "uncertainty_abs": 0.15,
            "uncertainty_rel": 0.001,
            "description": "water temperature",
            "unit": "degree Celsius"}
    },
    "WM": {
        "hws": {
            "uncertainty_rel": 0.015, #%RMS (<20 m/s)',
            "description": "horizontal wind speed",
            "uncertainty_abs": np.sqrt(0.01 ** 2 + 0.01 ** 2), # ‘resolution’ + ‘offset’
            "range": [  0.,  60.],
            "resolution": 0.01,
            "sampling_frequency": 4,
            "samples": 2400,
            "unit": "m/s"},
        "vws": {
            "uncertainty_rel": 0.015, #%RMS (<20 m/s)',
            "description": "vertical wind speed",
            "uncertainty_abs": np.sqrt(0.01 ** 2 + 0.01 ** 2), # ‘resolution’ + ‘offset’
            "range": [-60., 60.],
            "resolution": 0.01,
            "sampling_frequency": 4,
            "samples": 2400,
            "unit": "m/s"},
        "wd": {
            "uncertainty_abs": 2, # <25 m/s
            "description": "wind direction",
            "note": (
                "The majority of missing (Max, Min)-value pairs with "
                "non-missing Mean value occur when Mean is close to 0 degree or "
                "360 degree, so this data is 'MNAR', missing not at random, and is "
                "non-ignorable. (Possibly, such value pairs were removed to "
                "mistakenly avoid Min-values smaller than Max-values, although "
                "that makes perfect sense for directional data.)"),
            "range": [0., 360.],
            "resolution": 1,
            "sampling_frequency": 4,
            "samples": 2400,
            "unit": "degree"}
    },
    "WV": {
        "wd": {
            "uncertainty_abs": np.sqrt((360 / 256) ** 2 + 0.3 ** 2), # ‘resolution’ + ‘code disc error’
            "description": "wind direction",
            "note": (
                "The majority of missing (Max, Min)-value pairs with "
                "non-missing Mean value occur when Mean is close to 0 degree or "
                "360 degree, so this data is 'MNAR', missing not at random, and is "
                "non-ignorable. (Possibly, such value pairs were removed to "
                "mistakenly avoid Min-values smaller than Max-values, although "
                "that makes perfect sense for directional data.) The recorded "
                "Min and Max values do not reflect the stated resolution, even "
                "if we allow for Gray-code-to-analog and analog-to-digital "
                "conversion; we assume that whatever signal processing gave "
                "rise to the recorded values did not significantly decrease "
                "the uncertainty."),
            "range": [0., 360.],
            "resolution": 360 / 256, # ~1.4° (8 bit, so 360°/256)
            "unit": "degree"}
    }
}

#%% CREATE CHANNELS GROUP IN HDF5 FILE AND ADD DATASETS

c = f.create_group("channels")
for channel, metadata in channel_metadata.items():
    cha = metadata["channel"]
    dev = metadata["device"]
    qty = metadata["quantity"]
    ch = c.create_group(cha)
    ch.attrs.update(metadata)
    resolution = quantity_metadata[dev][qty].get("resolution", 0)
    resolution_stats = quantity_metadata[dev][qty].get("resolution_stats", 0)
    ε_abs = quantity_metadata[dev][qty].get("uncertainty_abs", 0)
    ε_rel = quantity_metadata[dev][qty].get("uncertainty_rel", 0)
    samples = quantity_metadata[dev][qty].get("samples", np.inf)
    data_mean = df["Mean." + str(channel)].values
    sq_ε_mean = ε_abs ** 2 + (ε_rel * data_mean) ** 2
    del data_mean
    data_stddev = df["StdDev." + str(channel)].values
    sq_ε_stddev = (ε_rel * data_stddev) ** 2
    del data_stddev
    for statistic in ["Max", "Min", "Mean", "StdDev"]:
        data = df[statistic + "." + str(channel)].values
        if statistic in {"Max", "Min"}:
            ε_abs = np.fmax(resolution, ε_abs)
            accuracy = np.sqrt(ε_abs ** 2 + (data * ε_rel) ** 2)
        elif statistic in {"Mean", "StdDev"}:
            if statistic == "Mean":
                sq_accuracy = sq_ε_mean + sq_ε_stddev
            elif statistic == "StdDev":
                sq_accuracy = sq_ε_mean + 3 * sq_ε_stddev
            accuracy = np.fmax(resolution_stats, np.sqrt(sq_accuracy / samples))
        ch.create_dataset(statistic, data=np.float32(sf.encode(data, accuracy)),
                          shuffle=True, fletcher32=True,
                          compression="gzip", compression_opts=9)

#%% CREATE USER-FRIENDLY GROUP STRUCTURE AND HARDLINK DATASETS THERE

d = f.create_group("devices")
for device, metadata in device_metadata.items():
    dev = d.create_group(device)
    dev.attrs.update(metadata)
    for quantity, metadata in quantity_metadata[device].items():
        qty = dev.create_group(quantity)
        qty.attrs.update(metadata)

for channel, metadata in channel_metadata.items():
    # create hard link to datasets already stored in channels group
    d[metadata["device"]][metadata["quantity"]][metadata["location"]] = \
                                                         c[metadata["channel"]]

#%%

f.close()
