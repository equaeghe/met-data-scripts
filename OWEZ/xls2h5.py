import zipfile as zf
import pandas as pd
import numpy as np
import h5py as h5
import datetime as dt
import pytz
from ruamel.yaml import YAML as yaml
import signifl as sf  # https://github.com/equaeghe/signifl


def yaml_load(f):
    return yaml(typ='safe').load(f)


# %% READ DATA FROM EXCEL FILES INTO PANDAS DATAFRAME

# The input file must be a zip file containing all the xls files available at
# https://www.noordzeewind.nl/kennis/meteogegevens.html
input_file = "OWEZ_M_181.xls-files.zip"

# The output file is an HDF5-file
output_file = "OWEZ.h5"

dsd = {}
with zf.ZipFile(input_file) as owezz:
    path, *files = owezz.namelist()
    for f in files:
        key = f.split("/")[1].split(".")[0]  # foo/key.bar
        dsd[key] = pd.read_excel(owezz.open(f), skiprows=[0],
                                 index_col=[0,1,2,3,4,5], na_values=-99999)
    del key

df = pd.concat([dsd[key] for key in sorted(dsd.keys())])
del dsd

channels = list(range(0, 39)) + list(range(50, 59))
df.columns = [col + "." + str(channel) for channel in channels
              for col in ["Channel", "Max", "Min", "Mean", "StdDev"]]
df.index.names = ["year", "month", "day", "hour", "minute", "seconds"]


# %% CREATE HDF5 FILE AND ADD GENERAL METADATA

f = h5.File(output_file, "w")
with open('metadata-global.yaml') as g:
    for attribute, value in yaml_load(g).items():
        f.attrs[attribute] = value


# %% ADD TIME INFORMATION TO GROUP IN HDF5 FILE

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


# %% PREPARE METADATA FOR DATASETS

with open('metadata-channels.yaml') as g:
    channel_metadata = yaml_load(g)
with open('metadata-devices.yaml') as g:
    device_metadata = yaml_load(g)
with open('metadata-quantities.yaml') as g:
    quantity_metadata = yaml_load(g)


# %% CREATE CHANNELS GROUP IN HDF5 FILE AND ADD DATASETS

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

# %% CREATE USER-FRIENDLY GROUP STRUCTURE AND HARDLINK DATASETS THERE

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

# %%

f.close()
