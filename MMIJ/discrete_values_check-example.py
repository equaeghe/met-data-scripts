# This file contains an example of a brief sequence of commands typically
# used in a live session, to check whether the finite set of values appearing
# for some series actually correspond to valid values

import netCDF4 as nc  # library to manipulate netCDF4 files
import numpy as np  # library for scientific computing


# %% load our netCDF4 version of the MMIJ dataset
f = nc.Dataset('MMIJ_compound-sem.nc', 'r')

# %% extract precipitation sensor Max values
x = f['prec/PM/Synop'][:]['max']

# %% extract list of all the different non-NaN values in the dataset
np.unique(x[~np.isnan(x)])

# Output:
# array([-998., -997., -953., -952., -950., -900., -176.,  -16.,    0.,
#          51.,   53.,   55.,   58.,   59.,   61.,   63.,   65.,   68.,
#          69.,   71.,   73.,   75.,   77.,   87.,   88.,   89.,   90.,
#         108.], dtype=float32)
