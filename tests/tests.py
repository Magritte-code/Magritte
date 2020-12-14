
import numpy   as np
print('numpy  ', np.__version__)

import magritte.core as core
print(core.AMU)

import h5py

with h5py.File('models/test.hdf5') as file:
    file['data'] = np.array([0.0])
