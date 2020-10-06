import os
import sys
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/../../../')

import numpy          as np
import magritte.setup as setup

from magritte.tools import dtypeSize, dtypeReal
from magritte.core  import Model, CC

dimension = 1
npoints   = 100
nrays     = 2
nspecs    = 5
nlspecs   = 1
nquads    = 100

dens = 1.0E+12        # [m^-3]
abun = 1.0E+08        # [m^-3]
temp = 4.5E+01        # [K]
turb = 0.0E+00        # [m/s]
dx   = 1.0E+04        # [m]
dv   = 1.0E+02 / CC   # [fraction of speed of light]


modelName    = 'constant_velocity_gradient_1d.hdf5'
linedataFile = '../../data/test.txt'


model = Model ()
model.parameters.set_model_name(modelName)
model.parameters.set_dimension (dimension)
model.parameters.set_npoints   (npoints)
model.parameters.set_nrays     (nrays)
model.parameters.set_nspecs    (nspecs)
model.parameters.set_nlspecs   (nlspecs)
model.parameters.set_nquads    (nquads)

model.geometry.points.position = np.array([[i*dx, 0, 0] for i in range(npoints)], dtypeReal)
model.geometry.points.velocity = np.array([[i*dv, 0, 0] for i in range(npoints)], dtypeReal)

model = setup.set_Delaunay_neighbor_lists (model)
model = setup.set_Delaunay_boundary       (model)
model = setup.set_uniform_rays            (model)
model = setup.set_linedata_from_LAMDA_file(model, linedataFile)
model = setup.set_quadrature              (model)

model.chemistry.species.abundance = [[     0.0,   abun,  dens,  0.0,      1.0] for _ in range(npoints)]
model.chemistry.species.symbol    =  ['dummy0', 'test',  'H2', 'e-', 'dummy1']

model.thermodynamics.temperature.gas   = temp * np.ones(npoints, dtypeReal)
model.thermodynamics.turbulence.vturb2 = turb * np.ones(npoints, dtypeReal)


model.write()
model.read ()

model.compute_spectral_discretisation ()
model.compute_LTE_level_populations   ()
model.compute_inverse_line_widths     ()
model.compute_radiation_field         ()

print(model.radiation.I)
