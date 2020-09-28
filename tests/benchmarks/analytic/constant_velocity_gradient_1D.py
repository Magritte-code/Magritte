import os
import sys
sys.path.append(f'{os.path.dirname(os.path.realpath(__file__))}/../../../')

from magritte.core import Model, Io

npoints   = 100
nrays     = 2
nspecs    = 5
nlspecs   = 1
nquads    = 100

CC = 299792458.0

dens = 1.0E+12        # [m^-3]
abun = 1.0E+08        # [m^-3]
temp = 4.5E+01        # [K]
turb = 0.0E+00        # [m/s]
dx   = 1.0E+04        # [m]
dv   = 1.0E+02 / CC   # [fraction of speed of light]

model = Model ()

model.parameters.set_npoints (npoints)
model.parameters.set_nrays   (nrays)
model.parameters.set_nspecs  (nspecs)
model.parameters.set_nlspecs (nlspecs)
model.parameters.set_nquads  (nquads)


print(model.thermodynamics.temperature.gas)
# print(model.geometry.points.position)

# model.geometry.points.x  = np.array([i*dx for i in range(npoints)])
# model.geometry.points.y  = np.array([0.0  for i in range(npoints)])
# model.geometry.points.z  = np.array([0.0  for i in range(npoints)])
#
# model.geometry.points.vx = np.array([i*dv for i in range(npoints)])
# model.geometry.points.vy = np.array([0.0  for i in range(npoints)])
# model.geometry.points.vz = np.array([0.0  for i in range(npoints)])
