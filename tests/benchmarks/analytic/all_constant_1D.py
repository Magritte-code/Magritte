import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f'{curdir}/../../../')


from tools import LTEpop, profile, lineEmissivity, lineOpacity, lineSource, I_CMB


import numpy             as np
import magritte.setup    as setup
import magritte.core     as magritte


dimension = 1
npoints   = 1000
nrays     = 2
nspecs    = 5
nlspecs   = 1
nquads    = 1

dens = 1.0E+12                 # [m^-3]
abun = 1.0E+10                 # [m^-3]
temp = 4.5E+00                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+04                 # [m]
dv   = 0.0E+00 / magritte.CC   # [fraction of speed of light]

modelFile = f'{curdir}/all_constant_1D.hdf5'
lamdaFile = f'{curdir}/../../data/test.txt'

model = magritte.Model ()
model.parameters.set_model_name(modelFile)
model.parameters.set_dimension (dimension)
model.parameters.set_npoints   (npoints)
model.parameters.set_nrays     (nrays)
model.parameters.set_nspecs    (nspecs)
model.parameters.set_nlspecs   (nlspecs)
model.parameters.set_nquads    (nquads)

model.geometry.points.position.set([[i*dx, 0, 0] for i in range(npoints)])
model.geometry.points.velocity.set([[i*dv, 0, 0] for i in range(npoints)])

model.chemistry.species.abundance = [[     0.0,   abun,  dens,  0.0,      1.0] for _ in range(npoints)]
model.chemistry.species.symbol    =  ['dummy0', 'test',  'H2', 'e-', 'dummy1']

model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

model = setup.set_Delaunay_neighbor_lists (model)
model = setup.set_Delaunay_boundary       (model)
model = setup.set_boundary_condition_CMB  (model)
model = setup.set_uniform_rays            (model)
model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
model = setup.set_quadrature              (model)

model.write()
model.read ()

model.compute_spectral_discretisation ()
model.compute_inverse_line_widths     ()
model.compute_LTE_level_populations   ()
model.compute_radiation_field         ()


x  = np.array(model.geometry.points.position)[:,0]
nu = np.array(model.radiation.frequencies.nu)
I  = np.array(model.radiation.I)


ld = model.lines.lineProducingSpecies[0].linedata

k = 0

frq = ld.frequency[k]
pop = LTEpop         (ld, temp) * abun
phi = profile        (ld, k, temp, (turb/magritte.CC)**2, frq)
eta = lineEmissivity (ld, pop)[k] * phi
chi = lineOpacity    (ld, pop)[k] * phi
src = lineSource     (ld, pop)[k]
bdy = I_CMB          (frq)


def I_0 (x):
    return src + (bdy-src)*np.exp(-chi*x)


def I_1 (x):
    return src + (bdy-src)*np.exp(-chi*(x[-1]-x))


def relative_error (a,b):
    return 2.0*(a-b)/(a+b)


error_I_0 = relative_error (I_0(x), I[0,:,0])
error_I_1 = relative_error (I_1(x), I[1,:,0])

print('max error in I(0) =', np.max(error_I_0))
print('max error in I(1) =', np.max(error_I_1))
