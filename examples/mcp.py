import sys
magritteFolder = '/home/frederik/Dropbox/GitHub/Magritte/'
sys.path.append(magritteFolder)

datdir = f'{magritteFolder}/tests/data/'
moddir = f'{magritteFolder}/tests/models/'
resdir = f'{magritteFolder}/tests/results/'

import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte

from astropy import constants

# Path to MCP input file
mcp_file = '/home/frederik/IvS/Sofia/IRAS_22035-MCPinput.txt'

# Load the MCP input file
(Rs, T, nH2, X, Vr, Vt) = np.loadtxt(mcp_file, skiprows=15, unpack=True, usecols=[1,2,4,5,6,7])

# Convert units to SI (where necessary)
Rs  *= 1.0e-2   # convert cm   to m
nH2 *= 1.0e+6   # convert cm-3 to m-3
Vr  *= 1.0e+3   # convert km/s to m/s
Vt  *= 1.0e+3   # convert km/s to m/s

# Convert velocities to fractions of the speed of light
Vr /= constants.c.si.value
Vt /= constants.c.si.value

# Define model file and line data file
modelFile = f'{moddir}/example_mcp.hdf5'
lamdaFile = f'{datdir}/co.txt'

# Define model parameters
dimension = 1         # effective spatial dimension (1 for spherical symmetry)
npoints   = len(Rs)   # number of spatial points
nrays     = 200       # number of rays to trace from each point
nspecs    = 5         # number of chemical species (minimum 5)
nlspecs   = 1         # number of chemical species we consider in line RT
nquads    = 11        # number of roots/weights to use in Gauss-Hermite quadrature


def create_model():
    # Setup magritte model
    model = magritte.Model ()
    model.parameters.set_model_name         (modelFile)
    model.parameters.set_dimension          (dimension)
    model.parameters.set_spherical_symmetry (True)
    model.parameters.set_npoints            (npoints)
    model.parameters.set_nrays              (nrays)
    model.parameters.set_nspecs             (nspecs)
    model.parameters.set_nlspecs            (nlspecs)
    model.parameters.set_nquads             (nquads)
    model.parameters.set_pop_prec           (1.0e-6)

    model.geometry.points.position.set([[r, 0, 0] for r in Rs])
    model.geometry.points.velocity.set([[v, 0, 0] for v in Vr])

    model.chemistry.species.symbol    =  ['dummy0', 'CO', 'H2', 'e-', 'dummy1']
    model.chemistry.species.abundance = [[     0.0, x*h2,   h2,  0.0,      1.0] for (x,h2) in zip(X,nH2)]

    model.thermodynamics.temperature.gas  .set(T    )
    model.thermodynamics.turbulence.vturb2.set(Vt**2)

    model = setup.set_Delaunay_neighbor_lists  (model)
    model = setup.set_Delaunay_boundary        (model)
    model = setup.set_boundary_condition_CMB   (model)
    model = setup.set_rays_spherical_symmetry  (model)
    model = setup.set_linedata_from_LAMDA_file (model, lamdaFile)
    model = setup.set_quadrature               (model)

    model.write()

    return


def run_model():
    # Load magritte model
    model = magritte.Model (modelFile)

    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()

    model.compute_level_populations (True, 100)


if (__name__ == '__main__'):
    create_model ()
    run_model    ()
