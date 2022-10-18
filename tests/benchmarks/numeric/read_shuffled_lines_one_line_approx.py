#by computing the radiation field in the OneLine approximation, we make sure that the indices of the frequencies and the lines match
#Made test by simplifying vanZadelhoff 1a model a lot
import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'

import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte

from scipy.interpolate import interp1d


dimension = 1
npoints   = 25
nrays     = 20
nspecs    = 5
nlspecs   = 1
nquads    = 11

r_in   = 1.0E13   # [m]
r_out  = 7.8E16   # [m]
nH2_in = 2.0E13   # [m^-3]
temp   =  20.00   # [K]
turb   = 150.00   # [.]
T_bound= 1e-8

get_X_mol = 1.0E-8

rs = np.logspace (np.log10(r_in), np.log10(r_out), npoints, endpoint=True)


def create_model():
    """
    Create a model file for the van Zadelhoff 1 benchmark in 1D.
    """

    modelName = f'read_shuffled_lines_1D'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test_unsorted.txt'

    X_mol = get_X_mol

    def nH2 (r):
        return nH2_in * np.power(r_in/r, 2.0)

    def nTT (r):
        return X_mol  * nH2(r)

    model = magritte.Model ()
    model.parameters.set_spherical_symmetry(True)
    model.parameters.set_model_name        (modelFile)
    model.parameters.set_dimension         (dimension)
    model.parameters.set_npoints           (npoints)
    model.parameters.set_nrays             (nrays)
    model.parameters.set_nspecs            (nspecs)
    model.parameters.set_nlspecs           (nlspecs)
    model.parameters.set_nquads            (nquads)

    model.geometry.points.position.set([[r, 0, 0] for r in rs])
    model.geometry.points.velocity.set(np.zeros((npoints, 3)))

    model.chemistry.species.abundance = [[     0.0, nTT(r), nH2(r),  0.0,      1.0] for r in rs]
    model.chemistry.species.symbol    =  ['dummy0', 'test',   'H2', 'e-', 'dummy1']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    model = setup.set_Delaunay_neighbor_lists (model)
    model = setup.set_Delaunay_boundary       (model)
    #boundary intensity should be approx 0 for the frequencies we care about
    model = setup.set_boundary_condition_1D   (model, T_bound, T_bound)
    model = setup.set_rays_spherical_symmetry (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    model = setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model ():

    modelName = f'read_shuffled_lines_1D'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    #for sake of consistency, the timers are included. However, they are not useful in the slightest
    # as this test is too small
    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()

    timer3 = tools.Timer('computing intensities without approx')
    timer3.start()
    model.compute_radiation_field_feautrier_order_2()
    timer3.stop()

    Jtrue = np.array(model.radiation.J)

    model.parameters.one_line_approximation = True

    timer4 = tools.Timer('computing intensities with approx')
    timer4.start()
    model.compute_radiation_field_feautrier_order_2()
    timer4.stop()
    Japprox = np.array(model.radiation.J)

    #Testing whether these results are close enough; as the lines are far enough away from eachother, they should give exactly the same result
    rel_error = np.abs((Japprox-Jtrue)/Jtrue)

    if (len(rel_error[rel_error>1.0e-10])>0):
        print("Difference between OneLine approximation and full computation is too much")
        return False

    return True


def run_test (n):

    create_model ()
    run_model    ()

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test ()
