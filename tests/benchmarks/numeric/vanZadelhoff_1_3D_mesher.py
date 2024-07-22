import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'

import numpy             as np
import scipy             as sp
import healpy            as hp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.mesher   as mesher
import magritte.core     as magritte

from scipy.interpolate import interp1d


dimension = 3
nrays     = 12*3**2
nspecs    = 3
nlspecs   = 1
nquads    = 11

r_in   = 1.0E13   # [m]
r_out  = 7.8E16   # [m]
nH2_in = 2.0E13   # [m^-3]
temp   = 20.0     # [K]
turb   = 150.00   # [.]

get_X_mol = {
    'a' : 1.0E-8,
    'b' : 1.0E-6
}

scale_max = 0.05 * r_out
scale_min = 0.05 * r_in
scale_cte = 0.05 * r_in
scale_fun = f'{scale_cte / r_in**2} * (x*x + y*y + z*z)'

meshName = f'{moddir}/vanZadelhoff_1_3D_mesher.vtk'

mesher.create_mesh_from_function(
    meshName       = meshName,
    boundary       = mesher.boundary_sphere_in_sphere(
                         radius_in  = r_in,
                         radius_out = r_out),
    scale_min      = scale_min,
    scale_max      = scale_max,
    scale_function = scale_fun )

mesh = mesher.Mesh(meshName)

npoints = len(mesh.points)
# nbs     = [n for sublist in mesh.neighbors for n in sublist]
# n_nbs   = [len(sublist) for sublist in mesh.neighbors]

rs = np.linalg.norm(mesh.points, axis=1)


def create_model (a_or_b):
    """
    Create a model file for the van Zadelhoff 1 benchmark in 3D.
    """

    modelName = f'vanZadelhoff_1{a_or_b}_3D_mesher'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'

    X_mol = get_X_mol[a_or_b]

    def nH2 (r):
        return nH2_in * np.power(r_in/r, 2.0)

    def nTT (r):
        return X_mol * nH2(r)

    model = magritte.Model ()
    model.parameters.set_spherical_symmetry(False)
    model.parameters.set_model_name        (modelFile)
    model.parameters.set_dimension         (dimension)
    model.parameters.set_npoints           (npoints)
    model.parameters.set_nrays             (nrays)
    model.parameters.set_nspecs            (nspecs)
    model.parameters.set_nlspecs           (nlspecs)
    model.parameters.set_nquads            (nquads)

    model.geometry.points.position.set(mesh.points)
    model.geometry.points.velocity.set(np.zeros((npoints, 3)))

#     model.geometry.points.  neighbors.set(  nbs)
#     model.geometry.points.n_neighbors.set(n_nbs)

    model.chemistry.species.abundance = [[nTT(r), nH2(r), 0.0] for r in rs]
    model.chemistry.species.symbol    = ['test', 'H2', 'e-']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    model = setup.set_Delaunay_neighbor_lists (model)

    model.parameters.set_nboundary(len(mesh.boundary))
    model.geometry.boundary.boundary2point.set(mesh.boundary)

    model = setup.set_boundary_condition_CMB  (model)
    model = setup.set_uniform_rays            (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    model = setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model (a_or_b, nosave=False):

    modelName = f'vanZadelhoff_1{a_or_b}_3D_mesher'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    model.parameters.min_line_opacity = 1.0e-13

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()

    timer3 = tools.Timer('running model')
    timer3.start()
    model.compute_level_populations (True, 100)
    timer3.stop()

    pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), 2))
    abun = np.array(model.chemistry.species.abundance)[:,0]
    rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)

    (i,ra,rb,nh,tk,nm,vr,db,td,lp0,lp1) = np.loadtxt (f'{curdir}/Ratran_results/vanZadelhoff_1{a_or_b}.out', skiprows=14, unpack=True)

    interp_0 = interp1d(0.5*(ra+rb), lp0, fill_value='extrapolate')
    interp_1 = interp1d(0.5*(ra+rb), lp1, fill_value='extrapolate')

    error_0 = np.abs(tools.relative_error(pops[:,0]/abun, interp_0(rs)))
    error_1 = np.abs(tools.relative_error(pops[:,1]/abun, interp_1(rs)))

    result  = f'--- Benchmark name -----------------------\n'
    result += f'{modelName                               }\n'
    result += f'--- Parameters ---------------------------\n'
    result += f'dimension = {model.parameters.dimension()}\n'
    result += f'npoints   = {model.parameters.npoints  ()}\n'
    result += f'nrays     = {model.parameters.nrays    ()}\n'
    result += f'nquads    = {model.parameters.nquads   ()}\n'
    result += f'--- Accuracy -----------------------------\n'
    result += f'max error in (0) = {np.max(error_0[1:])  }\n'
    result += f'max error in (1) = {np.max(error_1[1:])  }\n'
    result += f'--- Timers -------------------------------\n'
    result += f'{timer1.print()                          }\n'
    result += f'{timer2.print()                          }\n'
    result += f'{timer3.print()                          }\n'
    result += f'------------------------------------------\n'

    print(result)

    if not nosave:
        with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
            log.write(result)

        plt.title(modelName)
        plt.scatter(rs, pops[:,0]/abun, s=0.5, label='i=0', zorder=1)
        plt.scatter(rs, pops[:,1]/abun, s=0.5, label='i=1', zorder=1)
        plt.plot(ra, lp0, c='lightgray', zorder=0)
        plt.plot(ra, lp1, c='lightgray', zorder=0)
        plt.legend()
        plt.xscale('log')
        plt.xlabel('r [m]')
        plt.ylabel('fractional level populations [.]')
        plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

    return


def run_test (nosave=False):

    create_model ('a')
    run_model    ('a', nosave)

    create_model ('b')
    run_model    ('b', nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
