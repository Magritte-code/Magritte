import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/data/'
moddir = f'{curdir}/models/'
resdir = f'{curdir}/results/'

import numpy             as np
import scipy             as sp
import healpy            as hp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.mesher   as mesher
import magritte.core     as magritte
# import plotneighbors.py

from scipy.interpolate import interp1d

# dimension = 3
# nrays     = 12*3**2
# nspecs    = 5
# nlspecs   = 1
# nquads    = 11
#
# r_in   = 1.0E13   # [m]
# r_out  = 7.8E16   # [m]
# nH2_in = 2.0E13   # [m^-3]
# temp   = 20.0     # [K]
# turb   = 150.00   # [.]
#
# get_X_mol = {
#     'a' : 1.0E-8,
#     'b' : 1.0E-6
# }
#
# scale_max = 0.05 * r_out
# scale_min = 0.05 * r_in
# scale_cte = 0.05 * r_in
# scale_fun = f'{scale_cte / r_in**2} * (x*x + y*y + z*z)'
#
# meshName = f'{moddir}/vanZadelhoff_1_3D_mesher.vtk'
#
# mesher.create_mesh_from_function(
#     meshName       = meshName,
#     boundary       = mesher.boundary_sphere_in_sphere(
#                          radius_in  = r_in,
#                          radius_out = r_out),
#     scale_min      = scale_min,
#     scale_max      = scale_max,
#     scale_function = scale_fun )
#
# mesh = mesher.Mesh(meshName)
#
# npoints = len(mesh.points)
# nbs     = [n for sublist in mesh.neighbors for n in sublist]
# n_nbs   = [len(sublist) for sublist in mesh.neighbors]
#
# rs = np.linalg.norm(mesh.points, axis=1)


# def create_model (a_or_b):
#     """
#     Create a model file for the van Zadelhoff 1 benchmark in 1D.
#     """
#
#     modelName = f'vanZadelhoff_1{a_or_b}_3D_mesher'
#     modelFile = f'{moddir}{modelName}.hdf5'
#     lamdaFile = f'{datdir}test.txt'
#
#     X_mol = get_X_mol[a_or_b]
#
#     def nH2 (r):
#         return nH2_in * np.power(r_in/r, 2.0)
#
#     def nTT (r):
#         return X_mol * nH2(r)
#
#     model = magritte.Model ()
#     model.parameters.set_spherical_symmetry(False)
#     model.parameters.set_pop_prec          (1.0e-6)
#     model.parameters.set_model_name        (modelFile)
#     model.parameters.set_dimension         (dimension)
#     model.parameters.set_npoints           (npoints)
#     model.parameters.set_nrays             (nrays)
#     model.parameters.set_nspecs            (nspecs)
#     model.parameters.set_nlspecs           (nlspecs)
#     model.parameters.set_nquads            (nquads)
#
#     model.geometry.points.position.set(mesh.points)
#     model.geometry.points.velocity.set(np.zeros((npoints, 3)))
#
#     model.geometry.points.multiscale.set_all_neighbors(n_nbs,nbs)
#     # model.geometry.points.  neighbors.set(  nbs)
#     # model.geometry.points.n_neighbors.set(n_nbs)
#
#     model.chemistry.species.abundance = [[     0.0, nTT(r), nH2(r),  0.0,      1.0] for r in rs]
#     model.chemistry.species.symbol    =  ['dummy0', 'test',   'H2', 'e-', 'dummy1']
#
#     model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
#     model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))
#
#     model.parameters.set_nboundary(len(mesh.boundary))
#     model.geometry.boundary.boundary2point.set(mesh.boundary)
#
#     model = setup.set_boundary_condition_CMB  (model)
#     model = setup.set_uniform_rays            (model)
#     model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
#     model = setup.set_quadrature              (model)
#
#     model.write()
#
#     return #magritte.Model (modelFile)


def run_model (nosave=False):

    modelName = f'model_jan'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    nlevels=2;#should be coarsest level; misleading name
    #2 multigrid levels, minimum 1 point remaining, 0.1 as tolerance, mgImplementation=1 (Naive,Vcycle,Wcycle)
    # model.setup_multigrid(2,nlevels,0.1,1);
    timer2.stop()

    print("point 152833 lies not on boundary? "+str(model.geometry.not_on_boundary(152833)))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # neighbors after coarsening twice
    test_neighbors=[152399,152768,152769,152770,152776,152778,152832,152834,152857,152858,152896,152897,152898,152922,152970,   152833]
    # test_neighbors=model.geometry.points.multiscale.get_neighbors(152833,2);
    for neighbor in test_neighbors:
        ax.scatter(np.array(model.geometry.points.position)[neighbor][0], np.array(model.geometry.points.position)[neighbor][1], np.array(model.geometry.points.position)[neighbor][2])
        ax.text(np.array(model.geometry.points.position)[neighbor][0], np.array(model.geometry.points.position)[neighbor][1], np.array(model.geometry.points.position)[neighbor][2], neighbor, color='black')
    # ax.scatter(model.geometry.points.position[152833].x(), model.geometry.points.position[152833].y(), model.geometry.points.position[152833].z())

    # plt.show()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection='3d')
    # neighbors before coarsening
    test_neighbors=[152768,152769,152770,152777,152832,152834,152840,152841,152897,152898,152905,   152833]
    # test_neighbors=model.geometry.points.multiscale.get_neighbors(152833,2);
    for neighbor in test_neighbors:
        ax2.scatter(np.array(model.geometry.points.position)[neighbor][0], np.array(model.geometry.points.position)[neighbor][1], np.array(model.geometry.points.position)[neighbor][2])
        ax2.text(np.array(model.geometry.points.position)[neighbor][0], np.array(model.geometry.points.position)[neighbor][1], np.array(model.geometry.points.position)[neighbor][2], neighbor, color='black')
    # ax.scatter(model.geometry.points.position[152833].x(), model.geometry.points.position[152833].y(), model.geometry.points.position[152833].z())

    plt.show()

    # timer3 = tools.Timer('running model')
    # timer3.start()
    # model.compute_level_populations(True, 100)
    # # model.compute_level_populations_multigrid(True, 20)
    # timer3.stop()

    # pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), model.lines.lineProducingSpecies[0].linedata.nlev))
    # abun = np.array(model.chemistry.species.abundance)[:,1]
    # rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)
    #
    # # (i,ra,rb,nh,tk,nm,vr,db,td,lp0,lp1) = np.loadtxt (f'{curdir}/Ratran_results/vanZadelhoff_1{a_or_b}.out', skiprows=14, unpack=True)
    #
    # # interp_0 = interp1d(0.5*(ra+rb), lp0, fill_value='extrapolate')
    # # interp_1 = interp1d(0.5*(ra+rb), lp1, fill_value='extrapolate')
    # #
    # # error_0 = tools.relative_error(pops[:,0]/abun, interp_0(rs))
    # # error_1 = tools.relative_error(pops[:,1]/abun, interp_1(rs))
    #
    # result  = f'--- Benchmark name -----------------------\n'
    # result += f'{modelName                               }\n'
    # result += f'--- Parameters ---------------------------\n'
    # result += f'dimension = {model.parameters.dimension()}\n'
    # result += f'npoints   = {model.parameters.npoints  ()}\n'
    # result += f'nrays     = {model.parameters.nrays    ()}\n'
    # result += f'nquads    = {model.parameters.nquads   ()}\n'
    # # result += f'--- Accuracy -----------------------------\n'
    # # result += f'max error in (0) = {np.max(error_0[1:])  }\n'
    # # result += f'max error in (1) = {np.max(error_1[1:])  }\n'
    # result += f'--- Timers -------------------------------\n'
    # result += f'{timer1.print()                          }\n'
    # result += f'{timer2.print()                          }\n'
    # result += f'{timer3.print()                          }\n'
    # result += f'------------------------------------------\n'
    #
    # print(result)
    #
    # # chosen_lines=np.arange(0,model.lines.lineProducingSpecies[0].linedata.nlev,5)
    # chosen_lines=[0,1,2,3,4]#first five lines are the most important
    #
    # if not nosave:
    #     with open(f'{resdir}{modelName}-{timestamp}_non_multigrid.log' ,'w') as log:
    #         log.write(result)
    #
    #     plt.title(modelName)
    #     for levidx in range(len(chosen_lines)):
    #         plt.scatter(rs, pops[:,chosen_lines[levidx]]/abun, s=0.5, label=str('i='+str(chosen_lines[levidx])), zorder=1)
    #     # plt.scatter(rs, pops[:,1]/abun, s=0.5, label='i=1', zorder=1)
    #     # plt.plot(ra, lp0, c='lightgray', zorder=0)
    #     # plt.plot(ra, lp1, c='lightgray', zorder=0)
    #     plt.legend()
    #     plt.xscale('log')
    #     plt.xlabel('r [m]')
    #     plt.ylabel('fractional level populations [.]')
    #     plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

    return


def run_test (nosave=False):
# for simplicity, we only try the first model
    # create_model ('a')
    run_model    (nosave)

    # create_model ('b')
    # run_model    ('b', nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)