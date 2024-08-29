# Test derived from 'all_constant_single_ray.py'
# In order to make sure that the solver does not crash if we add points with zero densities in our model,
# this test file will be run using zero densities on a ray
import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'

import numpy             as np
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte


dimension = 1
npoints   = 50
nrays     = 2
nspecs    = 3
nlspecs   = 1
nquads    = 1

nH2  = 0.0                     # [m^-3]
nTT  = 0.0                     # [m^-3]
temp = 4.5E+00                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+12                 # [m]
dv   = 1.5E+02 / magritte.CC   # [fraction of speed of light]


def create_model ():
    """
    Create a model file for the all_constant benchmark, single ray.
    """

    modelName = f'all_zero_single_ray'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'

    #to make sure that we have a region with constant and one with large velocity gradient
    def velocity(i):
        if (i<=npoints/2):
            return 0
        else:
            return (i-npoints/2)*dv

    model = magritte.Model ()
    model.parameters.set_spherical_symmetry(False)
    model.parameters.set_model_name        (modelFile)
    model.parameters.set_dimension         (dimension)
    model.parameters.set_npoints           (npoints)
    model.parameters.set_nrays             (nrays)
    model.parameters.set_nspecs            (nspecs)
    model.parameters.set_nlspecs           (nlspecs)
    model.parameters.set_nquads            (nquads)

    model.geometry.points.position.set([[i*dx, 0, 0] for i in range(npoints)])
    model.geometry.points.velocity.set([[velocity(i), 0, 0] for i in range(npoints)])

    model.chemistry.species.abundance = [[nTT,  nH2, 0.0] for _ in range(npoints)]
    model.chemistry.species.symbol    = ['test', 'H2', 'e-']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    setup.set_Delaunay_neighbor_lists (model)
    setup.set_Delaunay_boundary       (model)
    setup.set_boundary_condition_CMB  (model)
    setup.set_uniform_rays            (model)
    setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model (nosave=False):

    modelName = f'all_zero_single_ray'
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


    timer3 = tools.Timer('shortchar 0  ')
    timer3.start()
    model.compute_radiation_field_shortchar_order_0 ()
    timer3.stop()
    u_0s = np.array(model.radiation.u)
    print(u_0s)

    timer4 = tools.Timer('feautrier 2  ')
    timer4.start()
    model.compute_radiation_field_feautrier_order_2 ()
    timer4.stop()
    u_2f = np.array(model.radiation.u)

    x  = np.array(model.geometry.points.position)[:,0]
    v  = np.array(model.geometry.points.velocity)[:,0]
    nu = np.array(model.radiation.frequencies.nu)

    ld = model.lines.lineProducingSpecies[0].linedata

    k = 0

    frq = ld.frequency[k]
    pop = 0.0
    phi = tools.profile        (ld, k, temp, (turb/magritte.CC)**2, frq)
    eta = 0.0
    chi = 0.0
    src = 0.0
    bdy1 = tools.I_CMB          (frq * (1-v))
    bdyn = tools.I_CMB          (frq * (1+v))

    def I_0 (x):
        return src + (bdy1-src)*np.exp(-chi*x)

    def I_1 (x):
        return src + (bdyn-src)*np.exp(-chi*(x[-1]-x))

    def u_ (x):
        return 0.5 * (I_0(x) + I_1(x))

    #should compute boundary intensity

    error_u_0s = np.abs(tools.relative_error (u_(x), u_0s[0,:,0]))
    error_u_2f = np.abs(tools.relative_error (u_(x), u_2f[0,:,0]))

    result  = f'--- Benchmark name ----------------------------\n'
    result += f'{modelName                                    }\n'
    result += f'--- Parameters --------------------------------\n'
    result += f'dimension = {model.parameters.dimension()     }\n'
    result += f'npoints   = {model.parameters.npoints  ()     }\n'
    result += f'nrays     = {model.parameters.nrays    ()     }\n'
    result += f'nquads    = {model.parameters.nquads   ()     }\n'
    result += f'--- Accuracy ----------------------------------\n'
    result += f'max error in shortchar 0 = {np.max(error_u_0s)}\n'
    result += f'max error in feautrier 2 = {np.max(error_u_2f)}\n'
    result += f'--- Timers ------------------------------------\n'
    result += f'{timer1.print()                               }\n'
    result += f'{timer2.print()                               }\n'
    result += f'{timer3.print()                               }\n'
    result += f'{timer4.print()                               }\n'
    result += f'-----------------------------------------------\n'

    print(result)

    if not nosave:
        with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
            log.write(result)

        fig = plt.figure(dpi=150)
        plt.title(modelName)
        plt.scatter(x, u_0s[0,:,0], s=0.5, label='0s', zorder=1)
        plt.scatter(x, u_2f[0,:,0], s=0.5, label='2f', zorder=1)
        plt.plot(x, u_(x), c='lightgray', zorder=0)
        plt.legend()
        plt.xscale('log')
        plt.xlabel('r [m]')
        plt.ylabel('Mean intensity [W/m$^{2}$]')
        # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)
        plt.show()

    #error bounds are chosen somewhat arbitrarily, based on previously obtained results; this should prevent serious regressions.
    FIRSTORDER_AS_EXPECTED=(np.max(error_u_0s)<6.42e-6)
    FEAUTRIER_AS_EXPECTED=(np.max(error_u_2f)<6.41e-6)

    if not FIRSTORDER_AS_EXPECTED:
        print("First order solver max error too large: ", np.max(error_u_0s))
    if not FEAUTRIER_AS_EXPECTED:
        print("Feautrier solver max error too large: ", np.max(error_u_2f))

    #The main thing we care about is this test not crashing, but at the same time,
    # it would be nice of zero density point do not change the incomping intensity significantly
    return (FEAUTRIER_AS_EXPECTED&FIRSTORDER_AS_EXPECTED)



def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
