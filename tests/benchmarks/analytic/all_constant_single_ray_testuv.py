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

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+03                 # [m^-3]
temp = 4.5E+00                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+12                 # [m]
dv   = 0.0E+00 / magritte.CC   # [fraction of speed of light]


def create_model ():
    """
    Create a model file for the all_constant benchmark, single ray.
    """

    modelName = f'all_constant_single_ray_testuv'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'


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
    model.geometry.points.velocity.set([[i*dv, 0, 0] for i in range(npoints)])

    model.chemistry.species.abundance = [[nTT, nH2, 0.0] for _ in range(npoints)]
    model.chemistry.species.symbol    = ['test', 'H2', 'e-']

    # model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    #Changing the temperature such that the source function is not constant
    model.thermodynamics.temperature.gas  .set( temp * (1.0+10*np.arange(npoints)/npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    model = setup.set_Delaunay_neighbor_lists (model)
    model = setup.set_Delaunay_boundary       (model)
    model = setup.set_boundary_condition_CMB  (model)
    model = setup.set_uniform_rays            (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    model = setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model (nosave=False):

    modelName = f'all_constant_single_ray_testuv'
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
    timer2.stop()

    timer3 = tools.Timer('shortchar 0  ')
    timer3.start()
    model.compute_radiation_field_shortchar_order_0 ()
    timer3.stop()
    u_0s = np.array(model.radiation.u)
    I_0s = np.array(model.radiation.I)

    timer4 = tools.Timer('feautrier 2  ')
    timer4.start()
    model.compute_radiation_field_feautrier_order_2 ()
    timer4.stop()
    u_2f = np.array(model.radiation.u)

    timer5 = tools.Timer('feautrier 2 uv')
    timer5.start()
    model.compute_radiation_field_feautrier_order_2_uv ()
    timer5.stop()
    u_2f_uv = np.array(model.radiation.u)
    v_2f_uv = np.array(model.radiation.v)

    x  = np.array(model.geometry.points.position)[:,0]
    nu = np.array(model.radiation.frequencies.nu)

    ld = model.lines.lineProducingSpecies[0].linedata

    k = 0

    frq = ld.frequency[k]
    pop = tools.LTEpop         (ld, temp) * nTT
    phi = tools.profile        (ld, k, temp, (turb/magritte.CC)**2, frq)
    eta = tools.lineEmissivity (ld, pop)[k] * phi
    chi = tools.lineOpacity    (ld, pop)[k] * phi
    src = tools.lineSource     (ld, pop)[k]
    bdy = tools.I_CMB          (frq)

    # def I_0 (x):
    #     return src + (bdy-src)*np.exp(-chi*x)

    # def I_1 (x):
    #     return src + (bdy-src)*np.exp(-chi*(x[-1]-x))

    # def u_ (x):
    #     return 0.5 * (I_0(x) + I_1(x))

    # error_u_0s = np.abs(tools.relative_error (u_(x), u_0s[0,:,0]))
    # error_u_2f = np.abs(tools.relative_error (u_(x), u_2f[0,:,0]))
    error_I_0_2f_uv = np.abs(tools.relative_error (I_0s[0,:,0], u_2f_uv[0,:,0]+v_2f_uv[0,:,0]))
    error_I_1_2f_uv = np.abs(tools.relative_error (I_0s[1,:,0], u_2f_uv[0,:,0]-v_2f_uv[0,:,0]))

    result  = f'--- Benchmark name ----------------------------\n'
    result += f'{modelName                                    }\n'
    result += f'--- Parameters --------------------------------\n'
    result += f'dimension = {model.parameters.dimension()     }\n'
    result += f'npoints   = {model.parameters.npoints  ()     }\n'
    result += f'nrays     = {model.parameters.nrays    ()     }\n'
    result += f'nquads    = {model.parameters.nquads   ()     }\n'
    result += f'--- Accuracy ----------------------------------\n'
    result += f'max error in I_0 2f uv = {np.max(error_I_0_2f_uv)}\n'
    result += f'max error in I_1 2f uv = {np.max(error_I_1_2f_uv)}\n'
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
        # plt.plot(x, u_(x), c='lightgray', zorder=0)
        plt.legend()
        # plt.xscale('log')
        plt.xlabel('r [m]')
        plt.ylabel('Mean intensity [W/m$^{2}$]')
        # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)
        plt.figure()
        plt.plot(x, I_0s[0,:,0], label='I shortchar 0')
        plt.plot(x, u_2f_uv[0,:,0]+v_2f_uv[0,:,0], label='I feautrier uv')
        plt.legend()
        plt.figure()
        plt.plot(x, I_0s[1,:,0], label='I shortchar 1')
        plt.plot(x, u_2f_uv[0,:,0]-v_2f_uv[0,:,0], label='I feautrier uv')
        plt.legend()
        plt.figure()
        plt.plot(x, I_0s[0,:,0]-(u_2f_uv[0,:,0]+v_2f_uv[0,:,0]), label='diff')
        plt.figure()
        plt.plot(x, (I_0s[0,:,0]-I_0s[1,:,0])/2.0, label='shortchar diff')
        plt.plot(x, v_2f_uv[0,:,0], label='feautrier v')
        plt.legend()
        plt.show()

    #error bounds are chosen somewhat arbitrarily, based on previously obtained results; this should prevent serious regressions.
    FEAUTRIER_UV_AS_EXPECTED=(np.mean(error_I_0_2f_uv)<0.021) and (np.mean(error_I_1_2f_uv)<0.01)

    if not FEAUTRIER_UV_AS_EXPECTED:
        print("Feautrier solver with uv mean error too large: ", np.mean(error_I_0_2f_uv), np.mean(error_I_1_2f_uv))


    return (FEAUTRIER_UV_AS_EXPECTED)


def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
