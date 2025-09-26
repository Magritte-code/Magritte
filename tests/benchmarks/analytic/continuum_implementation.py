#simple 1D test to check whether adding a continuum background to a model works
# Does not include any line contributions - just a continuum source function

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
import astropy.units  as units


dimension = 1
npoints   = 50
nrays     = 2
nspecs    = 3
nlspecs   = 1
nquads    = 1

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+03                 # [m^-3]
fraction_density_dust = 1.0E-4 # [.] fraction of mass in dust
temp = 4.5E+00                 # [K]
dust_temp = 1.5E+02            # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+12                 # [m]
dv   = 0.0E+00 / magritte.CC   # [fraction of speed of light]


def create_model ():
    """
    Create a model file for the all_constant benchmark, single ray.
    """

    modelName = f'all_constant_single_ray_continuum'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'
    continuumFile = f'{datdir}fe50o_henning1995.txt'
    continuum_ref_density = 4.9 * units.g/(units.cm)**3 # [g/cm3]; will be converted to the correct SI units below
    #corresponds Fe_(0.5)Mg_(0.5)O from Henning 1995, reference density somewhere between 4.79 g/cm3 and 5.05 g/cm3, so Ill assume 4.9 g/cm3


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

    model.chemistry.species.abundance = [[0*nTT, nH2, 0.0] for _ in range(npoints)]#disable line emission by setting density to 0
    #TODO: check if magritte can technically handle a model with only continuum, and no lines present at all
    model.chemistry.species.symbol    = ['test', 'H2', 'e-']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    frequency_grid, complex_refractive_index = tools.read_dust_opacity_table(continuumFile, 1e-6*units.m)
    frequency_grid_Hz, absorption_coeff_per_density = tools.convert_dust_opacity_table_to_SI(frequency_grid, complex_refractive_index, continuum_ref_density)

    model.dust.dust_opacities.set((absorption_coeff_per_density * fraction_density_dust * nH2 * 2.016 / 6.022e23).value[None, :]*np.ones(npoints)[:, None])#H2: 2.016 g/mol, 6.022e23 particles/mol -> particles/m3 to g/m3
    model.dust.dust_frequencies.set((frequency_grid_Hz).value)
    model.dust.dust_temperature.set(dust_temp * np.ones(npoints))

    setup.set_Delaunay_neighbor_lists (model)
    setup.set_Delaunay_boundary       (model)
    setup.set_boundary_condition_CMB  (model)
    setup.set_uniform_rays            (model)
    setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model (nosave=False):

    modelName = f'all_constant_single_ray_continuum'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    magritte.pcmt_set_n_threads_avail(1)

    dust_frequencies = np.array(model.dust.dust_frequencies)
    chi = np.array(model.dust.dust_opacities)[0,:][None, :]#same for all points in this model
    #conveniently, these frequencies correspond to the frequencies used to tabulate the continuum opacity, so no interpolation is needed
    dust_temp = np.array(model.dust.dust_temperature)[0]#same for all points in this model

    timer2 = tools.Timer('setting model')
    timer2.start()
    # model.compute_spectral_discretisation ()
    model.set_custom_spectral_discretization(dust_frequencies)
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()

    print("after setup")

    timer3 = tools.Timer('shortchar 0  ')
    timer3.start()
    model.compute_radiation_field_shortchar_order_0 ()
    timer3.stop()
    u_0s = np.array(model.radiation.u)

    print("after shortchar")

    timer4 = tools.Timer('feautrier 2 ')
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

    x  = np.array(model.geometry.points.position)[:,0][:, None]
    # nu = np.array(model.radiation.frequencies.nu)

    # ld = model.lines.lineProducingSpecies[0].linedata

    # k = 0

    

    src = tools.planck(dust_temp, dust_frequencies)[None, :]#[None, ]


    # frq = ld.frequency[k]
    # pop = tools.LTEpop         (ld, temp) * nTT
    # phi = tools.profile        (ld, k, temp, (turb/magritte.CC)**2, frq)
    # eta = tools.lineEmissivity (ld, pop)[k] * phi
    # chi = tools.lineOpacity    (ld, pop)[k] * phi
    # src = tools.lineSource     (ld, pop)[k]
    bdy = tools.I_CMB          (dust_frequencies)[None, :]

    print("check sizes")
    print("src", np.shape(src))
    print("bdy", np.shape(bdy))
    print("chi", np.shape(chi))

    def I_0 (x):
        return src + (bdy-src)*np.exp(-chi*x)

    def I_1 (x):
        return src + (bdy-src)*np.exp(-chi*(x[-1]-x))

    def u_ (x):
        return 0.5 * (I_0(x) + I_1(x))

    error_u_0s = np.abs(tools.relative_error (u_(x), u_0s[0,:,:]))
    error_u_2f = np.abs(tools.relative_error (u_(x), u_2f[0,:,:]))
    error_I_0_2f_uv = np.abs(tools.relative_error (I_0(x), u_2f_uv[0,:,:]+v_2f_uv[0,:,:]))
    error_I_1_2f_uv = np.abs(tools.relative_error (I_1(x), u_2f_uv[0,:,:]-v_2f_uv[0,:,:]))

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
        plt.plot(x, u_(x), c='lightgray', zorder=0)
        plt.legend()
        plt.xscale('log')
        plt.xlabel('r [m]')
        plt.ylabel('Mean intensity [W/m$^{2}$]')
        plt.show()
        # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

    #returning whether output is as expected (not too far from the input)
    # max_diff=max(u_(x))-min(u_(x))
    # print(max_diff)
    #error should at max be proportional to max diff? only maybe for testing non analytic models

    #error bounds are chosen somewhat arbitrarily, based on previously obtained results; this should prevent serious regressions.
    FEAUTRIER_AS_EXPECTED=(np.max(error_u_2f)<1.7e-4)
    FIRSTORDER_AS_EXPECTED=(np.max(error_u_0s)<1.7e-7)
    FEAUTRIER_UV_AS_EXPECTED=(np.max(error_I_0_2f_uv)<2.0e-4) and (np.max(error_I_1_2f_uv)<2.0e-4)

    if not FIRSTORDER_AS_EXPECTED:
        print("First order solver max error too large: ", np.max(error_u_0s))
    if not FEAUTRIER_AS_EXPECTED:
        print("Feautrier solver max error too large: ", np.max(error_u_2f))
    if not FEAUTRIER_UV_AS_EXPECTED:
        print("Feautrier solver with uv max error too large: ", np.max(error_I_0_2f_uv), np.max(error_I_1_2f_uv))


    return (FEAUTRIER_AS_EXPECTED&FIRSTORDER_AS_EXPECTED&FEAUTRIER_UV_AS_EXPECTED)


def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
