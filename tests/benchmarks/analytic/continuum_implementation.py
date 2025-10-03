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
import scipy.integrate


dimension = 1
npoints   = 50
nrays     = 2
nspecs    = 3
nlspecs   = 1
nquads    = 1

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+03                 # [m^-3]
fraction_density_dust = 1.0E-5 # [.] fraction of mass in dust
temp = 4.5E+00                 # [K]
dust_temp = 1.5E+02            # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+12                 # [m]
dv   = 0.0E+03 / magritte.CC   # [fraction of speed of light]


def create_model ():
    """
    Create a model file for the all_constant benchmark, single ray.
    """

    modelName = f'all_constant_single_ray_continuum'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'
    continuumFile = f'{datdir}fe50o_henning1995.txt'
    continuum_ref_density = 4.9 * units.g/(units.cm)**3 # [g/cm3]; will be converted to the correct SI units below
    #corresponds Fe_(0.5)Mg_(0.5)O from Henning 1995, reference density somewhere between 4.79 g/cm3 and 5.05 g/cm3, so I'll assume 4.9 g/cm3


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

    #Note: Magritte currently cannot handle creating continuum-only models, so we include a dummy species with zero abundance to avoid issues
    model.chemistry.species.abundance = [[0*nTT, nH2, 0.0] for _ in range(npoints)]#disable line emission/absorption by setting density to 0
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
    positions = np.array(model.geometry.points.position)

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.set_custom_spectral_discretization(dust_frequencies)
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()

    #compute the image, to check whether the intensity is as expected
    timer3 = tools.Timer('Compute intensity image')
    timer3.start()
    model.compute_image_new(0,1,1)
    timer3.stop()

    timer4 = tools.Timer('Compute image optical depth')
    timer4.start()
    model.compute_image_optical_depth_new(0,1,1)
    timer4.stop()

    I_image = np.array(model.images[0].I)[0,:]
    tau_image = np.array(model.images[1].I)[0,:]

    def evaluate_dust_opacity(frequency, velocity):
        return np.interp(dust_frequencies[None, :]*(1-velocity[:, None]), frequency, chi[0,:])

    
    #compute optical depth
    chi_shifted = evaluate_dust_opacity(dust_frequencies, np.array(model.geometry.points.velocity)[:,0])
    tau_ref = np.trapz(chi_shifted, x=np.array(model.geometry.points.position)[:,0], axis=0)
    #TODO: fix the calculation of the reference intensity; currently it is wrong for any model with non-zero velocity field due to not taking into account the doppler shifts properly
    #test attempt below
    # def evaluate_source_function(frequency, velocity):
    #     return tools.planck(dust_temp[None, None], dust_frequencies[None, :]*(np.ones((1,1))-velocity[:, None]))
    # S_shifted = evaluate_source_function(dust_frequencies, np.array(model.geometry.points.velocity)[:,0])
    # cumsum_tau = np.zeros((chi_shifted.shape[0]+1, chi_shifted.shape[1]))
    # cumsum_tau[1:-1,:] = scipy.integrate.cumulative_trapezoid(chi_shifted, axis=0)*dx
    # cumsum_tau[-1,:] = tau_ref
    # intensity_contributions = S_shifted * (1 - np.exp(-dx*chi_shifted)) * np.exp(-cumsum_tau[:-1, :])
    # ref_intensity = np.sum(intensity_contributions, axis=0) + tools.I_CMB(dust_frequencies)*np.exp(-tau_ref)

    # Technically, the source function for the dust is not constant due to doppler shifts, therefore NO doppler shifts in this benchmark
    ref_intensity = tools.planck(dust_temp, dust_frequencies) * (1-np.exp(-tau_ref)) + tools.I_CMB(dust_frequencies)*np.exp(-tau_ref)

    reldiff_I = 2.0*np.abs((I_image - ref_intensity)/(ref_intensity+I_image))
    reldiff_tau = 2.0*np.abs((tau_image - tau_ref)/(tau_ref+tau_image))

    result  = f'--- Benchmark name ----------------------------\n'
    result += f'{modelName                                    }\n'
    result += f'--- Parameters --------------------------------\n'
    result += f'dimension = {model.parameters.dimension()     }\n'
    result += f'npoints   = {model.parameters.npoints  ()     }\n'
    result += f'nrays     = {model.parameters.nrays    ()     }\n'
    result += f'nquads    = {model.parameters.nquads   ()     }\n'
    result += f'--- Accuracy ----------------------------------\n'
    result += f'max error in I = {np.max(reldiff_I[:40])}      \n'
    result += f'max error in tau = {np.max(reldiff_tau)}       \n'
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

        plt.figure(dpi = 150)
        plt.title('Intensity')
        plt.plot(ref_intensity, label='ref')
        plt.plot(I_image, label='model')
        plt.axvline(x=40, color='gray', linestyle='--', label='end of checked range')
        plt.yscale('log')
        plt.legend()
        plt.savefig(f'{resdir}{modelName}-intensity-{timestamp}.png', dpi=150)
        plt.figure(dpi=150)
        plt.title('Optical depth')
        plt.plot(tau_ref, label='ref')
        plt.plot(tau_image, label='model')
        plt.yscale('log')
        plt.savefig(f'{resdir}{modelName}-optical_depth-{timestamp}.png', dpi=150)

    #error bounds are chosen somewhat arbitrarily, based on previously obtained results; this should prevent serious regressions.
    RELDIFF_I_AS_EXPECTED = (np.max(reldiff_I[:40])<1.7e-5)#well, the higher frequencies have very low intensities/sharp declines of the source function, so the relative error can be larger there
    RELDIFF_TAU_AS_EXPECTED = (np.max(reldiff_tau)<4.4e-11)

    if not RELDIFF_I_AS_EXPECTED:
        print("Continuum intensity max error too large: ", np.max(np.max(reldiff_I[:40])))
    if not RELDIFF_TAU_AS_EXPECTED:
        print("Continuum optical depth max error too large: ", np.max(np.max(reldiff_tau)))


    return (RELDIFF_I_AS_EXPECTED&RELDIFF_TAU_AS_EXPECTED)


def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
