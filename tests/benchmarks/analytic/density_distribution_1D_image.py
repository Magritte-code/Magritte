import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'

import numpy             as np
import matplotlib.pyplot as plt
import ipywidgets        as widgets
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte


dimension = 1
npoints   = 200
nrays     = 50
nspecs    = 3
nlspecs   = 1
nquads    = 50

r_in   = 1.0E13   # [m]
r_out  = 7.8E16   # [m]
nH2_in = 2.0E13   # [m^-3]
temp   = 2.0E+01  # [K]
turb   = 1.5E+02  # [m/s]

get_X_mol = {
    'a' : 1.0E-8,
    'b' : 1.0E-6
}

rs = np.logspace (np.log10(r_in), np.log10(r_out), npoints, endpoint=True)


def create_model (a_or_b):
    """
    Create a model file for the density distribution benchmark 1D.
    """

    modelName = f'density_distribution_VZ{a_or_b}_1D_image'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'

    X_mol = get_X_mol[a_or_b]

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
    model.geometry.points.velocity.set([[0, 0, 0] for i in range(npoints)])

    model.chemistry.species.abundance = [[nTT(r), nH2(r), 0.0] for r in rs]
    model.chemistry.species.symbol    = ['test', 'H2', 'e-']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    model = setup.set_Delaunay_neighbor_lists (model)
    model = setup.set_Delaunay_boundary       (model)
    model = setup.set_boundary_condition_CMB  (model)
    model = setup.set_rays_spherical_symmetry (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    model = setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model (a_or_b, nosave=False, use_widgets=True):

    modelName = f'density_distribution_VZ{a_or_b}_1D_image'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    X_mol = get_X_mol[a_or_b]

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    model.parameters.min_line_opacity = 1.0e-13

    fcen = model.lines.lineProducingSpecies[0].linedata.frequency[0]
    dd = 3.0e+3 / magritte.CC
    fmin = fcen - fcen*dd
    fmax = fcen + fcen*dd

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation (fmin, fmax)
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()

    npoints = model.parameters.npoints()
    nrays   = model.parameters.nrays  ()
    hnrays  = model.parameters.hnrays ()
    nfreqs  = model.parameters.nfreqs ()

    ray_nr = hnrays - 1

    timer3 = tools.Timer('image        ')
    timer3.start()
    model.compute_image (ray_nr)
    timer3.stop()

    nu = np.array(model.radiation.frequencies.nu)[0]
    rs = np.array(model.images[-1].ImX, copy=True)
    im = np.array(model.images[-1].I,   copy=True)

    ld = model.lines.lineProducingSpecies[0].linedata

    k = 0

    frq = ld.frequency[k]
    pop = tools.LTEpop         (ld, temp) * X_mol * nH2_in
    eta = tools.lineEmissivity (ld, pop)[k]
    chi = tools.lineOpacity    (ld, pop)[k]
    src = tools.lineSource     (ld, pop)[k]

    def phi (nu):
        return tools.profile (ld, k, temp, (turb/magritte.CC)**2, nu)

    def bdy (nu):
        return tools.I_CMB (nu)

    def tau (nu, r, theta):
        pref = chi * phi(nu) * r_in**2
        if (theta == 0.0):
            return pref * (2.0/r_in - 1.0/r - 1.0/r_out)
        if (theta == np.pi):
            return pref * (1.0/r - 1.0/r_out)
        else:
            factor = np.arccos(r*np.sin(theta)/r_out) + 0.5*np.pi - theta
            # if (theta < np.arcsin(r_in/r)):
                # factor = factor - np.arccos(r*np.sin(theta)/r_in)
            return pref / (r*np.sin(theta)) * factor

    def I_im (nu, r):
        return src + (bdy(nu)-src)*np.exp(-2.0*tau(nu, r, 0.5*np.pi))

    im_a = np.array([[I_im(f,r) for f in nu] for r in rs])
    fs   = (nu-frq)/frq*magritte.CC

    def plot (p):
        plt.figure(dpi=150)
        plt.plot(fs,   (    im  [p,:]-tools.I_CMB(nu)), marker='.', label='Magritte')
        # plt.plot(fs,   (0.5*im  [p,:]-2*tools.I_CMB(nu)), marker='.', label='Magritte')
        plt.plot(fs,   (im_a[p,:]-tools.I_CMB(nu)),             label='analytic')
        # plt.plot(fs,   (im_a[p,:]-tools.I_CMB(nu)),             label='analytic')
        plt.plot(fs,           tools.I_CMB(nu),                 label='CMB')
        plt.yscale('log')
        plt.legend()
        
    #during automated testing, the widgets only consume time to create
    if use_widgets:
        widgets.interact(plot, p=(0,npoints-1,1))

    error = np.abs(tools.relative_error(im, im_a))[:-1]

    log_err_min = np.log10(np.min([error]))
    log_err_max = np.log10(np.max([error]))

    bins = np.logspace(log_err_min, log_err_max, 100)

    result  = f'--- Benchmark name ------------------------------\n'
    result += f'{modelName                                      }\n'
    result += f'--- Parameters ----------------------------------\n'
    result += f'dimension = {model.parameters.dimension()       }\n'
    result += f'npoints   = {model.parameters.npoints  ()       }\n'
    result += f'nrays     = {model.parameters.nrays    ()       }\n'
    result += f'nquads    = {model.parameters.nquads   ()       }\n'
    result += f'--- Accuracy ------------------------------------\n'
    result += f'max error in imager = {np.max(error)            }\n'
    result += f'--- Timers --------------------------------------\n'
    result += f'{timer1.print()                                 }\n'
    result += f'{timer2.print()                                 }\n'
    result += f'{timer3.print()                                 }\n'
    result += f'-------------------------------------------------\n'

    print(result)

    if not nosave:
        with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
            log.write(result)

        plt.figure()
        plt.title(modelName)
        plt.hist(error.ravel(), bins=bins, histtype='step')
        plt.xscale('log')
        plt.savefig(f'{resdir}{modelName}_hist-{timestamp}.png', dpi=150)

        plt.figure()
        plt.hist(error.ravel(), bins=bins, histtype='step', cumulative=True)
        plt.hist(error.ravel(), bins=bins, histtype='step', cumulative=True)
        plt.xscale('log')
        plt.savefig(f'{resdir}{modelName}_cumu-{timestamp}.png', dpi=150)

    #setting 'b' not yet used for testing
    if a_or_b == 'a':
        #error bounds are chosen somewhat arbitrarily, based on previously obtained results; this should prevent serious regressions.
        FEAUTRIER_AS_EXPECTED=(np.mean(error)<2e-3)
        if not FEAUTRIER_AS_EXPECTED:
            print("Feautrier solver mean error too large: ", np.mean(error))

        return (FEAUTRIER_AS_EXPECTED)

    return


def run_test (nosave=False):

    create_model ('a')
    run_model    ('a', nosave)

    # create_model ('b')
    # run_model    ('b', nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
