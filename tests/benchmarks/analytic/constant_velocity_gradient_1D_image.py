import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'

import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import ipywidgets        as widgets
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte


dimension = 1
npoints   = 100
nrays     = 4
nspecs    = 5
nlspecs   = 1
nquads    = 100

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+08                 # [m^-3]
temp = 4.5E+01                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+04                 # [m]
dv   = 2.5E+02 / magritte.CC   # [fraction of speed of light]

L    = dx*npoints
vmax = dv*npoints


def create_model ():
    """
    Create a model file for the constant velocity gradient benchmark 1D.
    """

    modelName = f'constant_velocity_gradient_1D_image'
    modelFile = f'{moddir}{modelName}.hdf5'
    lamdaFile = f'{datdir}test.txt'

    model = magritte.Model ()
    model.parameters.set_spherical_symmetry(True)
    model.parameters.set_model_name        (modelFile)
    model.parameters.set_dimension         (dimension)
    model.parameters.set_npoints           (npoints)
    model.parameters.set_nrays             (nrays)
    model.parameters.set_nspecs            (nspecs)
    model.parameters.set_nlspecs           (nlspecs)
    model.parameters.set_nquads            (nquads)

    model.geometry.points.position.set([[(i+1)*dx, 0, 0] for i in range(npoints)])
    model.geometry.points.velocity.set([[(i+1)*dv, 0, 0] for i in range(npoints)])

    model.chemistry.species.abundance = [[     0.0,    nTT,  nH2,  0.0,      1.0] for _ in range(npoints)]
    model.chemistry.species.symbol    =  ['dummy0', 'test', 'H2', 'e-', 'dummy1']

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


def run_model (nosave=False, benchindex=0):

    modelName = f'constant_velocity_gradient_1D_image'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    #by default, benchmark uses default magritte settings for max_width_fraction; denotes when to switch over to more expensive, more accurate way of computing optical depths when encountering doppler shifts
    model.parameters.max_width_fraction=0.5;
    #for controlling accuracy of optical depth computation during benchmark runs
    if (benchindex==1):
        model.parameters.max_width_fraction=0.5;#fast, but slightly inaccurate
    elif (benchindex==2):
        model.parameters.max_width_fraction=0.35;#more reasonably accurate, but minor performance hit
    elif (benchindex==3):
        model.parameters.max_width_fraction=0.10;#major performance hit, very accurate

    fcen = model.lines.lineProducingSpecies[0].linedata.frequency[0]
    dd = 1.0e+4 / magritte.CC
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
    pop = tools.LTEpop         (ld, temp) * nTT
    eta = tools.lineEmissivity (ld, pop)[k]
    chi = tools.lineOpacity    (ld, pop)[k]
    src = tools.lineSource     (ld, pop)[k]
    dnu = tools.dnu            (ld, k, temp, (turb/magritte.CC)**2)

    r_in  = rs[ 0]
    r_out = rs[-1]

    def bdy (nu):
        return tools.I_CMB (nu)

    def z_max(r, theta):
        if (theta < np.arcsin(r_in/r)):
            return r * np.cos(theta) - np.sqrt(r_in**2 - (r*np.sin(theta))**2)
        else:
            return r * np.cos(theta) + np.sqrt(L   **2 - (r*np.sin(theta))**2)

    def tau(nu, r, theta):
        l   = float (z_max(r, theta))
        arg = float ((nu - frq) / dnu)
        fct = float (vmax * nu / dnu)
        tauup=chi*L / (fct*dnu) * 0.5 * (sp.special.erf(arg) - sp.special.erf(arg-fct*l/L))
        taudown=chi*L / (-fct*dnu) * 0.5 * (sp.special.erf(arg) - sp.special.erf(arg+fct*l/L))
        return tauup+taudown

    def I_im (nu, r):
        nu = nu + vmax/magritte.CC*frq
        return src + (bdy(nu)-src)*np.exp(-tau(nu, r, 0.5*np.pi))

    im_a = np.array([[I_im(f,r) for f in nu] for r in rs])
    fs   = (nu-frq)/frq*magritte.CC

    def plot (p):
        plt.figure(dpi=150)
        plt.plot(fs, im  [p,:], marker='.')
        plt.plot(fs, im_a[p,:])
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

    FEAUTRIER_AS_EXPECTED=True
    #if we are actually doing benchmarks, the accuracy will depend on how accurately we compute the optical depth
    if (benchindex==1):
        FEAUTRIER_AS_EXPECTED=(np.mean(error_u_2f)<1.8e-5)
    elif(benchindex==2):
        FEAUTRIER_AS_EXPECTED=(np.mean(error_u_2f)<8e-6)
    elif(benchindex==3):
        FEAUTRIER_AS_EXPECTED=(np.mean(error_u_2f)<2e-8)

    if not FEAUTRIER_AS_EXPECTED:
        print("Feautrier solver mean error too large: ", np.mean(error_u_2f), "bench nr: ", benchindex)

    return (FEAUTRIER_AS_EXPECTED)


def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
