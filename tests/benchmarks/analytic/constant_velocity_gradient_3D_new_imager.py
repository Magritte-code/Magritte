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
import ipywidgets        as widgets
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte


dimension = 3
nrays     = 12* 2**0
nshells   = 10

nspecs    = 5
nlspecs   = 1
nquads    = 10

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+08                 # [m^-3]
temp = 4.5E+01                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+04                 # [m]
dv   = 2.5E+02 / magritte.CC   # [fraction of speed of light]

L    = dx*nshells
vmax = dv*nshells

r_shell = np.linspace (dx, nshells*dx, nshells, endpoint=True)

npoints_in_shell = [hp.nside2npix(2+s) for s in range(nshells)]
npoints          = sum(npoints_in_shell)

xyz = np.empty(shape=[0, 3])
velocity = np.empty(shape=[0, 3])
corresp_r= np.empty(shape=[0])


for (r, n) in zip(r_shell, npoints_in_shell):
    pos = r*np.array(hp.pixelfunc.pix2vec(hp.npix2nside(n), range(n))).T
    pos = sp.spatial.transform.Rotation.random().apply(pos)
    xyz = np.concatenate((xyz, pos))

    i = r / dx
    vel = i * dv * pos / r
    #velocity times direction (unit vector)
    velocity = np.concatenate((velocity, vel))
    corresp_r = np.append(corresp_r, r)


def create_model ():
    """
    Create a model file for the constant velocity gradient benchmark 1D.
    """

    modelName = f'constant_velocity_gradient_3D_image_new'
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

    model.geometry.points.position.set(xyz)
    model.geometry.points.velocity.set(velocity)

    model.chemistry.species.abundance = [[     0.0,    nTT,  nH2,  0.0,      1.0] for _ in range(npoints)]
    model.chemistry.species.symbol    =  ['dummy0', 'test', 'H2', 'e-', 'dummy1']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    model = setup.set_Delaunay_neighbor_lists (model)

    inner_bdy = np.array([b for b in range(npoints_in_shell[0])])
    outer_bdy = np.array([b for b in range(npoints-npoints_in_shell[-1], npoints)])
    boundary2point = np.concatenate((inner_bdy, outer_bdy))

    model.geometry.boundary.boundary2point.set(boundary2point)
    model.parameters.set_nboundary(len(boundary2point))

    model = setup.set_boundary_condition_CMB  (model)
    model = setup.set_uniform_rays (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    model = setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model (nosave=False, use_widgets=True):

    modelName = f'constant_velocity_gradient_3D_image_new'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    #by default, benchmark uses default magritte settings for max_width_fraction; denotes when to switch over to more expensive, more accurate way of computing optical depths when encountering doppler shift
    model.parameters.max_width_fraction=0.5;
    #however, this benchmark shows how difficult it is for the imager to deal with non-negligible impact of the boundary points. Thus it is not very accurate.

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
    model.compute_image_new (ray_nr, 50, 50)
    timer3.stop()

    nu = np.array(model.radiation.frequencies.nu)[0]
    positions = np.array(model.geometry.points.position)
    im = np.array(model.images[-1].I,   copy=True)
    imx = np.array(model.images[-1].ImX, copy=True)
    imy = np.array(model.images[-1].ImY, copy=True)
    rs = np.sqrt(np.power(imx, 2) + np.power(imy, 2))

    ld = model.lines.lineProducingSpecies[0].linedata

    k = 0

    frq = ld.frequency[k]
    pop = tools.LTEpop         (ld, temp) * nTT
    eta = tools.lineEmissivity (ld, pop)[k]
    chi = tools.lineOpacity    (ld, pop)[k]
    src = tools.lineSource     (ld, pop)[k]
    dnu = tools.dnu            (ld, k, temp, (turb/magritte.CC)**2)

    r_in  = np.linalg.norm(positions[0])
    r_out = np.linalg.norm(positions[-1])

    def bdy (nu):
        return tools.I_CMB (nu)

    def z_max(r, theta):
        if (r>=r_out):#image point outside domain
            return 0.0
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
        fig = plt.figure(dpi=150)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(fs, im  [p,:], marker='.')
        ax.plot(fs, im_a[p,:])
        ax.text(.05, .95,"(x,y)={x},{y}".format(x=imx[p], y=imy[p]), transform=ax.transAxes)

    npix = np.size(imx)
    #during automated testing, the widgets only consume time to create
    if use_widgets:
        widgets.interact(plot, p=(0,npix-1))#This plot the benchmark results for each pixel seperately.

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

    #note: this benchmark is for making sure that the new imager does not produce nans, when applied to 3D models
    FEAUTRIER_AS_EXPECTED=(np.max(error)<1.5e-3)

    return (FEAUTRIER_AS_EXPECTED)


def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
