import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'
sys.path.append(f'{curdir}/../../../')


import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import ipywidgets        as widgets
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte


dimension = 1
npoints   = 100
nrays     = 100
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

    modelName = f'constant_velocity_gradient_1D'
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

    return


def run_model (nosave=False):

    modelName = f'constant_velocity_gradient_1D'
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

    npoints = model.parameters.npoints()
    nrays   = model.parameters.nrays  ()
    hnrays  = model.parameters.hnrays ()
    nfreqs  = model.parameters.nfreqs ()

    timer3 = tools.Timer('shortchar 0  ')
    timer3.start()
    model.compute_radiation_field_shortchar_order_0 ()
    timer3.stop()
    u_0s = np.array(model.radiation.u)

    timer4 = tools.Timer('feautrier 2  ')
    timer4.start()
    model.compute_radiation_field_feautrier_order_2 ()
    timer4.stop()
    u_2f = np.array(model.radiation.u)

    rs = np.array(model.geometry.points.position)[:,0]
    nu = np.array(model.radiation.frequencies.nu)[0]

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
        return chi*L / (fct*dnu) * 0.5 * (sp.special.erf(arg) + sp.special.erf(fct*l/L-arg))

    def I_ (nu, r, theta):
        return src + (bdy(nu)-src)*np.exp(-tau(nu, r, theta))

    def u_ (nu, r, theta):
        return 0.5 * (I_(nu, r, theta) + I_(nu, r, np.pi-theta))

    rx, ry, rz = np.array(model.geometry.rays.direction).T
    angles     = np.arctan2(ry,rx)
    angles     = angles[:hnrays]

    us = np.array([[[u_(f,r,a) for f in nu] for r in rs] for a in angles])
    fs = (nu-frq)/frq*magritte.CC

    def plot (r, p):
        plt.figure(dpi=150)
        plt.plot(fs, us  [r,p,:], marker='.')
        plt.plot(fs, u_2f[r,p,:])
    widgets.interact(plot, r=(0,hnrays-1,1), p=(0,npoints-1,1))

    error_u_0s = np.abs(tools.relative_error(us, u_0s))
    error_u_2f = np.abs(tools.relative_error(us, u_2f))

    log_err_min = np.log10(np.min([error_u_0s, error_u_2f]))
    log_err_max = np.log10(np.max([error_u_0s, error_u_2f]))

    bins = np.logspace(log_err_min, log_err_max, 100)

    result  = f'--- Benchmark name ------------------------------\n'
    result += f'{modelName                                      }\n'
    result += f'--- Parameters ----------------------------------\n'
    result += f'dimension = {model.parameters.dimension()       }\n'
    result += f'npoints   = {model.parameters.npoints  ()       }\n'
    result += f'nrays     = {model.parameters.nrays    ()       }\n'
    result += f'nquads    = {model.parameters.nquads   ()       }\n'
    result += f'--- Accuracy ------------------------------------\n'
    result += f'mean error in shortchar 0 = {np.mean(error_u_0s)}\n'
    result += f'mean error in feautrier 2 = {np.mean(error_u_2f)}\n'
    result += f'--- Timers --------------------------------------\n'
    result += f'{timer1.print()                                 }\n'
    result += f'{timer2.print()                                 }\n'
    result += f'{timer3.print()                                 }\n'
    result += f'{timer4.print()                                 }\n'
    result += f'-------------------------------------------------\n'

    print(result)

    if not nosave:
        with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
            log.write(result)

        plt.figure()
        plt.title(modelName)
        plt.hist(error_u_0s.ravel(), bins=bins, histtype='step', label='0s')
        plt.hist(error_u_2f.ravel(), bins=bins, histtype='step', label='2f')
        plt.xscale('log')
        plt.legend()
        plt.savefig(f'{resdir}{modelName}_hist-{timestamp}.png', dpi=150)

        plt.figure()
        plt.hist(error_u_0s.ravel(), bins=bins, histtype='step', label='0s', cumulative=True)
        plt.hist(error_u_2f.ravel(), bins=bins, histtype='step', label='2f', cumulative=True)
        plt.xscale('log')
        plt.legend()
        plt.savefig(f'{resdir}{modelName}_cumu-{timestamp}.png', dpi=150)

    return


def run_test (nosave=False):

    create_model     ()
    return run_model (nosave)


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
