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
import ipywidgets        as widgets


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

    modelName = f'density_distribution_VZ{a_or_b}_1D'
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

    model.chemistry.species.abundance = [[nTT(r), nH2(r),  0.0] for r in rs]
    model.chemistry.species.symbol    = ['test', 'H2', 'e-']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    setup.set_Delaunay_neighbor_lists (model)
    setup.set_Delaunay_boundary       (model)
    setup.set_boundary_condition_CMB  (model)
    setup.set_rays_spherical_symmetry (model)
    setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    setup.set_quadrature              (model)

    model.write()

    return #magritte.Model (modelFile)


def run_model (a_or_b, nosave=False, use_widgets=True):

    modelName = f'density_distribution_VZ{a_or_b}_1D'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    X_mol = get_X_mol[a_or_b]

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
            return pref * (1.0/r_in - 1.0/r)
        if (theta == np.pi):
            return pref * (1.0/r - 1.0/r_out)
        #if intersecting with the inner core
        if (theta < np.arcsin(r_in/r)):
            closest = r * np.sin(theta)
            in_vertical = np.sqrt(r_in**2 - closest**2)
            out_vertical = np.sqrt(r**2 - closest**2)
            return pref * (np.arctan(out_vertical/closest) - np.arctan(in_vertical/closest)) / closest
        else:
            #Other derivation: results in exactly the same error, so should be equivalent
            # closest = r * np.sin(theta)
            # out_vertical = np.sqrt(r_out**2 - closest**2)
            # start_vertical = np.sqrt(r**2 - closest**2)
            # if np.abs(theta)<np.pi/2:
            #     return pref * (np.arctan(out_vertical/closest) + np.arctan(start_vertical/closest)) / closest
            # else:
            #     return pref * (np.arctan(out_vertical/closest) - np.arctan(start_vertical/closest)) / closest

            #Old derivation
            factor = np.arccos(r*np.sin(theta)/r_out) + 0.5*np.pi - theta
            return pref / (r*np.sin(theta)) * factor

    def I_ (nu, r, theta):
        return src + (bdy(nu)-src)*np.exp(-tau(nu, r, theta))

    def u_ (nu, r, theta):
        return 0.5 * (I_(nu, r, theta) + I_(nu, r, np.pi+theta))

    rx, ry, rz = np.array(model.geometry.rays.direction).T
    angles     = np.arctan2(ry,rx)
    angles     = angles[:hnrays]

    us = np.array([[[u_(f,r,a) for f in nu] for r in rs] for a in angles])

    fs = (nu-frq)/frq*magritte.CC
    #hmm, do we need this widget if we are not plotting it?
    def plot (r, p):
        plt.figure(dpi=150)
        plt.plot(fs, us  [r,p,:], marker='.')
        plt.plot(fs, u_2f[r,p,:])
        # plt.plot(fs, u_0s[r,p,:])

    #during automated testing, the widgets only consume time to create
    if use_widgets:
        widgets.interact(plot, r=(0,hnrays-1,1), p=(0,npoints-1,1))

    error_u_0s = np.abs(tools.relative_error(us, u_0s)[:-1, :-1, :])
    error_u_2f = np.abs(tools.relative_error(us, u_2f)[:-1, :-1, :])

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

    #setting 'b' not yet used for testing
    if a_or_b == 'a':
        #error bounds are chosen somewhat arbitrarily, based on previously obtained results; this should prevent serious regressions.
        FEAUTRIER_AS_EXPECTED=(np.mean(error_u_2f)<4.4e-4)
        FIRSTORDER_AS_EXPECTED=(np.mean(error_u_0s)<4.1e-4)

        if not FIRSTORDER_AS_EXPECTED:
            print("First order solver mean error too large: ", np.mean(error_u_0s))
        if not FEAUTRIER_AS_EXPECTED:
            print("Feautrier solver mean error too large: ", np.mean(error_u_2f))

        return (FEAUTRIER_AS_EXPECTED&FIRSTORDER_AS_EXPECTED)

    return


def run_test (nosave=False):

    create_model ('a')
    run_model    ('a', nosave)

    create_model ('b')
    run_model    ('b', nosave)


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
