import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'
sys.path.append(f'{curdir}/../../../')


import numpy             as np
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte


dimension = 1
npoints   = 5 #1000
nrays     = 2
nspecs    = 5
nlspecs   = 1
nquads    = 1

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+10                 # [m^-3]
temp = 4.5E+00                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+12                 # [m]
dv   = 0.0E+00 / magritte.CC   # [fraction of speed of light]


def create_model ():
    """
    Create a model file for the all_constant benchmark in 1D.
    """

    modelName = f'all_constant_1D'
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

    model.chemistry.species.abundance = [[     0.0,    nTT,  nH2,  0.0,      1.0] for _ in range(npoints)]
    model.chemistry.species.symbol    =  ['dummy0', 'test', 'H2', 'e-', 'dummy1']

    model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
    model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

    model = setup.set_Delaunay_neighbor_lists (model)
    model = setup.set_Delaunay_boundary       (model)
    model = setup.set_boundary_condition_CMB  (model)
    model = setup.set_uniform_rays            (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
    model = setup.set_quadrature              (model)

    model.write()

    return


def run_model ():

    modelName = f'all_constant_1D'
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
    model.compute_radiation_field_0th_short_characteristics ()
    timer3.stop()
    u_0s = np.array(model.radiation.u)

    timer4 = tools.Timer('feautrier 2  ')
    timer4.start()
    model.compute_radiation_field_2nd_order_Feautrier ()
    timer4.stop()
    u_2f = np.array(model.radiation.u)

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

    def I_0 (x):
        return src + (bdy-src)*np.exp(-chi*x)

    def I_1 (x):
        return src + (bdy-src)*np.exp(-chi*(x[-1]-x))

    def u_ (x):
        return 0.5 * (I_0(x) + I_1(x))

    error_u_0s = tools.relative_error (u_(x), u_0s[0,:,0])
    error_u_2f = tools.relative_error (u_(x), u_2f[0,:,0])

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

    with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
        log.write(result)

    return


if __name__ == '__main__':

    create_model ()
    run_model    ()