import os
import sys
#multiline
curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'

import numpy             as np
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte

#Note: the discretization of in the spatial dimension is important, as the frequency quadrature should be sufficiently overlapping from point to point
# TODO: turn this file into an automated test, checking the difference between comoving and the regular solvers.
dimension = 1
npoints   = 3*50
nrays     = 2
nspecs    = 5
nlspecs   = 1
nquads    = 25

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+03                 # [m^-3]
temp = 4.5E+00                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+12                 # [m]
# dv   = 0.0E+00 / magritte.CC   # [fraction of speed of light]
dv   = 1.0E+02 / magritte.CC   # [fraction of speed of light]
vmax = 1.5*1.0E+03 / magritte.CC
# vmax = 0.0E+03 / magritte.CC
period = 3*20
# max_line_width = 5.0


def create_model ():
    """
    Create a model file for the all_constant benchmark, single ray.
    """

    modelName = f'all_constant_single_ray_testing_comoving_multiline'
    modelFile = f'{moddir}{modelName}.hdf5'
    # lamdaFile = f'{datdir}co.txt'
    lamdaFile = f'{datdir}test2.txt'


    model = magritte.Model ()
    model.parameters.set_spherical_symmetry(False)
    model.parameters.set_model_name        (modelFile)
    model.parameters.set_dimension         (dimension)
    model.parameters.set_npoints           (npoints)
    model.parameters.set_nrays             (nrays)
    model.parameters.set_nspecs            (nspecs)
    model.parameters.set_nlspecs           (nlspecs)
    model.parameters.set_nquads            (nquads)
    # model.parameters.set_nlines            (2)
    # model.parameters.set_hnrays            (1)
    # model.parameters.set_nfreqs            (2*nquads)


    model.geometry.points.position.set([[i*dx, 0, 0] for i in range(npoints)])
    # model.geometry.points.velocity.set([[i*dv, 0, 0] for i in range(npoints)])
    vsin = [[vmax*np.sin(2.0*np.pi*i/period), 0, 0] for i in range (npoints)]
    model.geometry.points.velocity.set(vsin)
    # def wedgevelocity (i, period):
    #     modi = i%period
    #     if modi<period/2:
    #         return modi/(period)
    #     else:
    #         return 1.0 - modi/(period)


    # vwedge = [[vmax*wedgevelocity(i, period),0,0] for i in range (npoints)]
    # print(vwedge)
    # model.geometry.points.velocity.set(vwedge)



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
    # model = setup.set_quadrature_equidistant (model, max_line_width)

    model.write()

    return #magritte.Model (modelFile)


def run_model (nosave=False):


    modelName = f'all_constant_single_ray_testing_comoving_multiline'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    print("here")

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()
    magritte.pcmt_set_n_threads_avail(1)#for fairness, single threaded is a better comparison

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()
    inv_line_width = np.array(model.lines.inverse_width)
    line_freq = np.array(model.lines.line)
    max_line_width = np.max(np.array(model.lines.lineProducingSpecies[0].quadrature.roots))
    print("line width: ", 1.0/inv_line_width)
    print("inv width: ", inv_line_width[0])
    print("max velocity diff: ", 2.0 * vmax)
    print("relative width: ", line_freq[0]*inv_line_width[0])
    print("quad roots: ", np.array(model.lines.lineProducingSpecies[0].quadrature.roots))
    print("quad weights: ", np.array(model.lines.lineProducingSpecies[0].quadrature.weights))
    print("Icmb: ", tools.I_CMB(line_freq[0]), tools.I_CMB(line_freq[1]))
    print("dfreqs for line 0:", line_freq[0]* inv_line_width[0, 0], "dfreqs for line 1:", line_freq[1]*inv_line_width[0,1])
    print("in units of div c", 1/(line_freq[0]* inv_line_width[0, 0]), 1/(line_freq[1]*inv_line_width[0,1]))
    print("times max line width", 1/(line_freq[0]* inv_line_width[0, 0])*max_line_width, 1/(line_freq[1]*inv_line_width[0,1])*max_line_width)

    timer3 = tools.Timer('shortchar 0  ')
    timer3.start()
    model.compute_radiation_field_shortchar_order_0 ()
    model.compute_Jeff()
    timer3.stop()
    u_0s = np.array(model.radiation.u)
    J_0s = np.array(model.lines.lineProducingSpecies[0].Jlin)

    timer5 = tools.Timer('comoving  ')
    timer5.start()
    model.compute_radiation_field_comoving ()
    model.compute_Jeff_sparse()
    timer5.stop()
    J_co = np.array(model.lines.lineProducingSpecies[0].Jlin)
    print(timer5.print())

    timer6 = tools.Timer('comoving approx ')
    timer6.start()
    # model.compute_radiation_field_comoving_local_approx ()
    model.compute_Jeff_sparse()
    timer6.stop()
    J_coap = np.array(model.lines.lineProducingSpecies[0].Jlin)
    print(timer6.print())
    print(J_coap)

    timer4 = tools.Timer('feautrier 2  ')
    timer4.start()
    model.compute_radiation_field_feautrier_order_2 ()
    model.compute_Jeff();
    timer4.stop()
    u_2f = np.array(model.radiation.u)
    J_2f = np.array(model.lines.lineProducingSpecies[0].Jlin)
#TODO: figure out a way to compare this stuff

    # print(J_0s)
    # print(J_2f)
    # print(J_co)

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
    result += f'{timer5.print()                               }\n'
    result += f'{timer6.print()                               }\n'
    result += f'-----------------------------------------------\n'

    print(result)

    if not nosave:
        with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
            log.write(result)

        fig = plt.figure(dpi=150)
        plt.title(modelName+' line 0')
        # plt.scatter(x, u_0s[0,:,0], s=0.5, label='0s', zorder=1)
        # plt.scatter(x, u_2f[0,:,0], s=0.5, label='2f', zorder=1)
        plt.scatter(x, J_0s[:,0], s=0.5, label='0s', zorder=1)
        plt.scatter(x, J_2f[:,0], s=0.5, label='2f', zorder=1)
        plt.scatter(x, J_co[:,0], s=0.5, label='co', zorder=1)
        # plt.scatter(x, J_coap[:,0], s=0.5, label='co-ap', zorder=1)
        # plt.plot(x, u_(x), c='lightgray', zorder=0)
        plt.legend()
        # plt.xscale('log')
        plt.xlabel('r [m]')
        plt.ylabel('Mean intensity [W/m$^{2}$]')
        fig = plt.figure(dpi=150)
        plt.title(modelName+' line 1')
        # plt.scatter(x, u_0s[0,:,0], s=0.5, label='0s', zorder=1)
        # plt.scatter(x, u_2f[0,:,0], s=0.5, label='2f', zorder=1)
        plt.scatter(x, J_0s[:,1], s=0.5, label='0s', zorder=1)
        plt.scatter(x, J_2f[:,1], s=0.5, label='2f', zorder=1)
        plt.scatter(x, J_co[:,1], s=0.5, label='co', zorder=1)
        # plt.scatter(x, J_coap[:,0], s=0.5, label='co-ap', zorder=1)
        # plt.plot(x, u_(x), c='lightgray', zorder=0)
        plt.legend()
        # plt.xscale('log')
        plt.xlabel('r [m]')
        plt.ylabel('Mean intensity [W/m$^{2}$]')
        #TODO PLOT J instead
        plt.figure()
        plt.plot(np.array(model.geometry.points.velocity)[:,0])
        plt.plot(np.array(model.geometry.points.velocity)[:,0]+max_line_width / (line_freq[0] * inv_line_width[0, 0]))
        plt.plot(np.array(model.geometry.points.velocity)[:,0]-max_line_width / (line_freq[0] * inv_line_width[0, 0]))
        # plt.plot(np.array(model.geometry.points.velocity)[:,0]+max_line_width * inv_line_width[0, 0]/np.sqrt(np.pi))
        # plt.plot(np.array(model.geometry.points.velocity)[:,0]-max_line_width * inv_line_width[0, 0]/np.sqrt(np.pi))
        # plt.axhline(y=inv_line_width[0, 0], color='r', linestyle='-')
        # plt.axhline(y=-inv_line_width[0, 0], color='r', linestyle='-')


        plt.show()
        # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

    #returning whether output is as expected (not too far from the input)
    # max_diff=max(u_(x))-min(u_(x))
    # print(max_diff)
    #error should at max be proportional to max diff? only maybe for testing non analytic models

    #error bounds are chosen somewhat arbitrarily, based on previously obtained results; this should prevent serious regressions.
    #feautrier should be quite accurate, so we take 1e-5 as maximal max error
    FEAUTRIER_AS_EXPECTED=(np.max(error_u_2f)<1e-5)
    #not that well implmented, so less accurate
    FIRSTORDER_AS_EXPECTED=(np.max(error_u_0s)<0.1)

    if not FIRSTORDER_AS_EXPECTED:
        print("First order solver max error too large: ", np.max(error_u_0s))
    if not FEAUTRIER_AS_EXPECTED:
        print("Feautrier solver max error too large: ", np.max(error_u_2f))


    return (FEAUTRIER_AS_EXPECTED&FIRSTORDER_AS_EXPECTED)


def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
