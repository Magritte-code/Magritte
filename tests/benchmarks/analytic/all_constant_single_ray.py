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
npoints   = 40
nrays     = 2
nspecs    = 5
nlspecs   = 1
nquads    = 21

nH2  = 1.0E+12                 # [m^-3]
nTT  = 1.0E+03                 # [m^-3]
temp = 4.5E+00                 # [K]
turb = 0.0E+00                 # [m/s]
dx   = 1.0E+12                 # [m]
dv   = 0.0E+00 / magritte.CC   # [fraction of speed of light]
# dv   = 2.5E+01 / magritte.CC   # [fraction of speed of light]


def create_model ():
    """
    Create a model file for the all_constant benchmark, single ray.
    """

    modelName = f'all_constant_single_ray'
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

    return #magritte.Model (modelFile)


def run_model (nosave=False):

    modelName = f'all_constant_single_ray'
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
    model.compute_Jeff();
    u_0s = np.array(model.radiation.u)
    sum_J_0s=np.array(model.lines.lineProducingSpecies[0].Jlin)

    timer4 = tools.Timer('feautrier 2  ')
    timer4.start()
    model.compute_radiation_field_feautrier_order_2 ()
    timer4.stop()
    model.compute_Jeff();
    u_2f = np.array(model.radiation.u)
    J_2f = np.array(model.radiation.J)
    sum_J_2f=np.array(model.lines.lineProducingSpecies[0].Jlin)
    quads= np.array(model.lines.lineProducingSpecies[0].quadrature.weights)
    # sum_J_2f= np.sum(J_2f*quads, axis=1)
    #approx; because feautrier solver does not compute the intensities...
    # I_2f = np.array(model.radiation.I)*u_2f[0,:,0]/u_0s[0,:,0]

    timer5 = tools.Timer('collocation  ')
    timer5.start()
    # model.COMPUTING_SVD_AND_CONDITION_NUMBER=True
    model.compute_radiation_field_collocation ()
    timer5.stop()
    colf = np.array(model.lines.lineProducingSpecies[0].Jlin)
    # u_col= np.array(model.radiation.u)

    x  = np.array(model.geometry.points.position)[:,0]
    nu = np.array(model.radiation.frequencies.nu)
    print(colf[:,0])
    print(x)

    ld = model.lines.lineProducingSpecies[0].linedata

    k = 0

    frq = ld.frequency[k]
    pop = tools.LTEpop         (ld, temp) * nTT
    phi = tools.profile        (ld, k, temp, (turb/magritte.CC)**2, frq)
    eta = tools.lineEmissivity (ld, pop)[k] * phi
    chi = tools.lineOpacity    (ld, pop)[k] * phi
    src = tools.lineSource     (ld, pop)[k]
    bdy = tools.I_CMB          (frq)

    velocities=np.array(model.geometry.points.velocity);
    doppler_shifts=1+velocities[:,0];
    doppler_shifts_dist=(doppler_shifts*frq-frq)/tools.dnu(ld, k, temp, (turb/magritte.CC)**2)
    dopp_mult_factors=tools.profile(ld, k, temp, (turb/magritte.CC)**2, frq*doppler_shifts)/phi
    # dopp_exp_factor=np.exp(-(1-dopp_mult_factors)**2)
    # print("dopp exp factor", dopp_exp_factor)

    # print(dopp_mult_factors)
    # print(u_2f)
    # rev_dopp_mult_factors=dopp_exp_factor[::-1]
    #u is the sum of two intensities on the opposite directions, so we actually need to doppler shift them individually
    #also do not forget the background intensity; TODO;: figure out better formula
    # u_2f_dopp=(I_2f[0,:,0]*dopp_mult_factors+bdy*(1-dopp_mult_factors)+I_2f[1,:,0]*rev_dopp_mult_factors+bdy*(1-rev_dopp_mult_factors))/2
    # u_2f_dopp=u_2f[0,:,0]*dopp_mult_factors
    # u_2f_dopp=(I_2f[0,:,0]*dopp_mult_factors+I_2f[1,:,0]*rev_dopp_mult_factors)/2
    #if i accidentally multiplied the wrong thing
    # other_u_2f_dopp=(I_2f[0,:,0]*rev_dopp_mult_factors+bdy*(1-rev_dopp_mult_factors)+I_2f[1,:,0]*dopp_mult_factors+bdy*(1-dopp_mult_factors))/2
    # other_u_2f_dopp=(I_2f[0,:,0]*rev_dopp_mult_factors+I_2f[1,:,0]*dopp_mult_factors)/2

    # print("antipods", np.array(model.geometry.rays.antipod))
    # print(I_2f[:,:,0])
    # print(dopp_mult_factors)
    # print(u_2f[0,:,0]/u_0s[0,:,0])
    # print(bdy)
    # print(u_2f_dopp)
    # print(other_u_2f_dopp)
    print("quads: ",quads)
    print("J_2f: ",J_2f)
    print("integrated: ", sum_J_2f)
    print("quadroots", np.array(model.lines.lineProducingSpecies[0].quadrature.roots))
    print("doppler_shifts_dist", doppler_shifts_dist)

    def I_0 (x):
        return src + (bdy-src)*np.exp(-chi*x)

    def I_1 (x):
        return src + (bdy-src)*np.exp(-chi*(x[-1]-x))

    def u_ (x):
        return 0.5 * (I_0(x) + I_1(x))

    error_u_0s = tools.relative_error (u_(x), u_0s[0,:,0])
    error_u_2f = tools.relative_error (u_(x), u_2f[0,:,0])
    error_Jcol = tools.relative_error (u_(x), colf[:,0  ])

    result  = f'--- Benchmark name ----------------------------\n'
    result += f'{modelName                                    }\n'
    result += f'--- Parameters --------------------------------\n'
    result += f'dimension = {model.parameters.dimension()     }\n'
    result += f'npoints   = {model.parameters.npoints  ()     }\n'
    result += f'nrays     = {model.parameters.nrays    ()     }\n'
    result += f'nquads    = {model.parameters.nquads   ()     }\n'
    result += f'--- Accuracy ----------------------------------\n'
    result += f'max error in shortchar 0 = {np.max(abs(error_u_0s))}\n'
    result += f'max error in feautrier 2 = {np.max(abs(error_u_2f))}\n'
    result += f'max error in collocation = {np.max(abs(error_Jcol))}\n'
    result += f'--- Timers ------------------------------------\n'
    result += f'{timer1.print()                               }\n'
    result += f'{timer2.print()                               }\n'
    result += f'{timer3.print()                               }\n'
    result += f'{timer4.print()                               }\n'
    result += f'{timer5.print()                               }\n'
    result += f'-----------------------------------------------\n'

    print(result)

    mid=int(model.parameters.nquads()/2)
    print(mid)

    if not nosave:
        with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
            log.write(result)

        fig = plt.figure(dpi=150)
        # plt.title(modelName)
        plt.title("Small doppler shift")
        # plt.scatter(x, u_0s[0,:,mid], s=0.5, label='0s', zorder=1)
        # plt.scatter(x, u_2f[0,:,mid], s=0.5, label='2f', zorder=1)
        plt.scatter(x, sum_J_2f[:], s=0.5, label='J_2f', zorder=1)
        plt.scatter(x, colf[:,0  ], s=0.5, label='J_col', zorder=1)
        plt.scatter(x, sum_J_0s[:], s=0.5, label='J_0s', zorder=1)
        # plt.scatter(x, u_2f[0,:,0+2], s=0.5, label='2f+2', zorder=1)
        # plt.scatter(x, u_2f[0,:,0+3], s=0.5, label='2f+3', zorder=1)
        # plt.scatter(x, u_2f[0,:,0+4], s=0.5, label='2f+4', zorder=1)

        #these things are currently terrible approximations
        # plt.scatter(x, u_2f_dopp, s=0.5, label='u_2f_dopp', zorder=1)
        # plt.scatter(x, other_u_2f_dopp, s=0.5, label='u_2f_dopp', zorder=1)
        # plt.scatter(x, u_col[0,:,0], s=0.5, label='ucol', zorder=1)

        # plt.plot(x, u_(x), c='lightgray', zorder=0)

        # plt.plot(x, u_(x)/np.sqrt(np.sqrt(2)), c='lightgray', zorder=0)
        plt.legend()
        # plt.xscale('log')
        plt.xlabel('r [m]')
        # plt.ylabel('Mean intensity [W/m$^{2}$]/Mean effective intensity [W]')

        plt.ylabel('$J$ [J/s]')

        plt.show()
        # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

    return


def run_test (nosave=False):

    create_model ()
    run_model    (nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
