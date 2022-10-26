import os
import sys

curdir = os.path.dirname(os.path.realpath(__file__))
datdir = f'{curdir}/../../data/'
moddir = f'{curdir}/../../models/'
resdir = f'{curdir}/../../results/'

import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte
import random

from scipy.interpolate import interp1d


dimension = 1
npoints   = 20#250
nrays     = 10#50
nspecs    = 5
nlspecs   = 1
nquads    = 11#11

r_in   = 1.0E13   # [m]
r_out  = 7.8E16 * 100   # [m]
#for testing arbitrary density distributions
r_mid1 = 2.0E15
r_mid2 = 8.0E16
# nH2_in = 2.0E28 * 100  # [m^-3]
nH2_in = 8.0E13   # [m^-3]
temp   =  20.00   # [K]
# temp   =  200.00   # [K]
turb   = 150.00   # [.]
dv = 0.0
dv   = 1.5E+02 / magritte.CC   # [fraction of speed of light]

get_X_mol = {
    'a' : 1.0E-8,
    'b' : 1.0E-6
}

rs = np.logspace (np.log10(r_in), np.log10(r_out), npoints, endpoint=True)


def create_models (a_or_b):
    """
    Create a model file for the van Zadelhoff 1 benchmark in 1D.
    """

    for i in range(1):
        modelName = f'vanZadelhoff_1{a_or_b}_1D_testmol{i}temp'
        modelFile = f'{moddir}{modelName}.hdf5'
        # lamdaFile = f'{datdir}pH20_20.txt'
        # lamdaFile = f'{datdir}o-h2o_nu0_nu2_nu3.dat'
        lamdaFile = f'{datdir}o-h2o_nu0_nu2_nu3_pruned40.dat'
        # lamdaFile = f'{datdir}co.txt'
        # lamdaFile = f'{datdir}test.txt'
        # lamdaFile = f'{datdir}testmaser.txt'

        print("create model ", i)

        X_mol = get_X_mol[a_or_b]


        def Temp(r, i):
            return 2000 * np.sqrt(r_in/r)
            if (r<r_mid1):
                temparr = [100, 500, 2000]
                return temparr[i]
            if (r<r_mid2):
                temparr = [500, 2000, 100]
                return temparr[i]
            else:
                temparr = [2000, 100, 500]
                return temparr[i]

        def nH2 (r):
            return nH2_in * np.power(r_in/r, 2.0)
            # return nH2_in / 1e13

            if (r<r_mid2 and r>r_mid1):
                return nH2_in* np.power(r_in/r, 2.0) * 1e5
            else:
                return nH2_in* np.power(r_in/r, 2.0)

        def nTT (r):
            return X_mol  * nH2(r)
            # return X_mol  * nH2_in * np.power(r_in/r, 2.0)
            # return X_mol  * nH2_in / 1.0e10

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

        velocities = [[i*dv, 0, 0] for i in range(len(rs))]
        #testing random velocity profile
        # random.shuffle(velocities)
        # print(velocities)

        densities = [[     0.0, nTT(r), nH2(r),  0.0,      1.0] for r in rs]
        # random.shuffle(densities)

        model.geometry.points.velocity.set(velocities)

        # model.chemistry.species.abundance = [[     0.0, nTT(r), nH2(r),  0.0,      1.0] for r in rs]
        model.chemistry.species.abundance = densities
        model.chemistry.species.symbol    =  ['dummy0', 'oH2O',   'H2', 'e-', 'dummy1']
        # model.chemistry.species.symbol    =  ['dummy0', 'CO',   'H2', 'e-', 'dummy1']
        # model.chemistry.species.symbol    =  ['dummy0', 'test',   'H2', 'e-', 'dummy1']

        # model.parameters.one_line_approximation = True

        # model.thermodynamics.temperature.gas  .set( temp                 * np.ones(npoints))
        model.thermodynamics.temperature.gas  .set([Temp(r, i) for r in rs])
        model.thermodynamics.turbulence.vturb2.set((turb/magritte.CC)**2 * np.ones(npoints))

        model = setup.set_Delaunay_neighbor_lists (model)
        model = setup.set_Delaunay_boundary       (model)
        model = setup.set_boundary_condition_CMB  (model)
        # model = setup.set_boundary_condition_1D (model, Temp(r_in, i), Temp(r_out, i))
        model = setup.set_boundary_condition_1D (model, 2000, 50)
        model = setup.set_rays_spherical_symmetry (model)
        model = setup.set_linedata_from_LAMDA_file(model, lamdaFile)
        model = setup.set_quadrature              (model)

        model.write()



    return #magritte.Model (modelFile)


def run_models (a_or_b, nosave=False):

    temps = np.array([])
    cooling_rates = np.array([])
    abundances = np.array([])

    # magritte.pcmt_set_n_threads_avail(1)
    for i in range(1):
        modelName = f'vanZadelhoff_1{a_or_b}_1D_testmol{i}temp'
        modelFile = f'{moddir}{modelName}.hdf5'
        timestamp = tools.timestamp()

        timer1 = tools.Timer('reading model')
        timer1.start()
        model = magritte.Model (modelFile)
        timer1.stop()

        model.parameters.convergence_fraction=0.9995
        model.parameters.pop_prec=1e-7
        model.parameters.min_negative_dtau = -5.0;

        timer2 = tools.Timer('setting model')
        timer2.start()
        model.compute_spectral_discretisation ()
        model.compute_inverse_line_widths     ()
        model.compute_LTE_level_populations   ()
        timer2.stop()


        npoints = model.parameters.npoints()
        nlev    = model.lines.lineProducingSpecies[0].linedata.nlev
        LTEpops = np.array(model.lines.lineProducingSpecies[0].population).reshape((npoints, nlev))
        # pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((npoints, nlev))
        print(LTEpops)
        inv_widths = np.array(model.lines.inverse_width)
        print("inverse widths: ", inv_widths[1,:])
        lines = np.array(model.lines.line)
        print("line freqs: ", lines)
        print("minimal diff: ", np.diff(np.sort(lines)).min())

        timer3 = tools.Timer('running model')
        # model.parameters.one_line_approximation = True
        timer3.start()
        # model.compute_level_populations (True, 100)
        model.compute_level_populations_shortchar (True, 2000)
        # model.parameters.one_line_approximation = False
        # model.compute_level_populations_shortchar (True, 100)

        timer3.stop()

        Jshort = np.array(model.radiation.J)
        ushort = np.array(model.radiation.u)
        print("shortchar J: ", Jshort)
        print("shortchar negative J's: ", Jshort[Jshort<0])
        print("shortchar u: ", ushort)
        print("shortchar negative u's: ", ushort[ushort<0])



        model.compute_image_optical_depth(round(nrays/2-1));
        model.compute_image(round(nrays/2-1));

        optical_depths=np.array(model.images[0].I)
        us=np.array(model.images[1].I)

        print("opt depths:", optical_depths)
        print("negative opt depths:", optical_depths[optical_depths<0])
        plt.figure()
        plt.hist(np.log10(-optical_depths[optical_depths<0]))
        # plt.xscale('log')
        plt.xlabel('log τ [1]')
        print("u: ", us)
        print("u with negative optical depths: ", us[optical_depths<0])
        # print("maximum negative opt depth:", min(optical_depths[optical_depths<0]))
        emiss = np.array(model.lines.emissivity)
        print("emissivities: ", emiss)
        print("negative emissivities: ", emiss[emiss<0.0])

        J = np.array(model.radiation.J)
        # print("J:", J)
        print("negative J's image: ", J[J<0])
        model.compute_radiation_field_shortchar_order_0()
        I_shortchar = np.array(model.radiation.I)
        freqs = np.array(model.radiation.frequencies.nu)
        print("I shortchar: ", I_shortchar)
        print("negative I shortchar: ", I_shortchar[I_shortchar<0.0])
        plt.figure()
        # plt.scatter(freqs[0,:], I_shortchar[round(nrays/2-1),0,:])
        plt.scatter(freqs[0,:], I_shortchar[0,0,:])
        # plt.plot(freqs[0,:], Jshort[0,:])
        plt.xscale('log')
        plt.xlabel('ν [/s]')
        plt.yscale('log')
        plt.ylabel('I [J/(m^2)] per ν')

        plt.figure()
        # plt.plot(us[0,:])
        plt.title("u")
        plt.scatter(freqs[0,:], ushort[1,0,:])
        plt.xscale('log')
        plt.xlabel('ν [/s]')
        plt.yscale('log')

        plt.figure()
        plt.title("optical depths")
        plt.scatter(freqs[0,:],optical_depths[0,:])
        plt.yscale("symlog")
        plt.xscale('log')
        plt.xlabel('ν [/s]')

        # print("n off diag: ", model.parameters.n_off_diag)
        # print("lambda size: ", model.lines.lineProducingSpecies[0].Lambda.Ls)

        print("computing collisional")
        # now computing cooling as a test
        model.compute_cooling_collisional()
        coll_cooling_rate=np.array(model.cooling.cooling_rate);
        print("computing radiative")
        model.compute_cooling_radiative_shortchar()
        rad_cooling_rate_shortchar=np.array(model.cooling.cooling_rate);
        model.compute_cooling_radiative()
        rad_cooling_rate=np.array(model.cooling.cooling_rate);


        print("coll_cooling: ",coll_cooling_rate)
        print("rad_cooling: ", rad_cooling_rate)
        print("rad_cooling_sc: ", rad_cooling_rate_shortchar)


        pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((npoints, nlev))
        abun = np.array(model.chemistry.species.abundance)[:,1]
        rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)


        cooling_rates=np.append(cooling_rates, coll_cooling_rate)
        temps=np.append(temps, np.array(model.thermodynamics.temperature.gas))
        abundances=np.append(abundances, abun)

        print("levelpops last point: ", pops[-1,:])
        print("LTE levelpops last p: ", LTEpops[-1,:])
        print("negative levelpops: ", pops[pops<0.0])

        # plt.figure()
        # plt.plot((pops[0,:]-LTEpops[0,:])/(LTEpops[0,:]))
        # plt.title("relative difference levelpops vs LTE")

        (i,ra,rb,nh,tk,nm,vr,db,td,lp0,lp1) = np.loadtxt (f'{curdir}/Ratran_results/vanZadelhoff_1{a_or_b}.out', skiprows=14, unpack=True)

        interp_0 = interp1d(0.5*(ra+rb), lp0, fill_value='extrapolate')
        interp_1 = interp1d(0.5*(ra+rb), lp1, fill_value='extrapolate')

        error_0 = tools.relative_error(pops[:,0]/abun, interp_0(rs))
        error_1 = tools.relative_error(pops[:,1]/abun, interp_1(rs))

        result  = f'--- Benchmark name -----------------------\n'
        result += f'{modelName                               }\n'
        result += f'--- Parameters ---------------------------\n'
        result += f'dimension = {model.parameters.dimension()}\n'
        result += f'npoints   = {model.parameters.npoints  ()}\n'
        result += f'nrays     = {model.parameters.nrays    ()}\n'
        result += f'nquads    = {model.parameters.nquads   ()}\n'
        result += f'--- Accuracy -----------------------------\n'
        result += f'max error in (0) = {np.max(error_0[1:])  }\n'
        result += f'max error in (1) = {np.max(error_1[1:])  }\n'
        result += f'--- Timers -------------------------------\n'
        result += f'{timer1.print()                          }\n'
        result += f'{timer2.print()                          }\n'
        result += f'{timer3.print()                          }\n'
        result += f'------------------------------------------\n'

        print(result)

        if not nosave:

            # with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
                # log.write(result)

            # plt.figure(dpi=150)
            # plt.title(modelName)
            # for i in range(nlev):
            #     plt.scatter(rs, pops[:,i]/abun, s=0.5, label='i='+str(i), zorder=1)
            #
            # # plt.scatter(rs, pops[:,1]/abun, s=0.5, label='i=1', zorder=1)
            # # plt.scatter(rs, pops[:,2]/abun, s=0.5, label='i=2', zorder=1)
            # # plt.scatter(rs, pops[:,3]/abun, s=0.5, label='i=3', zorder=1)
            # # plt.scatter(rs, pops[:,4]/abun, s=0.5, label='i=4', zorder=1)
            # # plt.scatter(rs, pops[:,5]/abun, s=0.5, label='i=5', zorder=1)
            # # plt.plot(ra, lp0, c='lightgray', zorder=0)
            # # plt.plot(ra, lp1, c='lightgray', zorder=0)
            # # plt.legend()
            # plt.xscale('log')
            # plt.xlabel('r [m]')
            # plt.ylabel('fractional level populations [.]')
            #
            # plt.figure()
            # plt.title("cooling rates")
            # plt.scatter(rs, np.abs(coll_cooling_rate), s=0.5, label='collisional', zorder=1)
            # plt.scatter(rs, np.abs(rad_cooling_rate), s=0.5, label='radiative', zorder=1)
            # plt.legend()
            # plt.xscale('log')
            # plt.xlabel('r [m]')
            # plt.yscale('log')

            plt.figure()
            plt.title("cooling rates")
            plt.xscale('log')
            plt.xlabel('ρ H20 [1/m^3]')
            # plt.yscale('symlog')
            plt.yscale('log')
            plt.ylabel('Λ [J/(m^3 s)]')
            plt.scatter(abun, np.abs(coll_cooling_rate), s=0.5, label='collisional', zorder=1)
            plt.scatter(abun, np.abs(rad_cooling_rate), s=0.5, label='radiative', zorder=1)
            plt.scatter(abun, np.abs(rad_cooling_rate_shortchar), s=0.5, label='radiative sc', zorder=1)
            # plt.scatter(abun, coll_cooling_rate, s=0.5, label='collisional', zorder=1)
            # plt.scatter(abun, rad_cooling_rate, s=0.5, label='radiative', zorder=1)
            # plt.scatter(abun, rad_cooling_rate_shortchar, s=0.5, label='radiative sc', zorder=1)
            plt.legend()
            # plt.ylim(1e-27, 1e-5)

            # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

    # print(abundances)
    # print(cooling_rates)
    # print(temps)

    # plt.figure()
    # plt.title("T 100K")
    # plt.scatter(abundances[temps==100], cooling_rates[temps==100])
    # plt.xscale('log')
    # plt.xlabel('rho [todo]')
    # plt.yscale('log')
    #
    # plt.figure()
    # plt.title("T 500K")
    # plt.scatter(abundances[temps==500], cooling_rates[temps==500])
    # plt.xscale('log')
    # plt.xlabel('rho [todo]')
    # plt.yscale('log')
    #
    # plt.figure()
    # plt.title("T 2000K")
    # plt.scatter(abundances[temps==2000], cooling_rates[temps==2000])
    # plt.xscale('log')
    # plt.xlabel('rho [todo]')
    # plt.yscale('log')

    plt.show()
    # #maximal error lies mainly on the boundary (due to not exactly obeying boundary conditions) and in the regime change
    # #bound valid for both a and b benchmark
    # FEAUTRIER_AS_EXPECTED=((np.max(error_0[1:])<0.12)&(np.max(error_1[1:])<0.12))
    #
    # if not FEAUTRIER_AS_EXPECTED:
    #     print("Feautrier solver max error too large; [0]:", np.max(error_0[1:]), " [1]:", np.max(error_1[1:]))
    #
    # return (FEAUTRIER_AS_EXPECTED)
    return

def run_test (nosave=False):

    # create_model ('a')
    # run_model    ('a', nosave)

    create_models ('b')
    run_models    ('b', nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
