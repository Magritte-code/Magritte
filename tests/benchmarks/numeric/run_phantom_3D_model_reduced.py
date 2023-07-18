#this file checks whether we can run magritte on meshed data; it will not check whether the results are correct TODO: also check this maybe
#For the correctness of results, we have the analytic benchmarks

import os
import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte

path = os.path.dirname(os.path.realpath(__file__))

modeldir= path+"/../../models/"
datadir = path+"/../../data/"
result_dir= path+"/../../results/"
# model_file = os.path.join(modeldir, 'wind_00350.hdf5' )   # Resulting Magritte model
redux_file = os.path.join(modeldir, 'wind_red.hdf5' )   # Reduced Magritte model
# lamda_file = os.path.join(datadir, 'co.txt'                )   # Line data file


#something is wrong with this file/setup. We get singular eigen matrices


def run_model (nosave=True):

    modelName = f'wind_00350_red'
    # modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (redux_file)
    timer1.stop()


    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    timer2.stop()

    plt.figure()
    velocity = np.array(model.geometry.points.velocity)
    speed = np.sqrt(np.power(velocity[:,0], 2)+np.power(velocity[:,1], 2)+np.power(velocity[:,2], 2))
    inv_widths = np.array(model.lines.inverse_width)[:, 0]*np.array(model.lines.line)[0]
    plt.plot(speed*inv_widths)
    plt.show()
    # opacities=np.array(model.lines.opacity)
    # print("opacities: ",opacities)

    model.parameters.use_smoothing=False
    model.parameters.max_distance_opacity_contribution = 3.0
    # model.parameters.prune_zero_contribution_points=False
    timer3 = tools.Timer('running model')
    timer3.start()
    # model.compute_level_populations_sparse (True, 10)
    # sum_J_2f=np.array(model.lines.lineProducingSpecies[0].Jlin)

    # model = magritte.Model (modelFile)
    # model.compute_spectral_discretisation ()
    # model.compute_inverse_line_widths     ()
    # model.compute_LTE_level_populations   ()
    # # # model.COMPUTING_SVD_AND_CONDITION_NUMBER=True
    # # model.compute_level_populations_collocation(False, 1)
    # # colf = np.array(model.lines.lineProducingSpecies[0].Jlin)
    timer3.stop()

    pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), model.lines.lineProducingSpecies[0].linedata.nlev))
    abun = np.array(model.chemistry.species.abundance)[:,0]
    rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)

    np.save(f'pops_{modelName}_reduced_LTE.npy', pops)

    # (i,ra,rb,nh,tk,nm,vr,db,td,lp0,lp1) = np.loadtxt (f'{curdir}/Ratran_results/vanZadelhoff_1{a_or_b}.out', skiprows=14, unpack=True)

    # interp_0 = interp1d(0.5*(ra+rb), lp0, fill_value='extrapolate')
    # interp_1 = interp1d(0.5*(ra+rb), lp1, fill_value='extrapolate')
    #
    # error_0 = tools.relative_error(pops[:,0]/abun, interp_0(rs))
    # error_1 = tools.relative_error(pops[:,1]/abun, interp_1(rs))

    result  = f'--- Benchmark name -----------------------\n'
    result += f'{modelName                               }\n'
    result += f'--- Parameters ---------------------------\n'
    result += f'dimension = {model.parameters.dimension()}\n'
    result += f'npoints   = {model.parameters.npoints  ()}\n'
    result += f'nrays     = {model.parameters.nrays    ()}\n'
    result += f'nquads    = {model.parameters.nquads   ()}\n'
    # result += f'--- Accuracy -----------------------------\n'
    # result += f'max error in (0) = {np.max(error_0[1:])  }\n'
    # result += f'max error in (1) = {np.max(error_1[1:])  }\n'
    result += f'--- Timers -------------------------------\n'
    result += f'{timer1.print()                          }\n'
    result += f'{timer2.print()                          }\n'
    result += f'{timer3.print()                          }\n'
    result += f'------------------------------------------\n'

    print(result)

    # model.calculate_cooling_rates();

    if not nosave:
        # with open(f'{resdir}{modelName}-{timestamp}.log' ,'w') as log:
        #     log.write(result)

        plt.figure(dpi=150)
        plt.title(modelName)
        # plt.title("Constant opacity (optically thick)")
        plt.scatter(rs, pops[:,0]/abun, s=0.5, label='i=0', zorder=1)
        plt.scatter(rs, pops[:,1]/abun, s=0.5, label='i=1', zorder=1)
        # plt.scatter(rs, sum_J_2f, s=0.5, label='J_2f', zorder=1)
        # plt.plot(ra, lp0, c='lightgray', zorder=0)
        # plt.plot(ra, lp1, c='lightgray', zorder=0)
        plt.legend()
        plt.xscale('log')
        # plt.yscale('log')
        plt.xlabel('r [m]')

        plt.ylabel('$n$ []')

        # plt.ylabel('fractional level populations [.]')
        plt.show();
        # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

    # # return
    #
    # cooling_rates=np.array(model.lines.cooling_rates)
    # opacities=np.array(model.lines.opacity)
    # print("opacities: ",opacities)

    # emissivities=np.array(model.lines.emissivity)
    # Jlin=np.array(model.lines.lineProducingSpecies[0].Jlin)
    #
    # plt.figure()
    # plt.title(modelName)
    # plt.plot(rs, cooling_rates, label="cooling rate")
    # # plt.plot(rs, Jlin[:,0], label="Jlin")
    # # plt.plot(rs, opacities[:,0], label="opacity")
    # # plt.plot(rs, emissivities[:,0], label="emissivity")
    # # plt.scatter(rs, pops[:,1]/abun, s=0.5, label='i=1', zorder=1)
    # plt.legend()
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('r [m]')
    # plt.ylabel('cooling rates [J/(s*m^3)]')
    # plt.show()

    #this model has no analytic solution, so we will instead be checking whether all resulting level populations are positive
    #TODO: ADD THIS

    return #add result of test to return

def run_test (nosave=False):

    run_model    (nosave)

    # create_model ('b')
    # run_model    ('b', nosave)

    return

#this test should also be able to run on its own (not when imported)
if __name__ == '__main__':
    run_test()
