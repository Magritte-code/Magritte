#this file checks whether we can run magritte on meshed data
# It also compares the resulting image with a reference image, computed on a previous version of Magritte
# If the result is different, this could be due to the following causes:
# - changes in the mesher (fine, just create new reference image)
# - changes in the imager (might be fine, if it is just reordering the pixels (or similar))
# - changes in the solver (check the image itself whether it seems reasonable)

import os
import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte

VERSION = "0.7.2"#reference version of Magritte to compare against; should correspond to a commited intensity file

path = os.path.dirname(os.path.realpath(__file__))

modeldir= path+"/../../models/"
datadir = path+"/../../data/"
result_dir= path+"/../../results/"
# model_file = os.path.join(modeldir, 'wind_00350.hdf5' )   # Resulting Magritte model
redux_file = os.path.join(modeldir, 'wind_red.hdf5' )   # Reduced Magritte model
# lamda_file = os.path.join(datadir, 'co.txt'                )   # Line data file


def run_model (nosave=True):

    modelName = f'wind_00350_red'
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
    result += f'--- Timers -------------------------------\n'
    result += f'{timer1.print()                          }\n'
    result += f'{timer2.print()                          }\n'
    result += f'{timer3.print()                          }\n'
    result += f'------------------------------------------\n'
    
    model.compute_image_new(1.0,0.0,0.0, 256, 256)

    reference_intensity = np.load(f'{datadir}{modelName}_NLTE_intensity_magritte_{VERSION}.npy')
    reldiff = tools.relative_error(reference_intensity, np.array(model.images[0].I))

    print(result)
    print("maximum relative difference: ", np.max(reldiff))#only for a few points a significant difference can be seen
    print("mean absolute difference: ", np.mean(np.abs(reldiff)))

    # for debugging, plot of the level populations
    # if not nosave:
    #     plt.figure(dpi=150)
    #     plt.title(modelName)
    #     plt.scatter(rs, pops[:,0]/abun, s=0.5, label='i=0', zorder=1)
    #     plt.scatter(rs, pops[:,1]/abun, s=0.5, label='i=1', zorder=1)
    #     plt.legend()
    #     plt.xscale('log')
    #     plt.xlabel('r [m]')

    #     plt.ylabel('$n$ []')
    #     plt.show()

    return np.mean(np.abs(reldiff))<5e-5#actual measure difference of 2.2e-5 on the different github runners for different os. So I set the threshold a bit higher.

def run_test (nosave=False):

    run_model    (nosave)

    return

#this test should also be able to run on its own (not when imported)
if __name__ == '__main__':
    run_test()
