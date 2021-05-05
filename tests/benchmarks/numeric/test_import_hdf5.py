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
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.mesher   as mesher
import magritte.core     as magritte
import h5py

from scipy.interpolate import interp1d

def run_model (a_or_b, nosave=False):

    # modelName = f'model_Jan_reduced'
    # modelName = f'model_Jan'
    # modelName = f'model_Jan_reduced_coarsened_test'
    # modelName = f'vanZadelhoff_1a_3D_mesher'
    # modelName = f'vanZadelhoff_1a_3D_mesher_naive_multigrid_1lvl_01coars'
    modelName = f'vanZadelhoff_1b_3D_mesher_naive_multigrid_1lvl_01coars'
    modelName = f'vanZadelhoff_2a_3D_mesher_no_multigrid'
    # modelName = f'model_Jan_reduced_test_1it_3coars_split2'
    # modelName = f'model_Jan_reduced_final_test_1block'
    # modelName = f'model_Jan_reduced_naive_multigrid_1lvl_04coars'
    # modelName = f'all_constant_single_ray'

    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    file = h5py.File(modelFile, 'r')

    timer2 = tools.Timer('setting model')
    timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   () #uncomment when
    nlevels=1;#should be coarsest level; misleading name
    # 2 multigrid levels, minimum 1 point remaining, 0.1 as tolerance, mgImplementation=1 (Naive,Vcycle,Wcycle)
    model.setup_multigrid(nlevels,0.1,1,20);

    timer2.stop()

    print(file['lines/lineProducingSpecies_0'].keys())
    #
    # timer3 = tools.Timer('running model')
    # timer3.start()
    # # model.compute_level_populations_multigrid(True, 20)
    # timer3.stop()

    # pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), 2))
    abun = np.array(model.chemistry.species.abundance)[:,1]
    rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)
    coarser_points = np.array(model.geometry.points.multiscale.get_current_points_in_grid())
    print(model.parameters.n_off_diag)

    final_it_finest=33
    final_it_coarsest=135
    # final_it_finest=3
    # final_it_coarsest=19

    model.restart_from_iteration(final_it_coarsest,1);
    model.interpolate_levelpops_local(1);


    nlev=2#number of level populations of the species

    pops_interpolated = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), nlev))
    pops_full = np.array(file['lines/lineProducingSpecies_0/populationlvl'+str(0)+'it'+str(final_it_finest)])
    rel_diff_pops0=abs(pops_full[:,0]-pops_interpolated[:,0])/(pops_full[:,0])
    rel_diff_pops1=abs(pops_full[:,1]-pops_interpolated[:,1])/(pops_full[:,1])

    plt.title("vanZadelhoff_1b_3D_mesher")
            # plt.scatter(rs, abs(pops_full[:,0]-pops_incomplete[:,0])/(pops_full[:,0]), s=0.5, label='i=0', zorder=1)
            # plt.scatter(rs, abs(pops_full[:,1]-pops_incomplete[:,0])/(pops_full[:,1]), s=0.5, label='i=1', zorder=1)
    # print(len(rel_diff_pops0[coarser_points]))
    # plt.scatter(rs,rel_diff_pops0, s=0.5, label='i=0', zorder=1)
    # plt.scatter(rs,rel_diff_pops1, s=0.5, label='i=1', zorder=1)
    plt.scatter(rs[coarser_points],rel_diff_pops0[coarser_points], s=0.5, label='i=0', zorder=1)
    plt.scatter(rs[coarser_points],rel_diff_pops1[coarser_points], s=0.5, label='i=1', zorder=1)

            # plt.scatter(rs,pops_full[:,1]/abun, s=0.5, label='i=1', zorder=1)
            # plt.scatter(rs,pops_incomplete[:,1]/abun, s=0.5, label='i=1', zorder=1)
            # plt.plot(ra, lp0, c='lightgray', zorder=0)
            # plt.plot(ra, lp1, c='lightgray', zorder=0)
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('r [m]')
    plt.ylabel('level populations relative difference[.]')

    plt.show()
    #
    # (i,ra,rb,nh,tk,nm,vr,db,td,lp0,lp1) = np.loadtxt (f'{curdir}/Ratran_results/vanZadelhoff_1{a_or_b}.out', skiprows=14, unpack=True)

        # with open(f'{resdir}{modelName}-{timestamp}_multigrid_{nlevels}_lvls.log' ,'w') as log:
        #     log.write(result)

        # plt.title(modelName)
        # plt.scatter(rs, pops[:,0]/abun, s=0.5, label='i=0', zorder=1)
        # plt.scatter(rs, pops[:,1]/abun, s=0.5, label='i=1', zorder=1)
        # plt.plot(ra, lp0, c='lightgray', zorder=0)
        # plt.plot(ra, lp1, c='lightgray', zorder=0)
        # plt.legend()
        # plt.xscale('log')
        # plt.xlabel('r [m]')
        # plt.ylabel('fractional level populations [.]')
        # plt.savefig(f'{resdir}{modelName}-{timestamp}.png', dpi=150)

        # with np.load('test.npy', 'r') as pops:
            # pops=np.load(f, pops)
    # last_iteration=1 #iterations from 1 till this number
    print(file['lines/lineProducingSpecies_0'].keys())

    # plotstoshow=[0,1,2,3,4,5]
    plotstoshow=[1,2,3,4,5,6,7,8,9]
    # plotstoshow=[1,7,12,50,56]

    nlev=2#41

    lvltoplot=1
    previtpops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), nlev))

    for it in plotstoshow:

        relative_differences=np.empty((0,len(coarser_points)))
        print("it",it)
        try:
            if it!=0:
                previtpops=np.array(file['lines/lineProducingSpecies_0/populationlvl'+str(lvltoplot)+'it'+str(it)])
            curritpops=np.array(file['lines/lineProducingSpecies_0/populationlvl'+str(lvltoplot)+'it'+str(it+1)])
            # print(curritpops)
            # print(previtpops)
        except Exception as e:#probably not able to read populations, because forgot to save after ng acceleration....
            print("Skipping plot:...")
            print(e)
            continue

        fig, ax = plt.subplots()
        plt.title(str(modelName+"lvl"+str(lvltoplot)+"it"+str(it)))

        print(max(abs(np.ones(len(coarser_points))-np.sum(np.transpose(np.abs(curritpops[coarser_points,:]))/abun[coarser_points], axis=0))))
        print(np.min(np.transpose(curritpops[coarser_points,:])/abun[coarser_points]))
        print("print min levelpop",np.min(np.transpose(curritpops[coarser_points,:])))

        rel_diff_pops=abs(curritpops[coarser_points,:]-previtpops[coarser_points,:])/(curritpops[coarser_points,:])

        ax.boxplot(rel_diff_pops)
        plt.yscale('log')
        plt.ylabel('Relative differences level populations')
        plt.xlabel('Levels')

        plt.show(block=False)



    lvltoplot=0
    model.geometry.points.multiscale.set_curr_coars_lvl(lvltoplot)

    coarser_points = np.array(model.geometry.points.multiscale.get_current_points_in_grid())
    previtpops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), nlev))

    for it in plotstoshow:

        relative_differences=np.empty((0,len(coarser_points)))
        print("it",it)
        try:
            if it!=0:
                previtpops=np.array(file['lines/lineProducingSpecies_0/populationlvl'+str(lvltoplot)+'it'+str(it)])
            curritpops=np.array(file['lines/lineProducingSpecies_0/populationlvl'+str(lvltoplot)+'it'+str(it+1)])
            # print(curritpops)
            # print(previtpops)
        except Exception as e:#probably not able to read populations, because forgot to save after ng acceleration....
            print("Skipping plot:...")
            print(e)
            continue

        fig, ax = plt.subplots()
        plt.title(str(modelName+"lvl"+str(lvltoplot)+"it"+str(it)))

        print(max(abs(np.ones(len(coarser_points))-np.sum(np.transpose(np.abs(curritpops[coarser_points,:]))/abun[coarser_points], axis=0))))
        print(np.min(np.transpose(curritpops[coarser_points,:])/abun[coarser_points]))
        print("print min levelpop",np.min(np.transpose(curritpops[coarser_points,:])))



        # for level in range(41):
        # rel_diff_pops=abs(np.transpose(curritpops[coarser_points,:])-np.transpose(previtpops[coarser_points,:]))/(np.transpose(curritpops[coarser_points,:]))
        rel_diff_pops=abs(curritpops[coarser_points,:]-previtpops[coarser_points,:])/(curritpops[coarser_points,:])
        # rel_diff_pops=curritpops[coarser_points,:]/abun[coarser_points, np.newaxis]
        # wrong_points=coarser_points[np.max(rel_diff_pops,axis=1)>1.0]
        # print(wrong_points)
        # for point in wrong_points:
        #     print("point",point)
        #     print("value", rel_diff_pops[point,:])
        #     print("curritpops",curritpops[point,:])
        #     print("sum curritpos",np.sum(curritpops[point,:]))
        #     print("sum abs curritpops",np.sum(np.abs(curritpops[point,:])))
        #     print("abun",abun[point])
        # print(rel_diff_pops)
        # print("min abun",np.min(abun))
        # print("max abun",np.max(abun))
        # print("divided",np.max(abun)/np.min(abun))
            # print("ref diff pops",rel_diff_pops.shape)
            # print("curritpops",curritpops.shape)
            # print("abun",abun.shape)
        #rel_diff_pops=abs(curritpops[coarser_points,line]-previtpops[coarser_points,line])/(abun[coarser_points])
        # relative_differences=np.append(relative_differences,np.array([rel_diff_pops]),axis=0)


        # previtpops=curritpops
        # print(rel_diff_pops)
            # print(relative_differences.shape)

        ax.boxplot(rel_diff_pops)
        plt.yscale('log')
        plt.ylabel('Relative differences level populations')
        plt.xlabel('Levels')

        plt.show(block=False)

        # previtpops=curritpops
    # for it in range(last_iteration):
        # previtpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it+1)])
    #previtpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it)])
    #curritpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it+1)])#note: it is kind of slow to read almost everything twice, but no memory issues and a simple implementation
        # for line in range(model.parameters.nlspec()):


    # ax.hist2d

    plt.show()

    # pops_full = np.load('vanZadelhoff_full_multigird_1_lvl.npy')
    # pops_incomplete = np.load('vanZadelhoff_incomplete_multigird_1_lvl.npy')
    # rel_diff_pops0=abs(pops_full[:,0]-pops_incomplete[:,0])/(pops_full[:,0])
    # rel_diff_pops1=abs(pops_full[:,1]-pops_incomplete[:,1])/(pops_full[:,1])


    # plt.scatter(rs, abs(pops_full[:,0]-pops_incomplete[:,0])/(pops_full[:,0]), s=0.5, label='i=0', zorder=1)
    # plt.scatter(rs, abs(pops_full[:,1]-pops_incomplete[:,0])/(pops_full[:,1]), s=0.5, label='i=1', zorder=1)
    # print(len(rel_diff_pops0[coarser_points]))

    # plt.scatter(rs[coarser_points],rel_diff_pops1[coarser_points], s=0.5, label='i=1', zorder=1)

    # plt.legend()
    # plt.xscale('log')
    # plt.xlabel('r [m]')
    # plt.ylabel('level populations relative difference[.]')
        # plt.show()
    # plt.savefig(f'{resdir}{modelName}-{timestamp}-relative_difference.png', dpi=150)

    return


def run_test (nosave=False):
# for simplicity, we only try the first model
    # create_model ('a')
    run_model    ('a', nosave)

    # create_model ('b')
    # run_model    ('b', nosave)

    return


if __name__ == '__main__':

    nosave = (len(sys.argv) > 1) and (sys.argv[1] == 'nosave')

    run_test (nosave)
