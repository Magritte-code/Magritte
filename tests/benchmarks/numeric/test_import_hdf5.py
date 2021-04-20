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
    modelName = f'model_Jan_reduced_coarsened_test'
    # modelName = f'vanZadelhoff_1a_3D_mesher'
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
    nlevels=3;#should be coarsest level; misleading name
    # 2 multigrid levels, minimum 1 point remaining, 0.1 as tolerance, mgImplementation=1 (Naive,Vcycle,Wcycle)
    model.setup_multigrid(1,nlevels,0.2,1,20);

    timer2.stop()
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
    #
    # relative_differences=np.empty((0,len(coarser_points)))
    # print(relative_differences)
    # print(relative_differences.shape)
    #
    # J_eff=np.array(file['lines/lineProducingSpecies_0/J_effit'+str(2)])
    # wrong_points=coarser_points[np.min(J_eff,axis=1)<-0.000000001]
    # print("wrong Jeff points",wrong_points)
    #
    # # in iteration 7, check whether all jeff are positive
    # J_eff=np.array(file['lines/lineProducingSpecies_0/J_effit'+str(1)])
    # wrong_points=coarser_points[np.min(J_eff,axis=1)<-0.000000001]
    # print("wrong Jeff points",wrong_points)
    # print("Jeff point 0: ",J_eff[0,:])
    # print("max sum abs min sum divided by sum", np.max((np.sum(np.abs(J_eff),axis=1)-np.sum(J_eff,axis=1))/np.abs(np.sum(J_eff,axis=1))))
    #
    # fig, ax = plt.subplots()
    # plt.title(str(modelName+" "+"Jeff1"))
    #
    # ax.boxplot((np.sum(np.abs(J_eff),axis=1)-np.sum(J_eff,axis=1))/np.abs(np.sum(J_eff,axis=1)))
    # plt.yscale('log')
    # plt.ylabel('(sum abs - sum)/sum of J_eff per point')
    # # plt.xlabel('')
    #
    # plt.show(block=False)
    #
    # J_lin=np.array(file['lines/lineProducingSpecies_0/J_linit'+str(1)])
    # fig, ax = plt.subplots()
    # plt.title(str(modelName+" "+"Jlin1"))
    #
    # wrong_points=coarser_points[np.min(J_lin,axis=1)<-0.000000001]
    # print("wrong Jlin points",wrong_points)
    # print("Jlin point 387: ",J_lin[387,:])
    #
    # ax.boxplot((np.sum(np.abs(J_lin),axis=1)-np.sum(J_lin,axis=1))/np.abs(np.sum(J_lin,axis=1)))
    # plt.yscale('log')
    # plt.ylabel('(sum abs - sum)/sum of J_lin per point')
    # # plt.xlabel('')
    #
    # plt.show(block=False)
    #
    #
    # print("len J_eff",len(J_eff))
    # print("negative jeff vals",J_eff[J_eff<-0.000001], "length",len(J_eff[J_eff<-0.000000001]))
    # print("all J_eff values", J_eff, len(J_eff))
    #
    #
    #
    # fig, ax = plt.subplots()
    # plt.title(str(modelName+" "+"Jeff1/Jlin1"))
    #
    # toplot=(np.sum(np.abs(J_eff),axis=1))/np.abs(np.sum(J_lin,axis=1))
    # print("Big errors:",coarser_points[toplot>100])
    # print("values:",toplot[toplot>100])
    # print("Errors:",coarser_points[toplot>1])
    # print("values:",toplot[toplot>1])
    # ax.boxplot(toplot)
    # plt.yscale('log')
    # plt.ylabel('(sum abs J_eff)/sum of J_lin per point')
    #
    #
    # plt.show(block=False)
    #
    #
    # print("Jlin point 0:",J_lin[0,:])
    # print("Jeff point 0:",J_eff[0,:])
    #
    # fig, ax = plt.subplots()
    # plt.title(str(modelName+" "+"abs(Jeff1)/Jlin1 of point 0"))
    #
    # ax.boxplot(np.abs(J_eff[0,:])/J_lin[0,:])
    # plt.yscale('log')
    # plt.ylabel('(abs J_eff)/J_lin of point 0')
    #
    # plt.show(block=False)
    #
    #
    # print("Jlin point 564:",J_lin[564,:])
    # print("Jeff point 564:",J_eff[564,:])
    #
    # fig, ax = plt.subplots()
    # plt.title(str(modelName+" "+"abs(Jeff1)/Jlin1 of point 564"))
    #
    # ax.boxplot(np.abs(J_eff[564,:])/J_lin[564,:])
    # plt.yscale('log')
    # plt.ylabel('(abs J_eff)/J_lin of point 564')
    #
    # plt.show(block=False)
    #
    #
    # print("Jlin point 1162:",J_lin[1162,:])
    # print("Jeff point 1162:",J_eff[1162,:])
    #
    # fig, ax = plt.subplots()
    # plt.title(str(modelName+" "+"abs(Jeff1)/Jlin1 of point 1162"))
    #
    # ax.boxplot(np.abs(J_eff[1162,:])/J_lin[1162,:])
    # plt.yscale('log')
    # plt.ylabel('(abs J_eff)/J_lin of point 564')
    #
    # plt.show(block=False)
    #
    # previtpops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), model.lines.lineProducingSpecies[0].linedata.nlev))
    # # lastit=15;
    #
    # for point_to_test in coarser_points[toplot>1]:
    #     print("Point ",point_to_test, "lies on boundary ", not model.geometry.not_on_boundary(point_to_test))
    #
    #
    # #points with negative Jeff
    # negative_jeff_points=coarser_points[np.min(J_eff,axis=1)<0]
    # print("negative jeff points:",negative_jeff_points)
    #
    # for point_to_test in negative_jeff_points:
    #     print("Point ",point_to_test, "lies on boundary ", not model.geometry.not_on_boundary(point_to_test))


    # print("in iteration 1")
    # J_eff=np.array(file['lines/lineProducingSpecies_0/J_effit'+str(1)])
    # negative_jeff_points=coarser_points[np.min(J_eff,axis=1)<0]
    # print("len negative jeff points:",len(negative_jeff_points))
    # print("negative jeff points:",negative_jeff_points)
    #
    #
    # print("in iteration 2")
    # J_eff=np.array(file['lines/lineProducingSpecies_0/J_effit'+str(2)])
    # negative_jeff_points=coarser_points[np.min(J_eff,axis=1)<0]
    # print("len negative jeff points:",len(negative_jeff_points))
    # print("negative jeff points:",negative_jeff_points)

    # print("in iteration 90")
    # J_eff=np.array(file['lines/lineProducingSpecies_0/J_effit'+str(90)])
    # negative_jeff_points=coarser_points[np.min(J_eff,axis=1)<0]
    # print("len negative jeff points:",len(negative_jeff_points))
    # print("negative jeff points:",negative_jeff_points)

    # for point_to_test in negative_jeff_points:
    #     print("Point ",point_to_test, "lies on boundary ", not model.geometry.not_on_boundary(point_to_test))

    plotstoshow=[0,1,2,3,4,5]
    # plotstoshow=[1,2,3,4,5,6,7,8,9]
    # plotstoshow=[1,7,12,50,56]

    nlev=41
    previtpops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), nlev))

    for it in plotstoshow:

        relative_differences=np.empty((0,len(coarser_points)))
        print("it",it)
        try:
            if it!=0:
                previtpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it)])
            curritpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it+1)])
            # print(curritpops)
            # print(previtpops)
        except Exception as e:#probably not able to read populations, because forgot to save after ng acceleration....
            print("Skipping plot:...")
            print(e)
            continue

        fig, ax = plt.subplots()
        plt.title(str(modelName+" "+str(it)))

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
