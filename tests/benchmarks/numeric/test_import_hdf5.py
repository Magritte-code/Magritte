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

    modelName = f'model_Jan_reduced'
    modelFile = f'{moddir}{modelName}.hdf5'
    timestamp = tools.timestamp()

    timer1 = tools.Timer('reading model')
    timer1.start()
    model = magritte.Model (modelFile)
    timer1.stop()

    file = h5py.File(modelFile, 'r')

    # timer2 = tools.Timer('setting model')
    # timer2.start()
    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.compute_LTE_level_populations   ()
    # nlevels=4;#should be coarsest level; misleading name
    # #2 multigrid levels, minimum 1 point remaining, 0.1 as tolerance, mgImplementation=1 (Naive,Vcycle,Wcycle)
    # model.setup_multigrid(1,nlevels,0.1,1,20);

    # timer2.stop()
    #
    # timer3 = tools.Timer('running model')
    # timer3.start()
    # # model.compute_level_populations_multigrid(True, 20)
    # timer3.stop()

    # pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), 2))
    abun = np.array(model.chemistry.species.abundance)[:,1]
    rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)
    coarser_points = np.array(model.geometry.points.multiscale.get_current_points_in_grid())
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
    last_iteration=1 #iterations from 1 till this number
    print(file['lines/lineProducingSpecies_0'].keys())

    relative_differences=np.empty((0,len(coarser_points)))
    print(relative_differences)
    print(relative_differences.shape)

    fig, ax = plt.subplots()
    plt.title(modelName)

    #previtpops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), model.lines.lineProducingSpecies[0].linedata.nlev))

    it=1
    # for it in range(last_iteration):
        # previtpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it+1)])
    previtpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it)])
    curritpops=np.array(file['lines/lineProducingSpecies_0/populationit'+str(it+1)])#note: it is kind of slow to read almost everything twice, but no memory issues and a simple implementation
        # for line in range(model.parameters.nlspec()):
    print(curritpops)
    print(previtpops)
    for line in range(41):
        rel_diff_pops=abs(curritpops[coarser_points,line]-previtpops[coarser_points,line])/(curritpops[coarser_points,line])
        #rel_diff_pops=abs(curritpops[coarser_points,line]-previtpops[coarser_points,line])/(abun[coarser_points])
        relative_differences=np.append(relative_differences,np.array([rel_diff_pops]),axis=0)

        # previtpops=curritpops
        # print(rel_diff_pops)
        print(relative_differences.shape)
    # ax.hist2d
    ax.boxplot(np.transpose(relative_differences))
    plt.yscale('log')
    plt.ylabel('Relative differences level populations')
    plt.xlabel('Levels')

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
