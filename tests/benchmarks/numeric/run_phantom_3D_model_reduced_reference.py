#this python script will generate a reference image for the reduced 3D phantom model
#only run after major changes, and checking whether the result seems fine
# NOTE: the image should then be manually commited to the repository

import os
import numpy             as np
import scipy             as sp
import matplotlib.pyplot as plt
import magritte.tools    as tools
import magritte.setup    as setup
import magritte.core     as magritte
from importlib.metadata import version

path = os.path.dirname(os.path.realpath(__file__))

modeldir= path+"/../../models/"
datadir = path+"/../../data/"
result_dir= path+"/../../results/"
# model_file = os.path.join(modeldir, 'wind_00350.hdf5' )   # Resulting Magritte model
redux_file = os.path.join(modeldir, 'wind_red.hdf5' )   # Reduced Magritte model
# lamda_file = os.path.join(datadir, 'co.txt'                )   # Line data file

VERSION = version('magritte')

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

    timer3 = tools.Timer('running model')
    timer3.start()
    model.compute_level_populations_sparse (True, 10)
    timer3.stop()

    pops = np.array(model.lines.lineProducingSpecies[0].population).reshape((model.parameters.npoints(), model.lines.lineProducingSpecies[0].linedata.nlev))
    abun = np.array(model.chemistry.species.abundance)[:,0]
    rs   = np.linalg.norm(np.array(model.geometry.points.position), axis=1)

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

    reference_intensity = np.array(model.images[0].I)
    np.save(f'{datadir}{modelName}_NLTE_intensity_magritte_{VERSION}.npy', reference_intensity)

    print(result)

    return

def run_test (nosave=False):

    run_model    (nosave)

    return

#this test should also be able to run on its own (not when imported)
if __name__ == '__main__':
    run_test()
