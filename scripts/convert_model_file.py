import numpy as np
import h5py  as h5

from sys import argv


if len(argv) < 2:
    print('Please privide a model file to convert.')
if len(argv) > 1:
    fname_old = str(argv[1])
    fname_new = 'new_' + fname_old
if len(argv) > 2:
    fname_new = str(argv[2])


paths = {
    'geometry/points/position':
        'Geometry/Cells/position',
    'geometry/points/velocity':
        'Geometry/Cells/velocity',
    'geometry/points/neighbors':
        'Geometry/Cells/neighbors',
    'geometry/points/n_neighbors':
        'Geometry/Cells/n_neighbors',
    'geometry/boundary/boundary2point':
        'Geometry/Boundary/boundary2cell_nr',
    'geometry/rays/direction':
        'Geometry/Rays/rays',
    'geometry/rays/weight':
        'Geometry/Rays/weights',
    'chemistry/species/abundance':
        'Chemistry/Species/abundance',
    'chemistry/species/species':
        'Chemistry/Species/species',
    'thermodynamics/temperature/gas':
        'Thermodynamics/Temperature/gas',
    'thermodynamics/turbulence/vturb2':
        'Thermodynamics/Turbulence/vturb2',
    'lines/lineProducingSpecies_0/quadrature/roots':
        'Lines/LineProducingSpecies_0/Quadrature/roots',
    'lines/lineProducingSpecies_0/quadrature/weights':
        'Lines/LineProducingSpecies_0/Quadrature/weights',
    'lines/lineProducingSpecies_0/linedata/A':
        'Lines/LineProducingSpecies_0/Linedata/A',
    'lines/lineProducingSpecies_0/linedata/Ba':
        'Lines/LineProducingSpecies_0/Linedata/Ba',
    'lines/lineProducingSpecies_0/linedata/Bs':
        'Lines/LineProducingSpecies_0/Linedata/Bs',
    'lines/lineProducingSpecies_0/linedata/energy':
        'Lines/LineProducingSpecies_0/Linedata/energy',
    'lines/lineProducingSpecies_0/linedata/frequency':
        'Lines/LineProducingSpecies_0/Linedata/frequency',
    'lines/lineProducingSpecies_0/linedata/weight':
        'Lines/LineProducingSpecies_0/Linedata/weight',
    'lines/lineProducingSpecies_0/linedata/irad':
        'Lines/LineProducingSpecies_0/Linedata/irad',
    'lines/lineProducingSpecies_0/linedata/jrad':
        'Lines/LineProducingSpecies_0/Linedata/jrad',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_0/Cd':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_0/Cd',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_0/Ce':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_0/Ce',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_0/icol':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_0/icol',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_0/jcol':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_0/jcol',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_0/tmp':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_0/tmp',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_1/Cd':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_1/Cd',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_1/Ce':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_1/Ce',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_1/icol':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_1/icol',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_1/jcol':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_1/jcol',
    'lines/lineProducingSpecies_0/linedata/collisionPartner_1/tmp':
        'Lines/LineProducingSpecies_0/Linedata/CollisionPartner_0/tmp'
}


with h5.File(fname_old, 'r') as f_old:
    with h5.File(fname_new, 'w') as f_new:
        for (path_new, path_old) in paths.items():
            f_new[path_new] = np.array(f_old[path_old], dtype=f_old[path_old].dtype)
