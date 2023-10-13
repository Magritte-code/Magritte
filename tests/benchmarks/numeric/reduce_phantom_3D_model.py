import magritte.setup  as setup                        # Model setup
import magritte.core   as magritte                     # Core functionality
import magritte.mesher as mesher                       # Mesher
import numpy           as np                           # Data structures
import vtk                                             # Reading the model
import warnings                                        # Hide warnings
warnings.filterwarnings('ignore')                      # especially for yt
import yt                                              # 3D plotting
import os

from tqdm                   import tqdm                # Progress bars
from astropy                import constants           # Unit conversions
from vtk.util.numpy_support import vtk_to_numpy        # Converting data
from scipy.spatial          import Delaunay, cKDTree   # Finding neighbors
from yt.funcs               import mylog               # To avoid yt output
mylog.setLevel(40)                                     # as error messages


def reduce_phantom():

    path = os.path.dirname(os.path.realpath(__file__))

    modeldir= path+"/../../models"
    datadir = path+"/../../data"
    model_file = os.path.join(modeldir, 'model_Phantom_3D.hdf5' )   # Resulting Magritte model
    redux_file = os.path.join(modeldir, 'wind_red' )   # Reduced Magritte model (no extension!)
    lamda_file = os.path.join(datadir, 'co.txt'                )   # Line data file

    print("Reading model")

    model = magritte.Model(model_file)

    print("Setting arrays")

    position = np.array(model.geometry.points.position)
    velocity = np.array(model.geometry.points.velocity)
    boundary = np.array(model.geometry.boundary.boundary2point)
    nCO      = np.array(model.chemistry.species.abundance)[:,0]
    nH2      = np.array(model.chemistry.species.abundance)[:,1]
    tmp      = np.array(model.thermodynamics.temperature.gas)
    trb      = np.array(model.thermodynamics.turbulence.vturb2)

    positions_reduced, nb_boundary = mesher.remesh_point_cloud(position, nCO, max_depth = 12, threshold = 4e-1, hullorder = 3)

    print("Done meshing")

    npoints          = model.parameters.npoints()
    npoints_reduced  = len(positions_reduced)

    print('npoints original =', npoints)
    print('npoints reduced  =', npoints_reduced)
    print('reduction factor =', npoints/npoints_reduced)

    # Find closest points
    corresp_points = cKDTree(position).query(positions_reduced)[1]

    # Map data
    position_reduced = positions_reduced
    velocity_reduced = velocity[corresp_points]
    nCO_reduced      = nCO     [corresp_points]
    nH2_reduced      = nH2     [corresp_points]
    tmp_reduced      = tmp     [corresp_points]
    trb_reduced      = trb     [corresp_points]

    # Extract Delaunay vertices (= Voronoi neighbors)
    delaunay = Delaunay(position_reduced)
    (indptr, indices) = delaunay.vertex_neighbor_vertices
    neighbors = [indices[indptr[k]:indptr[k+1]] for k in range(npoints_reduced)]
    nbs       = [n for sublist in neighbors for n in sublist]
    n_nbs     = [len(sublist) for sublist in neighbors]

    # Convenience arrays
    zeros = np.zeros(npoints_reduced)
    ones  = np.ones (npoints_reduced)

    model = magritte.Model ()                                        # Create model object

    model.parameters.set_model_name         (f'{redux_file}.hdf5')   # Magritte model file
    model.parameters.set_spherical_symmetry (False)                  # No spherical symmetry
    model.parameters.set_dimension          (3)                      # This is a 3D model
    model.parameters.set_npoints            (npoints_reduced)        # Number of points
    model.parameters.set_nrays              (12)                     # Number of rays
    model.parameters.set_nspecs             (3)                      # Number of species (min. 5)
    model.parameters.set_nlspecs            (1)                      # Number of line species
    model.parameters.set_nquads             (21)                     # Number of quadrature points

    model.geometry.points.position.set(position_reduced)
    model.geometry.points.velocity.set(velocity_reduced)

    model.geometry.points.  neighbors.set(  nbs)
    model.geometry.points.n_neighbors.set(n_nbs)

    model.chemistry.species.abundance = np.array((nCO_reduced, nH2_reduced, zeros)).T
    model.chemistry.species.symbol    = ['CO', 'H2', 'e-']

    model.thermodynamics.temperature.gas  .set(tmp_reduced)
    model.thermodynamics.turbulence.vturb2.set(trb_reduced)

    model.parameters.set_nboundary(nb_boundary)
    model.geometry.boundary.boundary2point.set(range(nb_boundary))

    # direction = np.array([[0,0,+1], [0,0,-1]])            # Comment out to use all directions
    # model.geometry.rays.direction.set(direction)          # Comment out to use all directions
    # model.geometry.rays.weight   .set(0.5 * np.ones(2))   # Comment out to use all directions

    model = setup.set_uniform_rays            (model)   # Uncomment to use all directions
    model = setup.set_boundary_condition_CMB  (model)
    model = setup.set_linedata_from_LAMDA_file(model, lamda_file, {'considered transitions': [0]})
    # model = setup.set_linedata_from_LAMDA_file(model, lamda_file)   # Consider all transitions
    model = setup.set_quadrature              (model)

    model.write()

#this file should also be able to run on its own (not when imported)
if __name__ == '__main__':
    reduce_phantom()
