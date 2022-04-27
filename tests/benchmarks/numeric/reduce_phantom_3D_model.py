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
    model_file = os.path.join(modeldir, 'wind_00350.hdf5' )   # Resulting Magritte model
    redux_file = os.path.join(modeldir, 'wind_00350_red' )   # Reduced Magritte model (no extension!)
    lamda_file = os.path.join(datadir, 'co.txt'                )   # Line data file

    model = magritte.Model(model_file)

    position = np.array(model.geometry.points.position)
    velocity = np.array(model.geometry.points.velocity)
    boundary = np.array(model.geometry.boundary.boundary2point)
    nCO      = np.array(model.chemistry.species.abundance)[:,1]
    nH2      = np.array(model.chemistry.species.abundance)[:,2]
    tmp      = np.array(model.thermodynamics.temperature.gas)
    trb      = np.array(model.thermodynamics.turbulence.vturb2)

    delaunay = Delaunay(position)

    r_bdy = 0.99 * np.min(np.linalg.norm(position[boundary], axis=1))

    boundary_reduced = mesher.boundary_sphere(radius = r_bdy)


# %%capture
    mesher.reduce(
        meshName  = redux_file,         # Name for reduced model
        delaunay  = delaunay,           # Delaunay object of original mesh
        tracer    = nCO,                # Tracer to optimise for sampling for
        boundary  = boundary_reduced,   # Boundary of reduced mesh
        scale_max = 1.0e+99,            # Maximum scale parameter
        scale_min = 0.0e+00,            # Minimum scale parameter
        # threshold = 0.21,               # Threashold rel. diff. for coarsening
        # threshold = 0.25,               # Threashold rel. diff. for coarsening
        threshold = 0.31,               # Threashold rel. diff. for coarsening
        fmin      = 1.0,                # Don't allow refinenment
        # ftarget   = 2.15                # 10^(-1/3) for approx 10x fewer points
        # ftarget   = 3.11                # 30^(-1/3) for approx 50x fewer points
        ftarget   = 3.68                # 50^(-1/3) for approx 50x fewer points
        )


    # %%capture
    mesh = mesher.Mesh(f'{redux_file}.vtk')

    npoints          = model.parameters.npoints()
    npoints_reduced  = len(mesh.points)

    print('npoints original =', npoints)
    print('npoints reduced  =', npoints_reduced)
    print('reduction factor =', npoints/npoints_reduced)

    # Find closest points
    corresp_points = cKDTree(position).query(mesh.points)[1]

    # Map data
    position_reduced = position[corresp_points]
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
    model.parameters.set_nrays              (12)                      # Number of rays
    model.parameters.set_nspecs             (5)                      # Number of species (min. 5)
    model.parameters.set_nlspecs            (1)                      # Number of line species
    model.parameters.set_nquads             (31)                     # Number of quadrature points
    model.parameters.set_pop_prec           (1.0e-6)                 # Pops. convergence criterion

    model.geometry.points.position.set(position_reduced)
    model.geometry.points.velocity.set(velocity_reduced)

    model.geometry.points.  neighbors.set(  nbs)
    model.geometry.points.n_neighbors.set(n_nbs)

    model.chemistry.species.abundance = np.array((zeros, nCO_reduced, nH2_reduced, zeros, ones)).T
    model.chemistry.species.symbol    = ['dummy0', 'CO', 'H2', 'e-', 'dummy1']

    model.thermodynamics.temperature.gas  .set(tmp_reduced)
    model.thermodynamics.turbulence.vturb2.set(trb_reduced)

    model.parameters.set_nboundary(len(mesh.boundary))
    model.geometry.boundary.boundary2point.set(mesh.boundary)

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
