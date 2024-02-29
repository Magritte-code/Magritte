import magritte.setup as setup                # Model setup
import magritte.core  as magritte             # Core functionality
import numpy          as np                   # Data structures
import warnings                               # Hide warnings
warnings.filterwarnings('ignore')             # especially for yt
import plons                                  # Phantom import
import yt                                     # 3D plotting
import os
import subprocess

from tqdm          import tqdm                # Progress bars
from astropy       import constants, units    # Unit conversions
from scipy.spatial import Delaunay, cKDTree   # Finding neighbors
from yt.funcs      import mylog               # To avoid yt output
mylog.setLevel(40)                            # as error messages


def import_phantom():

    path = os.path.dirname(os.path.realpath(__file__))

    modeldir= path+"/../../models"
    datadir = path+"/../../data"

    dump_file  = os.path.join(datadir, 'model_Phantom_3D'     )   # Phantom full dump (snapshot)
    setup_file = os.path.join(datadir, 'wind.setup'           )   # Phantom setup file
    input_file = os.path.join(datadir, 'wind.in'              )   # Phantom input file
    model_file = os.path.join(modeldir, 'model_Phantom_3D.hdf5')   # Resulting Magritte model
    lamda_file = os.path.join(datadir, 'co.txt'               )   # Line data file

    dump_link  = "https://github.com/Ensor-code/phantom-models/raw/main/Malfait+2021/v05e50/wind_v05e50?download="
    setup_link = "https://raw.githubusercontent.com/Ensor-code/phantom-models/main/Malfait%2B2021/v05e50/wind.setup"
    input_link = "https://raw.githubusercontent.com/Ensor-code/phantom-models/main/Malfait%2B2021/v05e50/wind.in"
    lamda_link = "https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat"

    # %%capture
    subprocess.run(['wget', dump_link,  '--output-document', dump_file ])
    subprocess.run(['wget', setup_link, '--output-document', setup_file])
    subprocess.run(['wget', input_link, '--output-document', input_file])
    subprocess.run(['wget', lamda_link, '--output-document', lamda_file]) 
    # !wget $dump_link  --output-document $dump_file
    # !wget $setup_link --output-document $setup_file
    # !wget $input_link --output-document $input_file
    # !wget $lamda_link --output-document $lamda_file

    # Loading the data with plons
    setupData = plons.LoadSetup(datadir, "wind")
    dumpData  = plons.LoadFullDump(dump_file, setupData)

    position = dumpData["position"]*1e-2      # position vectors       [cm   -> m]
    velocity = dumpData["velocity"]*1e3      # velocity vectors        [km/s -> m/s]
    velocity = velocity/constants.c.si.value # velocity vectors        [m/s  -> 1/c]
    rho      = dumpData["rho"]               # density                 [g/cm^3]
    u        = dumpData["u"]                 # internal energy density [erg/g]
    tmp      = dumpData["Tgas"]              # temperature             [K]
    tmp[tmp<2.725] = 2.725                   # Cut-off temperatures below 2.725 K

    # Extract the number of points
    npoints = len(rho)

    # Convenience arrays
    zeros = np.zeros(npoints)
    ones  = np.ones (npoints)

    # Convert rho (total density) to abundances
    nH2 = rho * 1.0e+6 * constants.N_A.si.value / 2.02
    nCO = nH2 * 1.0e-4

    # Define turbulence at 150 m/s
    trb = (150.0/constants.c.si.value)**2 * ones

    # Extract Delaunay vertices (= Voronoi neighbors)
    delaunay = Delaunay(position)
    (indptr, indices) = delaunay.vertex_neighbor_vertices
    neighbors = [indices[indptr[k]:indptr[k+1]] for k in range(npoints)]
    nbs       = [n for sublist in neighbors for n in sublist]
    n_nbs     = [len(sublist) for sublist in neighbors]

    # Compute the indices of the boundary particles of the mesh, extracted from the Delaunay vertices
    boundary = set([])
    for i in tqdm(range(delaunay.neighbors.shape[0])):
        for k in range(4):
            if (delaunay.neighbors[i][k] == -1):
                nk1,nk2,nk3 = (k+1)%4, (k+2)%4, (k+3)%4 
                boundary.add(delaunay.simplices[i][nk1])
                boundary.add(delaunay.simplices[i][nk2])
                boundary.add(delaunay.simplices[i][nk3])
                
    boundary = list(boundary)
    boundary = np.array(boundary)

    # The above calculation turned out to be unsatisfactory.
    # Since the outer boundary is assumed to be a sphere,
    # we add all points which fall inside the boundary defined above.
    b_nms = np.linalg.norm(position[boundary], axis=1)
    p_nms = np.linalg.norm(position,           axis=1)
    boundary = np.array([i[0] for i in np.argwhere(p_nms >= np.min(b_nms))])

    ##creating model

    model = magritte.Model ()                              # Create model object

    model.parameters.set_model_name         (model_file)   # Magritte model file
    model.parameters.set_spherical_symmetry (False)        # No spherical symmetry
    model.parameters.set_dimension          (3)            # This is a 3D model
    model.parameters.set_npoints            (npoints)      # Number of points
    model.parameters.set_nrays              (12)            # Number of rays
    model.parameters.set_nspecs             (3)            # Number of species (min. 5)
    model.parameters.set_nlspecs            (1)            # Number of line species
    model.parameters.set_nquads             (31)           # Number of quadrature points

    model.geometry.points.position.set(position)
    model.geometry.points.velocity.set(velocity)

    model.geometry.points.  neighbors.set(  nbs)
    model.geometry.points.n_neighbors.set(n_nbs)

    model.chemistry.species.abundance = np.array((nCO, nH2, zeros)).T
    model.chemistry.species.symbol    = ['CO', 'H2', 'e-']

    model.thermodynamics.temperature.gas  .set(tmp)
    model.thermodynamics.turbulence.vturb2.set(trb)

    model.parameters.set_nboundary(boundary.shape[0])
    model.geometry.boundary.boundary2point.set(boundary)

    model.parameters.set_nboundary(boundary.shape[0])
    model.geometry.boundary.boundary2point.set(boundary)

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
    import_phantom()
