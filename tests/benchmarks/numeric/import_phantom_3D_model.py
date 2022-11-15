import magritte.setup as setup                # Model setup
import magritte.core  as magritte             # Core functionality
import numpy          as np                   # Data structures
import warnings                               # Hide warnings
warnings.filterwarnings('ignore')             # especially for yt
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

    #TODO maybe rename this ascii file (seems random number)

    input_file = os.path.join(modeldir, 'wind.ascii')   # Phantom snapshot
    model_file = os.path.join(modeldir, 'wind.hdf5' )   # Resulting Magritte model
    lamda_file = os.path.join(datadir, 'co.txt'                )   # Line data file

    input_link="https://owncloud.ster.kuleuven.be/index.php/s/Z7NfRBJrPeQa85g/download"
    lamda_link = "https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat"

    # %%capture
    subprocess.run(['wget', input_link, '--output-document', input_file])
    subprocess.run(['wget', lamda_link, '--output-document', lamda_file])
    # !wget $input_link --output-document $input_file
    # !wget $lamda_link --output-document $lamda_file

    #converting file to magritte model

    # Read the Phantom ascii file
    (x,y,z, h, rho, vx,vy,vz, u) = np.loadtxt(input_file, skiprows=14, usecols=(0,1,2,4,5,6,7,8,9), unpack=True)

    # Constants that can be read from ascii file
    # TODO: automate reading this is possible
    velocity_cte = 2.9784608e+06
    density_cte  = 5.9410314e-07
    energy_cte   = 8.8712277e+12

    nonzero_abundance=rho>0.0

    print("number points:", len(rho))
    print("number nonzero abundances:", len(nonzero_abundance[nonzero_abundance==True]))

    keep = np.logical_and(h>0.0, nonzero_abundance)

    # Exclude unphysical points and points with zero abundances
    x   = x  [keep]
    y   = y  [keep]
    z   = z  [keep]
    vx  = vx [keep]
    vy  = vy [keep]
    vz  = vz [keep]
    u   = u  [keep]
    rho = rho[keep]

    # x   = x  [h>0.0]
    # y   = y  [h>0.0]
    # z   = z  [h>0.0]
    # vx  = vx [h>0.0]
    # vy  = vy [h>0.0]
    # vz  = vz [h>0.0]
    # u   = u  [h>0.0]
    # rho = rho[h>0.0]

    # Extract the number of points
    npoints = len(x)

    # Convert rho (total density) to abundances
    nH2 = rho * density_cte * 1.0e+6 * constants.N_A.si.value / 2.02
    nCO = nH2 * 1.0e-4

    # Convenience arrays
    zeros = np.zeros(npoints)
    ones  = np.ones (npoints)

    position = np.array((x, y, z )).transpose()
    velocity = np.array((vx,vy,vz)).transpose()

    # Convert units
    position *= constants.au.si.value                    # Convert au to m
    velocity *= (velocity_cte / constants.c.cgs.value)   # cm/s to c fraction

    # Extract Delaunay vertices (= Voronoi neighbors)
    delaunay = Delaunay(position)
    (indptr, indices) = delaunay.vertex_neighbor_vertices
    neighbors = [indices[indptr[k]:indptr[k+1]] for k in range(npoints)]
    nbs       = [n for sublist in neighbors for n in sublist]
    n_nbs     = [len(sublist) for sublist in neighbors]

    # Compute the indices of the boundary particles of the mesh
    boundary = set([])
    for i in tqdm(range(delaunay.neighbors.shape[0])):
        m1  = (delaunay.neighbors[i] == -1)
        nm1 = np.sum(m1)
        if   (nm1 == 0):
            pass
        elif (nm1 == 1):
            for b in delaunay.simplices[i][m1]:
                boundary.add(b)
        elif (nm1 >= 2):
            for b in delaunay.simplices[i]:
                boundary.add(b)

    boundary = list(boundary)
    boundary = np.array(boundary)

    # The above calculation turned out to be unsatisfactory.
    # Since the outer boundary is assumed to be a sphere,
    # we add all points which fall inside the boundary defined above.
    b_nms = np.linalg.norm(position[boundary], axis=1)
    p_nms = np.linalg.norm(position,           axis=1)
    boundary = np.array([i[0] for i in np.argwhere(p_nms >= np.min(b_nms))])

    # Derive temperature from internal energy (assuming adiabatic heating/cooling)
    gamma = 1.2
    mu    = 2.381
    tmp   = mu * (gamma-1.0) * u * energy_cte * 1.00784 * (units.erg/units.g * constants.u/constants.k_B).to(units.K).value

    # Cut-off temperatures below 2.725 K
    tmp[tmp<2.725] = 2.725

    # Define turbulence at 150 m/s
    trb = (150.0/constants.c.si.value)**2 * ones


    ##creating model

    model = magritte.Model ()                              # Create model object

    model.parameters.set_model_name         (model_file)   # Magritte model file
    model.parameters.set_spherical_symmetry (False)        # No spherical symmetry
    model.parameters.set_dimension          (3)            # This is a 3D model
    model.parameters.set_npoints            (npoints)      # Number of points
    model.parameters.set_nrays              (12)            # Number of rays
    model.parameters.set_nspecs             (5)            # Number of species (min. 5)
    model.parameters.set_nlspecs            (1)            # Number of line species
    model.parameters.set_nquads             (31)           # Number of quadrature points

    model.geometry.points.position.set(position)
    model.geometry.points.velocity.set(velocity)

    model.geometry.points.  neighbors.set(  nbs)
    model.geometry.points.n_neighbors.set(n_nbs)

    model.chemistry.species.abundance = np.array((zeros, nCO, nH2, zeros, ones)).T
    model.chemistry.species.symbol    = ['dummy0', 'CO', 'H2', 'e-', 'dummy1']

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
