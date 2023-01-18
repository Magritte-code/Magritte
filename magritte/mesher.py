import sys, os
import meshio
import numpy     as np
import itertools as itt
import numba

from string                  import Template
from subprocess              import Popen, PIPE
from healpy                  import pixelfunc
from scipy.spatial           import Delaunay
from scipy.spatial.transform import Rotation
from numba                   import njit


def relocate_indices(arr, p):
    for i in range(len(arr)):
        for j in range(len(arr[i])):
            if (arr[i][j] > p):
                arr[i][j] -= 1
    return arr


class Mesh:
    """
    Magritte custom mesh class (that ignores unconnected points).
    """
    def __init__(self, meshFile):
        """
        Read the meshFile using meshio and remove unconnected points.

        Parameters
        ----------
        meshFile : str
            File name of the mesh file.
        """
        self.mesh      = meshio.read(meshFile)
        self.points    = self.mesh.points
        self.tetras    = self.mesh.cells_dict['tetra']
        self.edges     = self.get_edges()
        self.neighbors = self.get_neighbors()
        self.boundary  = self.get_boundary()

        # Remove non-connected points
        non_point = self.get_non_point()
        while not (non_point is None):
            print(non_point)
            self.del_non_point(non_point)
            non_point = self.get_non_point()

        return


    def get_edges(self):
        edges = set([])
        for tetra in self.tetras:
            for (A,B) in itt.combinations(tetra,2):
                edges.add((A,B))
        return np.array(list(edges))

    def get_non_point(self):
        for (p,nl) in enumerate(self.neighbors):
            if (len(nl) == 0):
                return p

    def del_non_point(self, p):
        # remove the non-connected point
        self.points    = np.delete(self.points,    p, 0)
        self.neighbors = np.delete(self.neighbors, p, 0)
        # Change all indices, since a point got removed
        self.tetras    = relocate_indices(self.tetras,    p)
        self.edges     = relocate_indices(self.edges,     p)
        # self.neighbors = relocate_indices(self.neighbors, p)

    # Does not work properly, crashes Magritte.
    def get_neighbors(self):
        neighbors = [[] for _ in range(len(self.points))]
        for edge in self.edges:
            neighbors[edge[0]].append(edge[1])
            neighbors[edge[1]].append(edge[0])
        return neighbors

    def get_tetra_volume(self, tetra):
        a = self.points[tetra[0]]
        b = self.points[tetra[1]]
        c = self.points[tetra[2]]
        d = self.points[tetra[3]]
        return np.abs(np.dot(a-d,np.cross(b-d,c-d)))/6

    def get_tetra_volumes(self):
        return [self.get_tetra_volume(tetra) for tetra in self.tetras]

    def get_tetra_position(self, tetra):
        a = self.points[tetra[0]]
        b = self.points[tetra[1]]
        c = self.points[tetra[2]]
        d = self.points[tetra[3]]
        return 0.25*(a+b+c+d)

    def get_tetra_positions(self):
        return [self.get_tetra_position(tetra) for tetra in self.tetras]

    def get_edge_length(self, edge):
        a = self.points[edge[0]]
        b = self.points[edge[1]]
        r = b - a
        return np.sqrt(np.dot(r,r))

    def get_edge_lengths(self):
        return [self.get_edge_length(edge) for edge in self.edges]

    def get_edge_position(self, edge):
        a = self.points[edge[0]]
        b = self.points[edge[1]]
        return 0.5*(a+b)

    def get_edge_positions(self):
        return [self.get_edge_position(edge) for edge in self.edges]

    def get_boundary(self):
        boundary = set([])
        for elem in self.mesh.cells_dict['triangle']:
            for p in elem:
                boundary.add(p)
        for elem in self.mesh.cells_dict['line']:
            for p in elem:
                boundary.add(p)
        for elem in self.mesh.cells_dict['vertex']:
            for p in elem:
                boundary.add(p)
        boundary = list(boundary)
        return boundary


def run(command):
    """
    Run command in shell and continuously print its output.

    Parameters
    ----------
    command : str
        The command to run in the shell.
    """
    # Run command pipe output
    process = Popen(command, stdout=PIPE, shell=True)

    # Continuously read output and print
    while True:
        # Extract output
        line = process.stdout.readline().rstrip()
        # Break if there is no more line, print otherwise
        if not line:
            break
        else:
            print(line.decode("utf-8"))

    return


def convert_msh_to_pos(meshName, replace:bool=False):
    """
    Convert a .msh file to a .pos file.

    Parameters
    ----------
    modelName : str
        Path to the .msh file to convert.
    replace : bool
        Whether or not to remove (and thus replace) the original .msh file.
    """
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)

    # create the converision gmsh script file
    conversion_script = f'{meshName}_convert_to_pos.geo'

    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))

    # Extract the template
    with open(f'{thisFolder}/templates/convert_to_pos.template', 'r') as file:
        template = Template(file.read())

    # Fill out the template and write the conversion script
    with open(conversion_script,                                 'w') as file:
        file.write(template.substitute(FILE_NAME=meshName))

    # Run gmsh in a subprocess to convert the background mesh to the .pos format
    run(f'gmsh -0 {conversion_script}')

    # Remove the auxiliary file that is created (geo_unrolled) and the script file
    os.remove(f'{meshName}_convert_to_pos.geo_unrolled')
    os.remove(conversion_script)
    if replace:
        os.remove(f"{meshName}.msh")

    return


def boundary_cuboid (minVec, maxVec):
    """
    Retuns the gmsh script for a cuboid element that can be used as boundary.

    Parameters
    ----------
    minVec : array_like
        (x_min, y_min, z_min) vector defining the cuboid.
    maxVec : array_like
        (x_max, y_max, z_max) vector defining the cuboid.

    Returns
    -------
    out : str
        A string containing the gmsh script defining the cuboid.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))

    # Extract the template
    with open(f'{thisFolder}/templates/cuboid.template', 'r') as file:
        cuboid = Template(file.read())

    # Return the filled out template
    return cuboid.substitute(
               I     = 1,
               X_MIN = minVec[0],
               X_MAX = maxVec[0],
               Y_MIN = minVec[1],
               Y_MAX = maxVec[1],
               Z_MIN = minVec[2],
               Z_MAX = maxVec[2] )


def boundary_sphere (centre=np.zeros(3), radius=1.0):
    """
    Returns the gmsh script for a sphere element that can be used as boundary.

    Parameters
    ----------
    centre : array_like
        (x, y, z) coordinates of the centre of the sphere.
    radius : float
        Radius of the sphere.

    Returns
    -------
    out : str
        A string containing the gmsh script defining the sphere.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))

    # Extract the template
    with open(f'{thisFolder}/templates/sphere.template', 'r') as file:
        sphere = Template(file.read())

    # Return the filled out template
    return sphere.substitute(
               I      = 1,
               CX     = centre[0],
               CY     = centre[1],
               CZ     = centre[2],
               RADIUS = radius    )


def boundary_sphere_in_sphere (centre_in =np.zeros(3), radius_in =1.0,
                               centre_out=np.zeros(3), radius_out=1.0 ):
    """
    Returns the gmsh script for a sphere inside a sphere that can be used as boundary.

    Parameters
    ----------
    centre_in : array_like
        (x, y, z) coordinates of the centre of the inner sphere.
    radius_in : float
        Radius of the inner sphere.
    centre_out : array_like
        (x, y, z) coordinates of the centre of the outer sphere.
    radius_out : float
        Radius of the outer sphere.

    Returns
    -------
    out : str
        A string containing the gmsh script defining the sphere in sphere.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))

    # Extract the template
    with open(f'{thisFolder}/templates/sphere_in_sphere.template', 'r') as file:
        sphere = Template(file.read())

    # Return the filled out template
    return sphere.substitute(
               CX_IN      = centre_in[0],
               CY_IN      = centre_in[1],
               CZ_IN      = centre_in[2],
               RADIUS_IN  = radius_in,
               CX_OUT     = centre_out[0],
               CY_OUT     = centre_out[1],
               CZ_OUT     = centre_out[2],
               RADIUS_OUT = radius_out    )


def boundary_sphere_in_cuboid (centre_in=np.zeros(3), radius_in=1.0,
                               minVec   =np.zeros(3), maxVec   =np.ones(3)):
    """
    Returns the gmsh script for a sphere inside a cuboid that can be used as boundary.

    Parameters
    ----------
    centre_in : array_like
        (x, y, z) coordinates of the centre of the sphere.
    radius_in : float
        Radius of the sphere.
    minVec : array_like
        (x_min, y_min, z_min) vector defining the cuboid.
    maxVec : array_like
        (x_max, y_max, z_max) vector defining the cuboid.

    Returns
    -------
    out : str
        A string containing the gmsh script defining the sphere in cuboid.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))

    # Extrace the template
    with open(f'{thisFolder}/templates/sphere_in_cuboid.template', 'r') as file:
        sphere = Template(file.read())

    # Return the filled out template
    return sphere.substitute(
        CX_IN      = centre_in[0],
        CY_IN      = centre_in[1],
        CZ_IN      = centre_in[2],
        RADIUS_IN  = radius_in,
        X_MIN      = minVec[0],
        X_MAX      = maxVec[0],
        Y_MIN      = minVec[1],
        Y_MAX      = maxVec[1],
        Z_MIN      = minVec[2],
        Z_MAX      = maxVec[2]   )


def create_mesh_from_background(meshName, boundary, scale_min, scale_max):
    """
    Creates a mesh with element sizes defined on a background mesh.

    Parameters
    ----------
    meshName : str
        File name of the mesh file.
    boundary : str
        Gmsh script containing the definition of the boundary.
    scale_min : float
        Minimum desired element size of the mesh.
    scale_max : float
        Maximum desired element size of the mesh.
    """
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)

    # create the mesh generating gmsh script file
    meshing_script = f'{meshName}.geo'
    background     = f'{meshName}.pos'
    resulting_mesh = f'{meshName}.vtk'

    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))

    # Extract the template
    with open(f'{thisFolder}/templates/mesh_from_background.template', 'r') as file:
        template = Template(file.read())

    # Write the gmsh script
    with open(meshing_script, 'w') as file:
        file.write(template.substitute(
            BOUNDARY   = boundary,
            SCALE_MIN  = scale_min,
            SCALE_MAX  = scale_max,
            BACKGROUND = background   ))

    # run gmsh in a subprocess to generate the mesh from the background
    run(f'gmsh {meshing_script} -3 -saveall -o {resulting_mesh}')

    # Remove the script file
    os.remove(meshing_script)

    return


def create_mesh_from_function(meshName, boundary, scale_min, scale_max, scale_function):
    """
    Creates a mesh with element sizes defined by a scale function.

    Parameters
    ----------
    meshName : str
        File name of the mesh file.
    boundary : str
        Gmsh script containing the definition of the boundary.
    scale_min : float
        Minimum desired element size of the mesh.
    scale_max : float
        Maximum desired element size of the mesh.
    scale_function : str
        Formula for the function describing the desired mesh element sizes.
    """
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)

    # create the mesh generating gmsh script file
    meshing_script = f'{meshName}.geo'
    resulting_mesh = f'{meshName}.vtk'

    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))

    # Extract the template
    with open(f'{thisFolder}/templates/mesh_from_function.template', 'r') as file:
        template = Template(file.read())

    # Write the gmsh script
    with open(meshing_script, 'w') as file:
        file.write(template.substitute(
            BOUNDARY       = boundary,
            SCALE_MIN      = scale_min,
            SCALE_MAX      = scale_max,
            SCALE_FUNCTION = scale_function ))

    # Run gmsh in a subprocess to generate the mesh from the background
    run(f'gmsh {meshing_script} -3 -saveall -o {resulting_mesh}')

    # Remove the script file
    os.remove(meshing_script)

    return


def generate_background_from_1D_data(meshName, R, data):
    """
    Generate a background mesh from 1D (radial) data.

    Parameters
    ----------
    meshName : str
        File name of the mesh file.
    R : array_type
        Array containing the radial positions at which the data is given.
    data : array_type
        Array containing the data to put a background.
    """
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)

    # Map the data to a set of concentric HEALPix spheres
    nsides = 8
    sphere = np.array(pixelfunc.pix2vec(nsides, range(12*nsides**2))).transpose()
    points = []
    data_s = []
    for (i,r) in enumerate(R):
        for s in Rotation.random().apply(sphere):
            points.append(r*s)
            data_s.append(data[i])

    # Delaunay tetrahedralise the point set
    delaunay = Delaunay(points)

    # Write the mesh
    meshio.write_points_cells(
        filename   = f'{meshName}.msh',
        points     =             delaunay.points,
        cells      = {'tetra'  : delaunay.simplices},
        point_data = {'weights': np.array(data_s)}   )

    return


@njit
def get_pfs(nns):
    pfs = np.zeros(nns.shape, dtype=np.int64)
    tot = 0
    for i, n in enumerate(nns):
        tot    += n
        pfs[i]  = tot
    return pfs


@njit
def get_nbs(nbs, pfs, i):
    if (i == 0):
        return nbs[0:pfs[0]]
    else:
        return nbs[pfs[i-1]:pfs[i]]


@njit
def norm(vec_arr):
    norms = np.zeros(vec_arr.shape[0])
    for i in range(vec_arr.shape[0]):
        for j in range(vec_arr.shape[1]):
            norms[i] = norms[i] + vec_arr[i,j]**2
    return np.sqrt(norms)


@njit
def get_Gs(tracer, nbs, pfs):
    """
    Get the maximum relative difference in the tracer function for each point.
    """
    Gs = np.zeros(tracer.shape[0])
    for i in range(tracer.shape[0]):
        tracer_i = tracer[i]
        tracer_n = tracer[get_nbs(nbs, pfs, i)]
        Gs[i] = np.max(np.abs(tracer_i - tracer_n) / (tracer_i + tracer_n))
    return Gs


@njit
def get_Ls(pos, nbs, pfs):
    """
    Get the element size distributions for each point.
    """
    Ls = np.zeros(pos.shape[0])
    for i in range(pos.shape[0]):
        Ls[i] = np.mean(norm(pos[i] - pos[get_nbs(nbs, pfs, i)]))
    return Ls


@njit
def get_weights(pos, nbs, nns, tracer, threshold=0.21, fmin=1.0, ftarget=2.15):
    """
    Get the weights (desired element sizes) for each point.
    """
    pfs = get_pfs (nns)
    GGs = get_Gs  (tracer, nbs, pfs)
    LLs = get_Ls  (pos,    nbs, pfs)

    L_min    = fmin    * LLs
    L_target = ftarget * LLs

    weights = L_target
    weights[GGs > threshold] = L_min[GGs > threshold]

    return weights


def reduce(
        meshName,
        delaunay,
        tracer,
        boundary,
        scale_max = 1.0e+99,  # Maximum value of the scale parameter
        scale_min = 0.0e+00,  # Minimum value of the scale parameter
        threshold = 0.21,     # Threashold relative difference for coarsening
        fmin      = 1.0,      # Don't allow refinenment
        ftarget   = 2.15      # Approx 10^(-1/3) for approx 10 times fewer points
    ):
    """
    Reduce a model mesh for a given tracer function.

    Parameters
    ----------
    meshName : str
        File name of the mesh file.
    delaunay :
        scipy.spatial.Delaunay object for the given mesh.
    tracer : array_like
        1D array containing the tabulated tracer function.
    boundary : str
        Gmsh script containing the definition of the boundary.
    scale_max : float
        Maximum desired element size of the mesh.
    scale_min : float
        Minimum desired element size of the mesh.
    threshold : float
        Threashold value for relative difference below which we coarsen.
    fmin : float
        Multiplication factor of the mesh element size where the relative difference is below the threashold.
    ftarget : float
        Multiplication factor of the mesh element size where the relative difference is below the threashold.
    """
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)

    # Extract mesh structure from delaunay object
    (indptr, indices) = delaunay.vertex_neighbor_vertices
    neighbors = [indices[indptr[k]:indptr[k+1]] for k in range(len(delaunay.points))]
    nbs       = np.array([n for sublist in neighbors for n in sublist])
    nns       = np.array([len(sublist) for sublist in neighbors])
    pos       = delaunay.points

    # Create a background mesh (in .msh format)
    meshio.write_points_cells(
        filename   = f'{meshName}.msh',
        points     = delaunay.points,
        cells      = {'tetra'  : delaunay.simplices},
        point_data = {'weights': get_weights(pos=pos,
                                             nbs=nbs,
                                             nns=nns,
                                             tracer=tracer,
                                             threshold=threshold,
                                             fmin=fmin,
                                             ftarget=ftarget)}
    )

    # Convert .msh to .pos mesh for Gmsh
    convert_msh_to_pos (meshName=meshName, replace=True)

    # Create a new mesh from the background mesh
    create_mesh_from_background(
        meshName  = meshName,
        boundary  = boundary,
        scale_min = scale_min,
        scale_max = scale_max
    )

    return

def remesh_point_cloud(positions, data, max_depth=9, threshold= 5e-2, hullorder = 3):
    '''
    Remeshing method by comparing the maximal variation of the data against the threshold.
    Uses a recursive method with maximal depth max_depth.
    The hullorder specifies the density of the generated uniform boundary.
    '''

    new_positions = np.zeros((len(positions), 3))#should be large enough to contain the new positions
    #in the worst case, the recursive remeshing procedure will return a point cloud with size of the old points
    remesh_nb_points = 0#helper index for where to put the generated new points

    xyz_min = np.min(positions, axis=0)
    xyz_max = np.max(positions, axis=0)

    #Recursive re-mesh procedure puts the new points into the new_positions vector (with remesh_nb_points now containing the size of the re-meshed point cloud)
    remesh_nb_points = get_recursive_remesh(positions, data, 0, max_depth, threshold, new_positions, remesh_nb_points)
    print("new interior points: ", remesh_nb_points)
    #shrink positions vector to the actual ones
    new_positions.resize((remesh_nb_points,3))

    hull = create_cubic_uniform_hull(xyz_min, xyz_max, order=hullorder)
    nb_boundary = hull.shape[0]
    print("number boundary points: ", nb_boundary)
    new_positions = np.concatenate((hull, new_positions), axis = 0)

    # The boundary hull is located at the first nb_boundary positions indices of the new_positions vector
    return (new_positions, nb_boundary)

@numba.njit(cache=True)
def create_cubic_uniform_hull(xyz_min, xyz_max, order=3):
    # hull = np.zeros((1000, 3))#TODO: actually fill in correct value
    nx, ny, nz = (2**order+1, 2**order+1, 2**order+1)
    x_vector = np.linspace(xyz_min[0], xyz_max[0], nx)
    y_vector = np.linspace(xyz_min[1], xyz_max[1], ny)
    z_vector = np.linspace(xyz_min[2], xyz_max[2], nz)

    #x plane does not yet intersect with other planes
    xmin_plane = grid3D(np.array([xyz_min[0]]), y_vector, z_vector)
    xmax_plane = grid3D(np.array([xyz_max[0]]), y_vector, z_vector)

    #y plane intersects with x plane, so using reduced vectors for x coordinate
    ymin_plane = grid3D(x_vector[1:nx-1], np.array([xyz_min[1]]), z_vector)
    ymax_plane = grid3D(x_vector[1:nx-1], np.array([xyz_max[1]]), z_vector)

    #z plane also intersects with x plane
    zmin_plane = grid3D(x_vector[1:nx-1], y_vector[1:ny-1], np.array([xyz_min[2]]))
    zmax_plane = grid3D(x_vector[1:nx-1], y_vector[1:ny-1], np.array([xyz_max[2]]))

    #At the edges, the hull will contain duplicate points. These need to be removed
    # hull = np.unique(np.concatenate((xmin_plane, xmax_plane, ymin_plane, ymax_plane, zmin_plane, zmax_plane), axis = 0), axis = 0)
    hull = np.concatenate((xmin_plane, xmax_plane, ymin_plane, ymax_plane, zmin_plane, zmax_plane), axis = 0)

    # print(hull)
    return hull

#Simple function for the outer product of 1D vectors
#Allows us to create surface point cloud by inserting 2 vectors (and a single value for the last coordinate)
@numba.njit(cache=True)
def grid3D(x, y, z):
    xyz = np.empty(shape=(x.size*y.size*z.size, 3))
    idx = 0
    for k in range(x.size):
        for j in range(y.size):
            for i in range(z.size):
                xyz[idx] = [x[k], y[j], z[i]]
                idx+=1
    return xyz


@numba.njit(cache=True)
def get_recursive_remesh(positions, data, depth, max_depth, threshold, remesh_points, remesh_nb_points):
    '''
    Uses recursion to remesh a given point cloud (uses all data to determine whether to recurse on a smaller scale),
    by evaluating the total variation in the data (compared against the threshold).
    '''

    #If no data is left, no point should be added
    if len(data)==0:
        return remesh_nb_points

    # print("pos: ", positions, " data: ", data, " depth: ", depth, max_depth, "rem pos: ", remesh_points, "rem nb_points: ", remesh_nb_points)

    #We crop the position grid; this results in a grid more centered around the original points
    #This incurs a minor computational cost (compared to inheriting coordinates from the larger box), but the grid looks way nicer at the edges.
    minx = np.min(positions[:,0])
    maxx = np.max(positions[:,0])
    miny = np.min(positions[:,1])
    maxy = np.max(positions[:,1])
    minz = np.min(positions[:,2])
    maxz = np.max(positions[:,2])
    #Defining the cropped box coordinates
    min_coord = np.array([minx, miny, minz])
    max_coord = np.array([maxx, maxy, maxz])


    # print("min_coord, max_coord: ", min_coord, max_coord)

    #If the subdivision has been going on for too long (only a single point remaining or too deep recursion), we replace this box with a point
    if len(data)==1 or depth==max_depth:
        #add this point to the list
        remesh_points[remesh_nb_points, :] = (min_coord + max_coord) / 2.0
        remesh_nb_points+=1
        return remesh_nb_points

    minval, maxval = np.min(data), np.max(data)

    #If the total variation in this box lies within the defined bounds, we can stop and approximate this box by a single point
    if (maxval-minval)<threshold*(minval+maxval):
        #add this point to the list
        remesh_points[remesh_nb_points, :] = (min_coord + max_coord) / 2.0
        remesh_nb_points+=1
        return remesh_nb_points

    else:
        #go and do some more recursive investigation
        middle = (max_coord + min_coord) / 2.0

        # recursive within recursive can't be cached by numba
        # remesh_nb_points = recurse_over_octants(positions, data, middlecoord, depth, max_depth, threshold, remesh_points, remesh_nb_points, idx=0)

        #the outer bounds of the grid are implicitly enforced by subdivision
        #naming scheme: sections are in each coordinate either 'u'pper or 'l'ower
        #dividing over x
        upidx = positions[:, 0] >= middle[0]
        loidx = positions[:, 0] < middle[0]

        u_positions = positions[upidx, :]
        u_data = data[upidx]
        l_positions = positions[loidx, :]
        l_data = data[loidx]

        #dividing over y
        upidx_ux = u_positions[:, 1] >= middle[1]
        loidx_ux = u_positions[:, 1] < middle[1]

        upidx_lx = l_positions[:, 1] >= middle[1]
        loidx_lx = l_positions[:, 1] < middle[1]

        uu_positions = u_positions[upidx_ux, :]
        uu_data = u_data[upidx_ux]
        ul_positions = u_positions[loidx_ux, :]
        ul_data = u_data[loidx_ux]

        lu_positions = l_positions[upidx_lx, :]
        lu_data = l_data[upidx_lx]
        ll_positions = l_positions[loidx_lx, :]
        ll_data = l_data[loidx_lx]

        #dividing over z
        upidx_uxuy = uu_positions[:, 2] >= middle[2]
        loidx_uxuy = uu_positions[:, 2] < middle[2]
        upidx_uxly = ul_positions[:, 2] >= middle[2]
        loidx_uxly = ul_positions[:, 2] < middle[2]
        upidx_lxuy = lu_positions[:, 2] >= middle[2]
        loidx_lxuy = lu_positions[:, 2] < middle[2]
        upidx_lxly = ll_positions[:, 2] >= middle[2]
        loidx_lxly = ll_positions[:, 2] < middle[2]

        uuu_positions = uu_positions[upidx_uxuy, :]
        uuu_data = uu_data[upidx_uxuy]

        remesh_nb_points = get_recursive_remesh(uuu_positions, uuu_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

        uul_positions = uu_positions[loidx_uxuy, :]
        uul_data = uu_data[loidx_uxuy]

        remesh_nb_points = get_recursive_remesh(uul_positions, uul_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

        ulu_positions = ul_positions[upidx_uxly, :]
        ulu_data = ul_data[upidx_uxly]

        remesh_nb_points = get_recursive_remesh(ulu_positions, ulu_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

        ull_positions = ul_positions[loidx_uxly, :]
        ull_data = ul_data[loidx_uxly]

        remesh_nb_points = get_recursive_remesh(ull_positions, ull_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

        luu_positions = lu_positions[upidx_lxuy, :]
        luu_data = lu_data[upidx_lxuy]

        remesh_nb_points = get_recursive_remesh(luu_positions, luu_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

        lul_positions = lu_positions[loidx_lxuy, :]
        lul_data = lu_data[loidx_lxuy]

        remesh_nb_points = get_recursive_remesh(lul_positions, lul_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

        llu_positions = ll_positions[upidx_lxly, :]
        llu_data = ll_data[upidx_lxly]

        remesh_nb_points = get_recursive_remesh(llu_positions, llu_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

        lll_positions = ll_positions[loidx_lxly, :]
        lll_data = ll_data[loidx_lxly]

        remesh_nb_points = get_recursive_remesh(lll_positions, lll_data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)


        # delta_coord = (max_coord - min_coord) / 2.0
        # #defining coordinates for sub-boxes
        # def_mincoord = min_coord
        # def_maxcoord = max_coord
        # #We define the maximum by starting from true maximum (if we start from the middle, rounding errors will occur)
        #
        # box1_mincoord = def_mincoord + np.array([0,0,0]) * delta_coord
        # box1_maxcoord = def_maxcoord - np.array([1,1,1]) * delta_coord
        # box2_mincoord = def_mincoord + np.array([0,0,1]) * delta_coord
        # box2_maxcoord = def_maxcoord - np.array([1,1,0]) * delta_coord
        # box3_mincoord = def_mincoord + np.array([0,1,0]) * delta_coord
        # box3_maxcoord = def_maxcoord - np.array([1,0,1]) * delta_coord
        # box4_mincoord = def_mincoord + np.array([0,1,1]) * delta_coord
        # box4_maxcoord = def_maxcoord - np.array([1,0,0]) * delta_coord
        # box5_mincoord = def_mincoord + np.array([1,0,0]) * delta_coord
        # box5_maxcoord = def_maxcoord - np.array([0,1,1]) * delta_coord
        # box6_mincoord = def_mincoord + np.array([1,0,1]) * delta_coord
        # box6_maxcoord = def_maxcoord - np.array([0,1,0]) * delta_coord
        # box7_mincoord = def_mincoord + np.array([1,1,0]) * delta_coord
        # box7_maxcoord = def_maxcoord - np.array([0,0,1]) * delta_coord
        # box8_mincoord = def_mincoord + np.array([1,1,1]) * delta_coord
        # box8_maxcoord = def_maxcoord - np.array([0,0,0]) * delta_coord
        #
        # # box1indices   = np.sum((positions>=box1_mincoord) & (positions<=box1_maxcoord), axis=1) == 3
        # box1indices   = ((positions[:,0]>=box1_mincoord[:,0]) & (positions[:,0]<=box1_maxcoord[:,0])#x
        #               & (positions[:,1]>=box1_mincoord[:,1]) & (positions[:,1]<=box1_maxcoord[:,1])#y
        #               & (positions[:,2]>=box1_mincoord[:,2]) & (positions[:,2]<=box1_maxcoord[:,2]))#z
        # box1positions = positions[box1indices,:]
        # box1data      = data[box1indices]
        # remesh_nb_points = get_recursive_remesh(box1positions, box1data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
        #
        # box2indices   = np.sum((positions>=box2_mincoord) & (positions<=box2_maxcoord), axis=1) == 3
        # box2positions = positions[box2indices,:]
        # box2data      = data[box2indices]
        # remesh_nb_points = get_recursive_remesh(box2positions, box2data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
        #
        # box3indices   = np.sum((positions>=box3_mincoord) & (positions<=box3_maxcoord), axis=1) == 3
        # box3positions = positions[box3indices,:]
        # box3data      = data[box3indices]
        # remesh_nb_points = get_recursive_remesh(box3positions, box3data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
        #
        # box4indices   = np.sum((positions>=box4_mincoord) & (positions<=box4_maxcoord), axis=1) == 3
        # box4positions = positions[box4indices,:]
        # box4data      = data[box4indices]
        # remesh_nb_points = get_recursive_remesh(box4positions, box4data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
        #
        # box5indices   = np.sum((positions>=box5_mincoord) & (positions<=box5_maxcoord), axis=1) == 3
        # box5positions = positions[box5indices,:]
        # box5data      = data[box5indices]
        # remesh_nb_points = get_recursive_remesh(box5positions, box5data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
        #
        # box6indices   = np.sum((positions>=box6_mincoord) & (positions<=box6_maxcoord), axis=1) == 3
        # box6positions = positions[box6indices,:]
        # box6data      = data[box6indices]
        # remesh_nb_points = get_recursive_remesh(box6positions, box6data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
        #
        # box7indices   = np.sum((positions>=box7_mincoord) & (positions<=box7_maxcoord), axis=1) == 3
        # box7positions = positions[box7indices,:]
        # box7data      = data[box7indices]
        # remesh_nb_points = get_recursive_remesh(box7positions, box7data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
        #
        # box8indices   = np.sum((positions>=box8_mincoord) & (positions<=box8_maxcoord), axis=1) == 3
        # box8positions = positions[box8indices,:]
        # box8data      = data[box8indices]
        # remesh_nb_points = get_recursive_remesh(box8positions, box8data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)

    return remesh_nb_points


#Err, recursive procedure within recursive is maybe a bit ridiculous. But this can guarantee that my subdivision is correct
@numba.njit#(cache=True)
def recurse_over_octants(points, data, middle, depth, max_depth, threshold, remesh_points, remesh_nb_points, idx=0):
    if idx < 3:
        posidx = points[:, idx] >= middle[idx]
        negidx = points[:, idx] < middle[idx]
        pos_points = points[posidx]
        neg_points = points[negidx]
        pos_data = data[posidx]
        neg_data = data[negidx]
        remesh_nb_points = recurse_over_octants(pos_points, pos_data, middle, depth, max_depth, threshold, remesh_points, remesh_nb_points, idx+1)
        return recurse_over_octants(neg_points, neg_data, middle, depth, max_depth, threshold, remesh_points, remesh_nb_points, idx+1)
    else:
        #we have divided the points into 3D octant
        #TODO apply recursive mesh stuff
        return get_recursive_remesh(points, data, depth+1, max_depth, threshold, remesh_points, remesh_nb_points)
