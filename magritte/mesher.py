import sys, os
import meshio
import numpy     as np
import itertools as itt

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
    Magritte custom mesh class.
    """
    def __init__(self, meshFile):
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
    :param command: the command to run in the shell.
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


def convert_msh_to_pos(meshName, replace: bool = False):
    """
    Convert .msh to .pos background  mesh file.
    :param modelName: path to the .msh file to convert.
    :param replace: remove the original .msh file if True.
    """
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)
    # create the converision gmsh script file
    conversion_script = f'{meshName}_convert_to_pos.geo'
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))
    with open(f'{thisFolder}/templates/convert_to_pos.template', 'r') as file:
        template = Template(file.read())
    with open(conversion_script,                                 'w') as file:
        file.write(template.substitute(FILE_NAME=meshName))
    # run gmsh in a subprocess to convert the background mesh to the .pos format
    run(f'gmsh -0 {conversion_script}')
    # Remove the auxiliary file that is created (geo_unrolled) and the script file
    os.remove(f'{meshName}_convert_to_pos.geo_unrolled')
    os.remove(conversion_script)
    if replace:
        os.remove(f"{meshName}.msh")


def boundary_cuboid (minVec, maxVec):
    """
    Create the gmsh script for a cuboid element.
    :param minVec: lower cuboid vector.
    :param maxVec: upper cuboid vector.
    :return: gmsh script string.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))
    with open(f'{thisFolder}/templates/cuboid.template', 'r') as file:
        cuboid = Template(file.read())
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
    Create the gmsh script for a sphere element.
    :param centre: centre of the sphere.
    :param radius: radius of the sphere.
    :return: gmsh script string.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))
    with open(f'{thisFolder}/templates/sphere.template', 'r') as file:
        sphere = Template(file.read())
    return sphere.substitute(
               I      = 1,
               CX     = centre[0],
               CY     = centre[1],
               CZ     = centre[2],
               RADIUS = radius    )


def boundary_sphere_in_sphere (centre_in =np.zeros(3), radius_in =1.0,
                               centre_out=np.zeros(3), radius_out=1.0 ):
    """
    Create the gmsh script for a cuboid element.
    :param centre_in  : centre of the inner sphere.
    :param radius_in  : radius of the inner sphere.
    :param centre_out : centre of the outer sphere.
    :param radius_out : radius of the outer sphere.
    :return: gmsh script string.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))
    with open(f'{thisFolder}/templates/sphere_in_sphere.template', 'r') as file:
        sphere = Template(file.read())
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
    Create the gmsh script for a cuboid element.
    :param centre_in : centre of the inner sphere.
    :param radius_in : radius of the inner sphere.
    :param minVec    : lower cuboid vector.
    :param maxVec    : upper cuboid vector.
    :return: gmsh script string.
    """
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))
    with open(f'{thisFolder}/templates/sphere_in_cuboid.template', 'r') as file:
        sphere = Template(file.read())
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
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)
    # create the mesh generating gmsh script file
    meshing_script = f'{meshName}.geo'
    background     = f'{meshName}.pos'
    resulting_mesh = f'{meshName}.vtk'
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))
    with open(f'{thisFolder}/templates/mesh_from_background.template', 'r') as file:
        template = Template(file.read())
    with open(meshing_script,                                          'w') as file:
        file.write(template.substitute(
            BOUNDARY   = boundary,
            SCALE_MIN  = scale_min,
            SCALE_MAX  = scale_max,
            BACKGROUND = background   ))
    # run gmsh in a subprocess to generate the mesh from the background
    run(f'gmsh {meshing_script} -3 -saveall -o {resulting_mesh}')
    # Remove the script file
    # os.remove(meshing_script)


def create_mesh_from_function(meshName, boundary, scale_min, scale_max, scale_function):
    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)
    # create the mesh generating gmsh script file
    meshing_script = f'{meshName}.geo'
    resulting_mesh = f'{meshName}.vtk'
    # Get the path to this folder
    thisFolder = os.path.dirname(os.path.abspath(__file__))
    with open(f'{thisFolder}/templates/mesh_from_function.template', 'r') as file:
        template = Template(file.read())
    with open(meshing_script,                                        'w') as file:
        file.write(template.substitute(
            BOUNDARY       = boundary,
            SCALE_MIN      = scale_min,
            SCALE_MAX      = scale_max,
            SCALE_FUNCTION = scale_function ))
    # run gmsh in a subprocess to generate the mesh from the background
    run(f'gmsh {meshing_script} -3 -saveall -o {resulting_mesh}')
    # Remove the script file
    # os.remove(meshing_script)


def generate_background_from_1D_data(meshName, R, data):

    # Remove extension from meshName
    meshName, extension = os.path.splitext(meshName)

    nsides = 8
    sphere = np.array(pixelfunc.pix2vec(nsides, range(12*nsides**2))).transpose()
    points = []
    data_s = []

    for (i,r) in enumerate(R):
        for s in Rotation.random().apply(sphere):
            points.append(r*s)
            data_s.append(data[i])

    delaunay = Delaunay(points)

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