import numpy as np
import scipy as sp
import healpy
import re
import astroquery.lamda as lamda

from magritte.core import LineProducingSpecies, vLineProducingSpecies,            \
                          CollisionPartner, vCollisionPartner, CC, HH, KB, T_CMB, \
                          BoundaryCondition


def check_if_1D(model):
    """
    Check if the point positions are 1D for Magritte, i.e only have a non-zero x-coordinate,
    raises a ValueError if not.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object of which to check the dimension.
    """
    for (x,y,z) in np.array(model.geometry.points.position):
        if (y != 0.0) or (z != 0.0):
            raise ValueError('Not 1D (only x-coordinate can be non-zero).')
    return


def check_if_ordered(arr):
    """
    Check if an array is ordered, raise a ValueError if not.

    Parameters
    ----------
    arr : array_like
        1D array to check for ordering.
    """
    if np.all(-arr == np.sort(-arr)) or np.all(arr == np.sort(arr)):
        pass
    else:
        raise ValueError('Not 1D (x-coordinates are not ordered).')
    return


def set_Delaunay_neighbor_lists (model):
    """
    Setter for the neighbor lists for each point, assuming they are the cell centers of a Voronoi tesselation.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    if (model.parameters.dimension() == 1):
        check_if_1D     (model)
        check_if_ordered(np.array(model.geometry.points.position)[:,0])
        # For the first point
        nbs   = [1]
        n_nbs = [1]
        # For the middle points
        for p in range(1, model.parameters.npoints()-1):
            nbs  .extend([p-1, p+1])
            n_nbs.append(2)
        # For the last point
        nbs  .append(model.parameters.npoints()-2)
        n_nbs.append(1)
    elif (model.parameters.dimension() == 2):
        raise ValueError ('Dimension = 2 is not supported.')
        #points  = [[cells.x[p], cells.y[p]] for p in range(ncells)]
        ## Make a Delaulay triangulation
        #delaunay = Delaunay (points)
        ## Extract Delaunay vertices (= Voronoi neighbors)
        #(indptr, indices) = delaunay.vertex_neighbor_vertices
        #cells.neighbors   = Long2 ([Long1 (indices[indptr[k]:indptr[k+1]]) for k in #range(ncells)])
        ## Extract the number of neighbors for each point
        #cells.n_neighbors = Long1 ([len (nList) for nList in cells.neighbors])
    elif (model.parameters.dimension() == 3):
        # Make a Delaulay triangulation
        delaunay = sp.spatial.Delaunay(np.array(model.geometry.points.position))
        # Extract Delaunay vertices (= Voronoi neighbors)
        # (indptr, indices) = delaunay.vertex_neighbor_vertices
        # nbs = [indices[indptr[k]:indptr[k+1]] for k in range(model.parameters.npoints())]
        # n_nbs = [len (nb) for nb in nbs]
        (indptr, nbs) = delaunay.vertex_neighbor_vertices
        # Extract the number of neighbors for each point
        n_nbs = [indptr[k+1]-indptr[k] for k in range(model.parameters.npoints())]
    else:
        raise ValueError ('Dimension should be 1 or 3.')
    # Cast to numpy arrays of appropriate type
    model.geometry.points.  neighbors.set(  nbs)
    model.geometry.points.n_neighbors.set(n_nbs)
    # Done
    return model


def create_rotation_matrix(n):
    '''
    Helper function to create rotation matrices.
    Returns a 3x3 orthonormal matrix around n.
    '''

    nx = np.array([1.0, 0.0, 0.0])
    ny = np.array([0.0, 1.0, 0.0])

    n1  = n
    n1 /= np.linalg.norm(n1)

    if np.linalg.norm(n1-nx) < 1.0e-6:
        n2 = np.cross(ny, n1)
    else:
        n2 = np.cross(nx, n1)
    n2 /= np.linalg.norm(n2)

    n3  = np.cross(n2, n)
    n3 /= np.linalg.norm(n3)

    return np.array([n1, n2, n3]).T


def load_balance_directions(directions, comm_size):
    """
    Basic attempt to imporved load balancing between different MPI
    processes by reordering the rays such that each process gets a
    similar set of direcitons.
    """
    # Initialize
    Rs = directions

    # Precompute distances between all directions
    distances = np.arccos(np.matmul(Rs, Rs.T))

    # Get antipode for each direction
    antipodes = np.argmax(distances, axis=0)

    # Restrict to one hemisphere (antipodes will follow)
    Rs = Rs[:len(Rs)//2]

    # Precompute distances between all directions
    distances = np.arccos(np.matmul(Rs, Rs.T))

    # Initialize index lists
    inds   = list(range(len(Rs)))
    inds_o = [[] for _ in range(comm_size)]

    while len(inds) >= comm_size:
        # Take the next index
        id0 = inds[0]
        inds     .remove(id0)
        inds_o[0].append(id0)
        # Append indices of close directions to other processes
        for i in range(1, comm_size):
            idi = inds[np.argmin(distances[id0][inds])]
            inds     .remove(idi)
            inds_o[i].append(idi)

    # Append final indices
    for i, idx in enumerate(inds):
        inds_o[i].append(idx)

    # Unravel
    inds_o = [j for sl in inds_o for j in sl]

    # Add antipodes
    inds_o.extend(antipodes[inds_o].tolist())

    return directions[inds_o]


def set_uniform_rays(model, randomize=False, first_ray=np.array([1.0, 0.0, 0.0])):
    """
    Setter for rays to uniformly distributed directions.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.
    randomized : bool
        Whether or not to randomize the directions of the rays.
    first_ray: array-like
        Direction vector of the first ray in the ray list.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    # Define nrays for convenience
    nrays = model.parameters.nrays()
    if (model.parameters.dimension() == 1):
        if (model.parameters.spherical_symmetry()):
            return set_rays_spherical_symmetry (model, uniform=True)
        else:
            if (nrays != 2):
                raise ValueError ('In 1D without spherical symmetry, nrays should be 2.')
            direction = [[-1, 0, 0], [1, 0, 0]]
    elif (model.parameters.dimension() == 2):
        raise ValueError ('Dimension = 2 is not supported.')
    elif (model.parameters.dimension() == 3):
        direction  = healpy.pixelfunc.pix2vec(healpy.npix2nside(nrays), range(nrays))
        direction  = np.array(direction).transpose()
        if randomize:
            direction = sp.spatial.transform.Rotation.random().apply(direction)

        # Rotate such that the first ray is in the given direction
        R1 = create_rotation_matrix(direction[0])
        R2 = create_rotation_matrix(first_ray)
        R  = R2 @ np.linalg.inv(R1)

        direction = direction @ R.T

    else:
        raise ValueError ('Dimension should be 1 or 3.')
    # Cast to numpy arrays of appropriate type and shape
    model.geometry.rays.direction.set(direction)
    model.geometry.rays.weight   .set((1.0/nrays) * np.ones(nrays))
    # Done
    return model


def set_rays_spherical_symmetry(model, uniform=True, nextra=0, step=1):
    """
    Setter for rays in a 1D spherically symmetric model.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.
    uniform : bool
        Whether or not to use uniformly distributed rays.
    nextra : int
        Number of extra rays to add.
    step : int
        Step size used to step through the points when shooting rays (step=1 uses all points).

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    if uniform:
        # Treat every ray as "extra", since "extra" rays are added uniformly anyway
        nrays  = model.parameters.nrays()
        nextra = nrays//2-1
        rs     = []
    else:
        rs    = np.array(model.geometry.points.position)[:,0]
        # Below we assume that rs is sorted and the last element is the largest, hence sort
        rs = np.sort(rs)
        R  = rs[-1]
    
    # Add the first ray, right through the centre, i.e. along the x-axis
    Rx = [1.0]
    Ry = [0.0]
    Rz = [0.0]
    
    # Add the other rays, such that for the outer shell each ray touches another shell
    for ri in rs[0:-1:step]:
        ry = ri / R
        rx = np.sqrt(1.0 - ry**2)
        Rx.append(rx)
        Ry.append(ry)
        Rz.append(0.0)
    
    # Determine up to which angle we already have rays
    angle_max   = np.arctan(Ry[-1] / Rx[-1])
    angle_extra = (0.5*np.pi - angle_max) / nextra
    
    # Fill the remaining angular space uniformly
    for k in range(1, nextra):
        Rx.append(np.cos(angle_max + k*angle_extra))
        Ry.append(np.sin(angle_max + k*angle_extra))
        Rz.append(0.0)
    
    # Add the last ray, orthogonal to the radial direction, i.e. along the y-axis
    Rx.append(0.0)
    Ry.append(1.0)
    Rz.append(0.0)
    
    # Determine the number of rays
    assert len(Rx) == len(Ry)
    assert len(Rx) == len(Rz)
    nrays = 2*len(Rx)
    
    # Compute the weights for each ray
    # \int_{half previous}^{half next} sin \theta d\theta
    # = cos(half previous angle) - cos(half next angle)
    Wt = []
    for n in range(nrays//2):
        if   (n == 0):
            upper_x, upper_y = 0.5*(Rx[n]+Rx[n+1]), 0.5*(Ry[n]+Ry[n+1])
            lower_x, lower_y = Rx[ 0], Ry[ 0]
        elif (n == nrays//2-1):
            upper_x, upper_y = Rx[-1], Ry[-1]
            lower_x, lower_y = 0.5*(Rx[n]+Rx[n-1]), 0.5*(Ry[n]+Ry[n-1])
        else:
            upper_x, upper_y = 0.5*(Rx[n]+Rx[n+1]), 0.5*(Ry[n]+Ry[n+1])
            lower_x, lower_y = 0.5*(Rx[n]+Rx[n-1]), 0.5*(Ry[n]+Ry[n-1])
    
        Wt.append(  lower_x / np.sqrt(lower_x**2 + lower_y**2)
                  - upper_x / np.sqrt(upper_x**2 + upper_y**2) )
    
    # Append the antipodal rays
    for n in range(nrays//2):
        Rx.append(-Rx[n])
        Ry.append(-Ry[n])
        Rz.append(-Rz[n])
        Wt.append( Wt[n])
    
    # Create numpy arrays
    direction = np.array((Rx,Ry,Rz)).transpose()
    weight    = np.array(Wt)
    
    # Normalize the weights
    weight /= np.sum(weight)
    
    # Set the direction and the weights in the Magritte model
    model.geometry.rays.direction.set(direction)
    model.geometry.rays.weight   .set(weight)
    
    # Set nrays in the model
    try:
        model.parameters.set_nrays(nrays)
    except:
        raise RuntimeError(
            f"The specified number of rays in the model (nrays={model.parameters.nrays()}) does not match\n"
            f"with the specified nextra={nextra} and step={step}. Either don't specify the number of rays\n"
            f"for the model (i.e. don't invoke model.parameters.set_nrays) so this function can\n"
            f"set nrays, or specify the right number, which is {nrays} in this particular case."
        )    
    
    # Done
    return model


def set_Delaunay_boundary (model):
    """
    Setter for the boundary, assuming all points are the cell centers of a Voronoi tesselation.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    if (model.parameters.dimension() == 1):
        check_if_1D     (model)
        check_if_ordered(np.array(model.geometry.points.position)[:,0])
        model.parameters.set_nboundary(2)
        model.geometry.boundary.boundary2point.set([0, model.parameters.npoints()-1])
    elif (model.parameters.dimension() == 2):
        raise ValueError ('Dimension = 2 is not supported.')
    elif (model.parameters.dimension() == 3):
        raise ValueError ('Dimension = 3 is not supported.')
    else:
        raise ValueError ('Dimension should be 1.')
    # Done
    return model


def set_boundary_condition_zero (model):
    """
    Setter for incoming zero boundary condition at each boundary point.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    for b in range(model.parameters.nboundary()):
        model.geometry.boundary.set_boundary_condition (b, BoundaryCondition.Zero)
    # Done
    return model


def set_boundary_condition_CMB (model):
    """
    Setter for incoming CMB boundary condition at each boundary point.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    for b in range(model.parameters.nboundary()):
        model.geometry.boundary.set_boundary_condition (b, BoundaryCondition.CMB)
    model.geometry.boundary.boundary_temperature.set([T_CMB for _ in range(model.parameters.nboundary())])
    # Done
    return model


def set_boundary_condition_1D (model, T_in=T_CMB, T_out=T_CMB):
    """
    Setter for incoming black body radiation boundary condition at each boundary point.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.
    T_in : float
        Boundary temperature at the inner boundary.
    T_out : float
        Boundary temperature at the outer boundary.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    if not (model.parameters.dimension() == 1):
        raise ValueError ('These boundary conditions only work for a 1D model.')
    if not (model.parameters.nboundary() == 2):
        raise ValueError ('A 1D model must have exactly 2 boundary points.')
    else:
        # Set all boundary conditions to Thermal
        for b in range(model.parameters.nboundary()):
            model.geometry.boundary.set_boundary_condition (b, BoundaryCondition.Thermal)
        # Set inner and outer temperature
        model.geometry.boundary.boundary_temperature.set([T_in, T_out])
    # Done
    return model


def set_quadrature(model):
    """
    Setter for the quadrature roots and weights for the Gauss-Hermite
    quadrature, used for integrating over (Gaussian) line profiles.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.
    """
    # Get (Gauss-Hermite) quadrature roots and weights
    (roots, weights) = np.polynomial.hermite.hermgauss(model.parameters.nquads())
    # Normalize weights
    weights = weights / np.sum(weights)
    # Set roots and weights
    for l in range(model.parameters.nlspecs()):
        model.lines.lineProducingSpecies[l].quadrature.roots  .set(roots)
        model.lines.lineProducingSpecies[l].quadrature.weights.set(weights)
    return model


def getProperName(name):
    '''
    Return the standard name for the species.
    '''
    if name in ['e']:
        return 'e-'
    if name in ['pH2', 'oH2', 'p-H2', 'o-H2', 'PH2', 'OH2']:
        return 'H2'
    # If none of the above special cases, it should be fine
    return name


def getSpeciesNumber (species, name):
    '''
    Returns number of species given by 'name'
    '''
    # Note that there are dummy species in Magritte at places 0 and NLSPEC
    if isinstance (name, list):
        return [getSpeciesNumber (species,elem) for elem in name]
    else:
        for i in range (len (species.symbol)):
            if (species.symbol[i] == getProperName (name)):
                return i
        return 0


def extractCollisionPartner (fileName, line, species, elem):
    '''
    Returns collision partner and whether it is ortho or para (for H2)
    '''
    with open (fileName) as dataFile:
        data = dataFile.readlines ()
    partner   = re.findall (elem.replace ('+','\+')+'\s*[\+\-]?\s*([\w\+\-]+)\s*', data[line])[0]
    excess    = re.findall ('[op]\-?', partner)
    if (len (excess) > 0):
        orthoPara = re.findall ('[op]', partner)[0]
        partner   = partner.replace (excess[0],'')
    else:
        orthoPara = 'n'
    return [getSpeciesNumber (species, partner), orthoPara]


#Astroquery actually reads the number in front of the line to determine which collider is used.
# Thus the names for H2 are limited to PH2, OH2.
def extractCollisionPartnerAstroquery (partner, species):
    orthopara = 'n'
    if partner == 'PH2':
        orthopara = 'p'
    elif partner == 'OH2':
        orthopara = 'o'

    return [getSpeciesNumber (species, partner), orthopara]


#Astroquery names the temperature colums for the collision rates in the following manner: "C_ij(T=temp)"
def convertTempToColumnNameAstroquery (temp):
    return "C_ij(T="+str(temp)+")"



class LamdaFileReader ():
    """
    Reader for LAMDA line data files.
    """

    def __init__ (self ,fileName):
        """
        Set the file name of the LAMDA file.

        Parameters
        ----------
        fileName : str
            Name of the LAMDA line data file.
        """
        self.fileName = fileName

    def readColumn (self, start, nElem, columnNr, type):
        '''
        Returns a column of data as a list
        '''
        with open (self.fileName) as dataFile:
            lineNr = 0
            column = []
            for line in dataFile:
                if (lineNr >= start) and (lineNr < start+nElem):
                    if type == 'float':
                        column.append (float(line.split()[columnNr]))
                    if type == 'int':
                        column.append (int  (line.split()[columnNr]))
                    if type == 'str':
                        column.append (str  (line.split()[columnNr]))
                lineNr += 1
        return column

    def extractCollisionPartner (self, line, species, elem):
        '''
        Returns collision partner and whether it is ortho or para (for H2)
        '''
        with open (self.fileName) as dataFile:
            data = dataFile.readlines ()
        partner   = re.findall (elem.replace ('+','\+')+'\s*[\+\-]?\s*([\w\+\-]+)\s*', data[line])[0]
        excess    = re.findall ('[op]\-?', partner)
        if (len (excess) > 0):
            orthoPara = re.findall ('[op]', partner)[0]
            partner   = partner.replace (excess[0],'')
        else:
            orthoPara = 'n'
        return [getSpeciesNumber (species, partner), orthoPara]



def set_linedata_from_LAMDA_file (model, fileNames, config={}):
    """
    Set line data by reading it from a data file in LAMDA format.

    Parameters
    ----------
    model : Magritte model object
        Magritte model object to set.
    fileNames : list of strings
        List of file names for the LAMDA line data files.
    config : dict
        Optionally specify specific radiative transitions to consider.

    Returns
    -------
    out : Magritte model object
        Updated Magritte object.

    Note
    ----
    Do not use the Magritte objects linedata etc. this will kill performance. Hence, the copies.
    """
    # Make sure fileNames is a list
    if not type(fileNames) is list:
        fileNames = [fileNames]
    # Make sure a file is provides for each species
    if len(fileNames) != model.parameters.nlspecs():
        raise ValueError('Number of provided LAMDA files != nlspecs')
    # Create lineProducingSpecies objects    
    model.lines.resize_LineProducingSpecies (len(fileNames))
    # Convenient name
    species = model.chemistry.species
    # Add data for each LAMDA file
    for lspec, fileName in enumerate(fileNames):
        #new astroquery file read
        collrates, radtransitions, enlevels = lamda.parse_lamda_datafile(fileName)
        sym = enlevels.meta['molecule']# is given after ! MOLECULE
        num          = getSpeciesNumber(species, sym)
        mass = enlevels.meta['molwt']
        inverse_mass = float (1.0/mass)
        nlev = enlevels.meta['nenergylevels']
        energy = enlevels['Energy']
        weight = enlevels['Weight']
        nrad = radtransitions.meta['radtrans']
        irad = radtransitions['Upper'] #upper transition
        jrad = radtransitions['Lower'] #lower transition
        A = radtransitions['EinsteinA']

        # Old reader object
        # # Create reader for data file
        # rd = LamdaFileReader(fileName)
        # # Read radiative data
        # sym          = rd.readColumn(start= 1,      nElem=1,    columnNr=0, type='str')[0]
        # num          = getSpeciesNumber(species, sym)
        # mass         = rd.readColumn(start= 3,      nElem=1,    columnNr=0, type='float')[0]
        # inverse_mass = float (1.0 / mass)
        # nlev         = rd.readColumn(start= 5,      nElem=1,    columnNr=0, type='int')[0]
        # energy       = rd.readColumn(start= 7,      nElem=nlev, columnNr=1, type='float')
        # weight       = rd.readColumn(start= 7,      nElem=nlev, columnNr=2, type='float')
        # nrad         = rd.readColumn(start= 8+nlev, nElem=1,    columnNr=0, type='int')[0]
        # irad         = rd.readColumn(start=10+nlev, nElem=nrad, columnNr=1, type='int')
        # jrad         = rd.readColumn(start=10+nlev, nElem=nrad, columnNr=2, type='int')
        # A            = rd.readColumn(start=10+nlev, nElem=nrad, columnNr=3, type='float')

        # Change index range from [1, nlev] to [0, nlev-1]
        for k in range(nrad):
            irad[k] += -1
            jrad[k] += -1
        # Convert to SI units
        for i in range(nlev):
            # Energy from [cm^-1] to [J]
            energy[i] *= 1.0E+2*HH*CC
        # Set data
        model.lines.lineProducingSpecies[lspec].linedata.sym          = sym
        model.lines.lineProducingSpecies[lspec].linedata.num          = num
        model.lines.lineProducingSpecies[lspec].linedata.inverse_mass = inverse_mass
        model.lines.lineProducingSpecies[lspec].linedata.nlev         = nlev
        model.lines.lineProducingSpecies[lspec].linedata.energy       = energy
        model.lines.lineProducingSpecies[lspec].linedata.weight       = weight

        ncolpar = len(collrates)

        # # Start reading collisional data
        # nlr = nlev + nrad
        # # Get number of collision partners
        # ncolpar = rd.readColumn(start=11+nlr, nElem=1,  columnNr=0, type='int')[0]
        # ind     = 13 + nlr

        # Set number of collision partners
        model.lines.lineProducingSpecies[lspec].linedata.ncolpar = ncolpar
        # Create list of CollisionPartners
        model.lines.lineProducingSpecies[lspec].linedata.colpar = vCollisionPartner([CollisionPartner() for _ in range(ncolpar)])
        # Loop over the collision partners
        # for c in range(ncolpar):
        for partner, c in zip(collrates.keys(), range(ncolpar)):
            num_col_partner, orth_or_para_H2 = extractCollisionPartnerAstroquery (partner, species)
            ncol = collrates[partner].meta['ntrans']
            ntmp = collrates[partner].meta['ntemp']
            icol = collrates[partner]['Upper']
            jcol = collrates[partner]['Lower']

            # num_col_partner = rd.extractCollisionPartner(line=ind, species=species, elem=sym)[0]
            # orth_or_para_H2 = rd.extractCollisionPartner(line=ind, species=species, elem=sym)[1]
            # ncol            = rd.readColumn(start=ind+2, nElem=1,    columnNr=0,   type='int')[0]
            # ntmp            = rd.readColumn(start=ind+4, nElem=1,    columnNr=0,   type='int')[0]
            # icol            = rd.readColumn(start=ind+8, nElem=ncol, columnNr=1,   type='int')
            # jcol            = rd.readColumn(start=ind+8, nElem=ncol, columnNr=2,   type='int')
            # Change index range from [1, nlev] to [0, nlev-1]
            for k in range(ncol):
                icol[k] += -1
                jcol[k] += -1
            tmp = []
            Cd  = []

            tmp = collrates[partner].meta['temperatures']

            for temp in tmp:
                Cd.append(collrates[partner][convertTempToColumnNameAstroquery(temp)])

            # for t in range (ntmp):
            #     tmp.append (rd.readColumn(start=ind+6, nElem=1,    columnNr=t,   type='float')[0])
            #     Cd .append (rd.readColumn(start=ind+8, nElem=ncol, columnNr=3+t, type='float'))

            # Convert to SI units
            for t in range(ntmp):
                for k in range(ncol):
                    # Cd from [cm^3] to [m^3]
                    Cd[t][k] *= 1.0E-6
            # Derive excitation coefficients
            Ce = [[0.0 for _ in range(ncol)] for _ in range(ntmp)]
            for t in range(ntmp):
                for k in range(ncol):
                    i = icol[k]
                    j = jcol[k]
                    Ce[t][k] = Cd[t][k] * weight[i]/weight[j] * np.exp( -(energy[i]-energy[j]) / (KB*tmp[t]) )
            # Set data
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].num_col_partner = num_col_partner
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].orth_or_para_H2 = orth_or_para_H2
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].ncol            = ncol
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].ntmp            = ntmp
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].icol            = icol
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].jcol            = jcol
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].tmp             = tmp
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].Cd              = Cd
            model.lines.lineProducingSpecies[lspec].linedata.colpar[c].Ce              = Ce
            # Increment index
            # ind += 9 + ncol
        # Limit to the specified lines if required
        if ('considered transitions' in config) and (config['considered transitions'] is not None):
            if not isinstance(config['considered transitions'], list):
                config['considered transitions'] = [config['considered transitions']]
            if (len(config['considered transitions']) > 0):
                print('Not considering all radiative transitions on the data file but only the specified ones!')
                nrad = len (config['considered transitions'])
                irad = [irad[k] for k in config['considered transitions']]
                jrad = [jrad[k] for k in config['considered transitions']]
                A    = [   A[k] for k in config['considered transitions']]
        # Set derived quantities
        Bs        = [0.0 for _ in range(nrad)]
        Ba        = [0.0 for _ in range(nrad)]
        frequency = [0.0 for _ in range(nrad)]
        for k in range(nrad):
            i = irad[k]
            j = jrad[k]
            frequency[k] = (energy[i]-energy[j]) / HH
            Bs[k]        = A[k] * CC**2 / (2.0*HH*(frequency[k])**3)
            Ba[k]        = weight[i]/weight[j] * Bs[k]
        # Set data
        model.lines.lineProducingSpecies[lspec].linedata.nrad      = nrad
        model.lines.lineProducingSpecies[lspec].linedata.irad      = irad
        model.lines.lineProducingSpecies[lspec].linedata.jrad      = jrad
        model.lines.lineProducingSpecies[lspec].linedata.A         = A
        model.lines.lineProducingSpecies[lspec].linedata.Bs        = Bs
        model.lines.lineProducingSpecies[lspec].linedata.Ba        = Ba
        model.lines.lineProducingSpecies[lspec].linedata.frequency = frequency
    # Done
    return model
