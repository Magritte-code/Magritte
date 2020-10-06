import numpy as np
import scipy as sp
import healpy
import re

from magritte.tools import dtypeSize, dtypeReal
from magritte.core  import LineProducingSpecies, vLineProducingSpecies, \
                           CollisionPartner, vCollisionPartner, CC, HH, KB


def check_if_1D(model):
    """
    Check if the point positions are 1D for Magritte,
    i.e only has a non-zero x-coordinate.
    """
    for (x,y,z) in model.geometry.points.position:
        if (y != 0.0) or (z != 0.0):
            raise ValueError('Not 1D (only x-coordinate can be non-zero).')
    return


def check_if_ordered(arr):
    if np.all(-arr == np.sort(-arr)) or np.all(arr == np.sort(arr)):
        pass
    else:
        raise ValueError('Not 1D (x-coordinates are not ordered).')
    return


def set_Delaunay_neighbor_lists (model):
    """
    Setter for the neighbor lists for each point, assuming
    they are the cell centers of a Voronoi tesselation.
    """
    if (model.parameters.dimension() == 1):
        check_if_1D     (model)
        check_if_ordered(model.geometry.points.position[:,0])
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
        delaunay = sp.spatial.Delaunay(model.geometry.points.position)
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
    model.geometry.points.  neighbors = np.array(nbs,   dtypeSize)
    model.geometry.points.n_neighbors = np.array(n_nbs, dtypeSize)
    # Done
    return model


def set_uniform_rays(model, randomize=False):
    """
    Setter for rays to uniformly distributed directions.
    """
    # Define nrays for convenience
    nrays = model.parameters.nrays()
    if (model.parameters.dimension() == 1):
        if (nrays != 2):
            raise ValueError ('In 1D, nrays should always be 2.')
        direction = np.array([[-1,0,0], [1,0,0]], dtypeReal)
    elif (model.parameters.dimension() == 2):
        raise ValueError ('Dimension = 2 is not supported.')
    elif (model.parameters.dimension() == 3):
        direction  = healpy.pixelfunc.pix2vec(healpy.npix2nside(nrays), range(nrays))
        direction  = np.array(direction, dtypeReal).transpose()
        if randomize:
            direction = sp.spatial.transform.Rotation.random().apply(rays)
    else:
        raise ValueError ('Dimension should be 1 or 3.')
    # Cast to numpy arrays of appropriate type and shape
    model.geometry.rays.direction = direction
    model.geometry.rays.weight    = (1.0/nrays) * np.ones(nrays, dtypeReal)
    # Done
    return model


def set_Delaunay_boundary (model):
    """
    Setter for the boundary, assuming all points points
    the cell centers of a Voronoi tesselation.
    """
    if (model.parameters.dimension() == 1):
        check_if_1D     (model)
        check_if_ordered(model.geometry.points.position[:,0])
        model.geometry.boundary.boundary2point = np.array([0, model.parameters.npoints()-1], dtypeSize)
    elif (model.parameters.dimension() == 2):
        raise ValueError ('Dimension = 2 is not supported.')
    elif (model.parameters.dimension() == 3):
        raise ValueError ('Dimension = 3 is not supported.')
    else:
        raise ValueError ('Dimension should be 1 or 3.')
    # Done
    return model


def Gauss_Hermite_roots(n):
    """
    Returns the roots for Gauss-Hermite quadrature.
    """
    coeffs  = [0.0 for _ in range(n)]
    coeffs += [1.0]
    roots = np.polynomial.hermite.hermroots(coeffs)
    return np.array(roots, dtypeReal)


def Gauss_Hermite_weights(n):
    """
    Returns the weights for Gauss-Hermite quadrature.
    """
    coeffs  = [0.0 for _ in range(n-1)]
    coeffs += [1.0]
    H = np.polynomial.hermite.hermval(Gauss_Hermite_roots(n), coeffs)
    weights = 2**(n-1) * np.math.factorial(n) / (n*H)**2
    return np.array(weights, dtypeReal)


def set_quadrature(model):
    """
    Setter for the quadrature roots and weights for the Gauss-Hermite
    quadrature, used for integrating over (Gaussian) line profiles.
    """
    for l in range(model.parameters.nlspecs()):
        model.lines.lineProducingSpecies[l].quadrature.roots   = Gauss_Hermite_roots   (model.parameters.nquads())
        model.lines.lineProducingSpecies[l].quadrature.weights = Gauss_Hermite_weights (model.parameters.nquads())
    return model


def get_proper_name(name):
    '''
    Return the standard name for the species
    '''
    if name in ['e']:
        return 'e-'
    if name in ['pH2', 'oH2', 'p-H2', 'o-H2']:
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


class LamdaFileReader ():

    def __init__ (self ,fileName):
        """
        Constructor
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



def set_linedata_from_LAMDA_file (model, fileName, config={}):
    """
    Set line data by reading it from a data file in LAMDA format.
    Note: Do not use the Magritte objects ld etc. this will kill performance. Hence, the copies.
    """
    # Create lineProducingSpecies object
    model.lines.lineProducingSpecies = vLineProducingSpecies([LineProducingSpecies()])
    # Convenient name
    species = model.chemistry.species
    # Create reader for data file
    rd = LamdaFileReader(fileName)
    # Read radiative data
    sym          = rd.readColumn(start= 1,      nElem=1,    columnNr=0, type='str')[0]
    num          = getSpeciesNumber(species, sym)
    mass         = rd.readColumn(start= 3,      nElem=1,    columnNr=0, type='float')[0]
    inverse_mass = float (1.0 / mass)
    nlev         = rd.readColumn(start= 5,      nElem=1,    columnNr=0, type='int')[0]
    energy       = rd.readColumn(start= 7,      nElem=nlev, columnNr=1, type='float')
    weight       = rd.readColumn(start= 7,      nElem=nlev, columnNr=2, type='float')
    nrad         = rd.readColumn(start= 8+nlev, nElem=1,    columnNr=0, type='int')[0]
    irad         = rd.readColumn(start=10+nlev, nElem=nrad, columnNr=1, type='int')
    jrad         = rd.readColumn(start=10+nlev, nElem=nrad, columnNr=2, type='int')
    A            = rd.readColumn(start=10+nlev, nElem=nrad, columnNr=3, type='float')
    # Change index range from [1, nlev] to [0, nlev-1]
    for k in range(nrad):
        irad[k] += -1
        jrad[k] += -1
    # Convert to SI units
    for i in range(nlev):
        # Energy from [cm^-1] to [J]
        energy[i] *= 1.0E+2*HH*CC
    # Set data
    model.lines.lineProducingSpecies[0].linedata.num          = num
    model.lines.lineProducingSpecies[0].linedata.inverse_mass = inverse_mass
    model.lines.lineProducingSpecies[0].linedata.nlev         = nlev
    model.lines.lineProducingSpecies[0].linedata.energy       = energy
    model.lines.lineProducingSpecies[0].linedata.weight       = weight
    # Start reading collisional data
    nlr = nlev + nrad
    # Get number of collision partners
    ncolpar = rd.readColumn(start=11+nlr, nElem=1,  columnNr=0, type='int')[0]
    ind     = 13 + nlr
    # Set number of collision partners
    model.lines.lineProducingSpecies[0].linedata.ncolpar = ncolpar
    # Create list of CollisionPartners
    model.lines.lineProducingSpecies[0].linedata.colpar = vCollisionPartner([CollisionPartner() for _ in range(ncolpar)])
    # Loop over the collision partners
    for c in range(ncolpar):
        num_col_partner = rd.extractCollisionPartner(line=ind, species=species, elem=sym)[0]
        orth_or_para_H2 = rd.extractCollisionPartner(line=ind, species=species, elem=sym)[1]
        ncol            = rd.readColumn(start=ind+2, nElem=1,    columnNr=0,   type='int')[0]
        ntmp            = rd.readColumn(start=ind+4, nElem=1,    columnNr=0,   type='int')[0]
        icol            = rd.readColumn(start=ind+8, nElem=ncol, columnNr=1,   type='int')
        jcol            = rd.readColumn(start=ind+8, nElem=ncol, columnNr=2,   type='int')
        # Change index range from [1, nlev] to [0, nlev-1]
        for k in range(ncol):
            icol[k] += -1
            jcol[k] += -1
        tmp = []
        Cd  = []
        for t in range (ntmp):
            tmp.append (rd.readColumn(start=ind+6, nElem=1,    columnNr=t,   type='float')[0])
            Cd .append (rd.readColumn(start=ind+8, nElem=ncol, columnNr=3+t, type='float'))
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
        model.lines.lineProducingSpecies[0].linedata.colpar[c].num_col_partner = num_col_partner
        model.lines.lineProducingSpecies[0].linedata.colpar[c].orth_or_para_H2 = orth_or_para_H2
        model.lines.lineProducingSpecies[0].linedata.colpar[c].ncol            = ncol
        model.lines.lineProducingSpecies[0].linedata.colpar[c].ntmp            = ntmp
        model.lines.lineProducingSpecies[0].linedata.colpar[c].icol            = icol
        model.lines.lineProducingSpecies[0].linedata.colpar[c].jcol            = jcol
        model.lines.lineProducingSpecies[0].linedata.colpar[c].tmp             = tmp
        model.lines.lineProducingSpecies[0].linedata.colpar[c].Cd              = Cd
        model.lines.lineProducingSpecies[0].linedata.colpar[c].Ce              = Ce
        # Increment index
        ind += 9 + ncol
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
    model.lines.lineProducingSpecies[0].linedata.nrad      = nrad
    model.lines.lineProducingSpecies[0].linedata.irad      = irad
    model.lines.lineProducingSpecies[0].linedata.jrad      = jrad
    model.lines.lineProducingSpecies[0].linedata.A         = A
    model.lines.lineProducingSpecies[0].linedata.Bs        = Bs
    model.lines.lineProducingSpecies[0].linedata.Ba        = Ba
    model.lines.lineProducingSpecies[0].linedata.frequency = frequency
    # Done
    return model
