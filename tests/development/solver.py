import numpy as np

from numba         import njit, jit
from scipy.spatial import cKDTree
from scipy         import sparse




class Solver:
    
    def __init__():




matsize = (nboundary + npoints) * nfreqs * nrays








# @njit
def get_expansion():
    """
    Setup the sparse matrix for the kernel expansion.
    """
    # initialize
    data = np.zeros(datsize_exp, dtype=np.float64)
    id_i = np.zeros(datsize_exp, dtype=np.int64)
    id_j = np.zeros(datsize_exp, dtype=np.int64)
    # set index
    index = 0
        
    # compute
    
    # For all boundary points
    for b2 in range(nboundary):
        p2 = boundary2point[b2]
        
        for f2 in range(nfreqs):
            for r2 in range(nrays):    

                # For all (general) points (close to the considered boundary point p1)
                for i1 in range(len(ind_all[p2])):
                    p1 = ind_all[p2][i1]
                    
                    for f1 in range(nfreqs):
                        for r1 in range(nrays):
            
                            kk = k(r1, p1, f1, r2, p2, f2)
                            j1 = r1 + nrays*(f1 + nfreqs*p1) 
                            j2 = r2 + nrays*(f2 + nfreqs*b2)
                
                            data[index] = kk
                            id_i[index] = j1 
                            id_j[index] = j2
                            index += 1    

    # For all (general) points
    for p1 in range(npoints):
        i1 = p1 + nboundary
        
        for f1 in range(nfreqs):
            for r1 in range(nrays):
            
                # For all (general) points (close to the considered general point p1)
                for i2 in range(min([len(ind_all[p1]), npoints])):
                    p2 = ind_all[p1][i2]
                    
                    for f2 in range(nfreqs):
                        for r2 in range(nrays):
                            
                            Lk = L_k(r1, p1, f1, r2, p2, f2)
                            j1 = r1 + nrays*(f1 + nfreqs* p1             ) 
                            j2 = r2 + nrays*(f2 + nfreqs*(p2 + nboundary))
                            
                            data[index] = Lk
                            id_i[index] = j1
                            id_j[index] = j2
                            index += 1

    return (data, (id_i, id_j))


# @njit
def get_covariance():
    # initialize
    data = np.zeros(datsize_cov, dtype=np.float64)
    id_i = np.zeros(datsize_cov, dtype=np.int64)
    id_j = np.zeros(datsize_cov, dtype=np.int64)
    # set index
    index = 0
    
    # compute     
    
    # For all boundary points
    for b1 in range(nboundary):
        p1 = boundary2point[b1]

        for f1 in range(nfreqs):
            for r1 in range(nrays):
                
                # For all boundary points (close to the considered boundary point p1)
                for i2 in range(min([len(ind_bdy[b1]), nboundary])):
                    b2 = ind_bdy[b1][i2]
                    p2 = boundary2point[b2]
                    
                    for f2 in range(nfreqs):
                        for r2 in range(nrays):
            
                            kk = k(r1, p1, f1, r2, p2, f2)
                            j1 = r1 + nrays*(f1 + nfreqs*b1) 
                            j2 = r2 + nrays*(f2 + nfreqs*b2)
                
                            data[index] = kk
                            id_i[index] = j1 
                            id_j[index] = j2
                            index += 1
                            
#                             print ('kk  ', p1, r1, ';', p2, r2, ' ', kk)
 
                # For all (general) points (close to the considered boundary point p1)
                for i2 in range(min([len(ind_all[p1]), npoints])):
                    p2 = ind_all[p1][i2]
                    
                    for f2 in range(nfreqs):
                        for r2 in range(nrays):
                    
                            Lk = L_k(r1, p1, f1, r2, p2, f2)
                            j1 = r1 + nrays*(f1 + nfreqs* b1             ) 
                            j2 = r2 + nrays*(f2 + nfreqs*(p2 + nboundary))

                            data[index] = Lk
                            id_i[index] = j1
                            id_j[index] = j2
                            index += 1
                            
                            data[index] = Lk
                            id_i[index] = j2
                            id_j[index] = j1 
                            index += 1
                        
#                             print ('Lk', Lk, p1, r1, p2, r2)
                            
    # For all (general) points
    for p1 in range(npoints):
        i1 = p1 + nboundary
        
        for f1 in range(nfreqs):
            for r1 in range(nrays):
            
                # For all (general) points (close to the considered general point p1)
                for i2 in range(min([len(ind_all[p1]), npoints])):
                    p2 = ind_all[p1][i2]
                    
                    for f2 in range(nfreqs):
                        for r2 in range(nrays):
                            
                            L2 = L2_k(r1, p1, f1, r2, p2, f2)
                            j1 = r1 + nrays*(f1 + nfreqs* i1             ) 
                            j2 = r2 + nrays*(f2 + nfreqs*(p2 + nboundary))
                            
                            data[index] = L2
                            id_i[index] = j1
                            id_j[index] = j2
                            index += 1
                            
#                             print ('L2  ', p1, r1, ';', p2, r2, ' ', L2)
                            
    return (data, (id_i, id_j))
 

# @njit
def get_condition():
    
    # initialize
    condition = np.zeros(matsize)
    # compute

    # For all boundary points
    for i1 in range(nboundary):
        p1 = boundary2point[i1]
        
        for f1 in range(nfreqs):
            for r1 in range(nrays):    
            
#                 if np.dot([1,0,0], Rs[r1]) > 0:
#                 if (i1 == 0):
            
                j1 = r1 + nrays*(f1 + nfreqs*i1)
                
#                 condition[j1] = boundary_condition[i1,f1]
                condition[j1] = 0.0

    # For all (general) points
    for p1 in range(npoints):
        i1 = p1 + nboundary
        
        for f1 in range(nfreqs):
            for r1 in range(nrays):

                j1 = r1 + nrays*(f1 + nfreqs* i1)
                
#                 condition[j1] = eta[p1,f1] / chi[p1,f1]
                condition[j1] = 1.0

    return condition