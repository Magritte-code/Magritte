import numpy as np
from numba import njit:q


@njit
def RBF_Lucy(r):
    """
    Lucy's smoothing kernel (Lucy 1977).
    """
    # Ensure to take the absolute value
    rr = np.abs(r)
    # Compute Lucy's kernel
    if (rr < 1.0):
        return (1.0 + 3.0*rr) * (1 - rr)**3
    else:
        return 0.0

    
@njit
def grad_RBF_Lucy(r):
    # Ensure to take the absolute value
    rr = np.abs(r)
    # Compute the gradient of Lucy's kernel
    if (rr < 1.0):
        return -12.0 * (rr - 1.0)**2 * rr * np.sign(r)
    else:
        return 0.0

@njit
def grad2_RBF_Lucy(r):
    # Ensure to take the absolute value
    rr = np.abs(r)
    # Compute the gradient^2 of Lucy's kernel
    if (rr < 1.0):
        return -12.0 * (rr - 1.0) * (3.0 * rr - 1.0)
    else:
        return 0.0
    

@njit
def kernel(i, j, xs, hs):
    """
    k(x,y) = K(|x-y|^2 / h(x)h(y))
    """
    d = xs[i] - xs[j]
    r = np.dot(d,d) / (hs[i] * hs[j])
    return KK(r)


@njit
def dj_kernel(i, j, xs, hs):
    """
    ∂k(x,y)/∂y
    """
    d = xs[i] - xs[j]
    a = 1.0 / (hs[i] * hs[j])
    r = np.dot(d,d) * a
    return (-2.0*d*a - r*dlog_hs[j]) * grad_KK(r)


@njit
def di_kernel(i, j):
    """
    ∂k(x,y)/∂x
    """
    d = xs[i] - xs[j]
    a = 1.0 / (hs[i] * hs[j])
    r = np.dot(d,d) * a
    return (+2.0*d*a - r*dlog_hs[i]) * grad_KK(r)


@njit
def dij_kernel(i, j):
    """
    ∂∂k(x,y)/∂x∂y
    """
    d = xs[i] - xs[j]
    a = 1.0 / (hs[i] * hs[j])
    r = np.dot(d,d) * a
    term1 = -2.0*a*(1.0 + d*(dlog_hs[j] - dlog_hs[i])) - r*dlog_hs[i]*dlog_hs[j] 
    term2 =  (+2.0*d*a - r*dlog_hs[i])
            *(-2.0*d*a - r*dlog_hs[j])
    return term1*grad_KK(r) + term2*grad2_KK(r)
