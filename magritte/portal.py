# (c) Boy Lankhaar, 2022
#     lankhaar@chalmers.se / lankhaar@strw.leidenuniv.nl
#
#     This set of routines compute the angular momentum coupling
#     coefficients that are relevant to the alginment specific
#     radiative interactions. See Eqs. (49)-(54) of Landi Degl'
#     Innocenti 1984 (Sol. Ph. 91 1L) for the related derivation
#     and expressions.
#
#     All these subroutines make use of the earlier initiated set
#     of Wigner matrices. These are made with the wigxjpf package.
#

import numpy     as np
import pywigxjpf as wig

from numba   import njit
from astropy import constants


wig.wig_table_init(2*100, 9)
wig.wig_temp_init (2*100)


@njit
def align_book(j_lev):
    """
    Set up bookkeeping array for the alignment states.

    Parameters
    ----------
    j_lev : array_like
        Array containing 2x angular momentum of the energy levels.

    Returns
    -------
    a_lev : array_like
        Bookkeeping array.
    sa : int
        Total number of alignment states.
    """

    nlev = j_lev.size

    a_lev = np.zeros((nlev,2), dtype=np.int64)

    sa = 0
    for s in range(nlev):
        km = 1
        if (j_lev[s] >= 2):
            km = 2
        for k in range(km):
            a_lev[s,k] = sa
            sa += 1

    return a_lev, sa


@njit
def get_fk(k1, k2, K):
    """
    input:
        k1, k2 : irreducible tensor rank of states 1 and 2
        K      : irreducible tensor rank of radiation field
    output:
        sqrt(3[k1][k2][K])
    """
    nk1 = 2 * k1 + 1
    nk2 = 2 * k2 + 1
    nK  = 2 * K  + 1
    return np.sqrt(3 * nk1 * nk2 * nK)


@njit
def get_fas(j1,j2):
    """
    input:
        j1, j2 : 2*angular momentum of states 1 and 2
    output:
        fas    : phase (-1)^(1 + j1 + j2)
    """
    fas = 1.0
    if (j1+j2 % 4) == 0:
        fas = -1.0
    return fas


def get_rp(j1, j2, k1, k2, K, Bt):
    """
    input:
        j1,j2      : 2*angular momentum of states 1 and 2
        k1,k2      : irreducible tensor rank of states 1 and 2
        K          : irreducible tensor rank of radiation field
        Bt         : einstein coefficient j1 -> j2 * [j1]
    output:
        coupling factor for abs/stim emission to state 1
    """
    fk = get_fk(k1, k2, K)

    w9 = wig.wig9jj(np.int64(2  ), np.int64(  j1), np.int64(  j2),
                    np.int64(2  ), np.int64(  j1), np.int64(  j2),
                    np.int64(2*K), np.int64(2*k1), np.int64(2*k2) )

    w3 = wig.wig3jj(np.int64(2*k1), np.int64(2*k2), np.int64(2*K),
                    np.int64(0),    np.int64(0),    np.int64(0)   )

    return Bt * fk * w9 * w3


def get_rm(j1, j2, k1, k2, K, Bt):
    """
    input:
        j1, j2 : 2*angular momentum of states 1 and 2
        k1, k2 : irreducible tensor rank of states 1 and 2
        K      : irreducible tensor rank of radiation field
        Bt     : einstein coefficient j2 -> j1 * [j2]
    output:
        coupling factor for abs/stim emission from state 1
    """
    fk  = get_fk (k1, k2, K)
    fas = get_fas(j1, j2)

    w6_1 = wig.wig6jj(np.int64(2   ), np.int64(2   ), np.int64(2*K),
                      np.int64(  j1), np.int64(  j1), np.int64( j2) )
    w6_2 = wig.wig6jj(np.int64(2*k1), np.int64(2*k2), np.int64(2*K),
                      np.int64(  j1), np.int64(  j1), np.int64( j1) )
    w3   = wig.wig3jj(np.int64(2*k1), np.int64(2*k2), np.int64(2*K),
                      np.int64(0   ), np.int64(0   ), np.int64(0  ) )

    return fk * fas * w6_1 * w6_2 * w3 * Bt


def get_tp(j1, j2, k1, A):
    """
    input:
        j1, j2 : 2*angular momentum of states 1 and 2
        k1     : irreducible tensor rank of states 1 and 2 (\delta_{k1,k2})
        A      : Einstein A-coefficient j2 -> j1
    output:
        coupling factor for spontaneous emission to state 1
    """

    nj2 =  j2 + 1
    At  = nj2 * A

    w6  = wig.wig6jj(np.int64(j1), np.int64(j1), np.int64(2*k1),
                     np.int64(j2), np.int64(j2), np.int64(2   ) )

    fas = get_fas(j1, j2)

    return At * fas * w6


@njit
def A_to_B(A, nu):
    """
    input:
        A  : Einstein A coefficient
        nu : frequency of the transition (Hz)
    ouput:
        Einstein B coefficient
    """
    return A / (1.4744994647625417e-50 * nu * nu * nu)


def setup_rt(t_ju, t_jd, t_A, nu):
    """
    input:
        t_ju/d : 2*angular momentum of up and down states of the transitions
        t_A    : Einstein A coefficients of the transitions
        nu     : frequency of the transitions (Hz)
    output:
        rp     : populating first order transition coupling coefficients
        rm     : depopulating first order transition coupling coefficients
        tp     : populating spontaneous emission transition coupling coefficients
    """
    t_B  = A_to_B(t_A, nu)
    nju  = t_ju + 1.0
    t_Bt = t_B * nju

    rp = np.zeros((2, t_ju.size, 2, 2, 2))
    rm = np.zeros((2, t_ju.size, 2, 2, 2))
    tp = np.zeros((   t_ju.size, 2))

    for t in range(t_ju.size):
        for ik1 in range(2):
            tp[t,ik1] = get_tp(t_jd[t],t_ju[t],2*ik1,t_A[t])                                   # spont em. to   // r4
            for ik2 in range(2):
                for iK in range(2):
                    rp[0,t,ik1,ik2,iK] = get_rp(t_ju[t], t_jd[t], 2*ik1, 2*ik2, 2*iK, t_Bt[t]) # + absorption   // r3
                    rm[0,t,ik1,ik2,iK] = get_rm(t_ju[t], t_jd[t], 2*ik1, 2*ik2, 2*iK, t_Bt[t]) # - stim. em.    // r2
                    rp[1,t,ik1,ik2,iK] = get_rp(t_jd[t], t_ju[t], 2*ik1, 2*ik2, 2*iK, t_Bt[t]) # + stim em.     // r5
                    rm[1,t,ik1,ik2,iK] = get_rm(t_jd[t], t_ju[t], 2*ik1, 2*ik2, 2*iK, t_Bt[t]) # - absorption   // r6

    return rp, rm, tp


def w_fac(j1, j2):
    """
    function   w_fac   -- function to compute the w_{j1,j2}^{(2)}
    radiative transition coupling coefficient for the anisotropic part

    input:
        j1, j2 : 2*angular momenta of the transition states
    output:
        w_{j1,j2}^{(2)} coupling factor
    """
    w6 = wig.wig6jj(np.int64(2 ), np.int64(2 ), np.int64(2*2),
                    np.int64(j1), np.int64(j1), np.int64(j2 ) )

    fac = np.sqrt(3.0*(j1+1.0))

    if (((j1+j2) % 4)==0):
        fac = -fac

    return fac * w6


def w_fac_list(j1, j2):
    """
    function   w_fac_list  -- function to call the w_fac
    for a list of transitions, with j1 > j2.
    input:
        j1, j2 : list of transition angular momenta j1>j2
    output:
        w1     : w_{j1,j2}^{(2)} coupling factor stimulated emission
        w2     : w_{j1,j2}^{(2)} coupling factor stimulated absorption
    """
    w1 = np.zeros(j1.size)
    w2 = np.zeros(j1.size)

    for i in range(j1.size):
        w1[i] = w_fac(j1[i],j2[i])
        w2[i] = w_fac(j2[i],j1[i])

    return w1, w2


def rad_fac_meter(A, nu):
    """
    function   rad_fac_meter   -- function to compute the radiation factor
    k_\lambda for all transitions. In SI units.

    input:
        A      : Einstein A-coefficent * gu
        nu     : Frequency in Hz
    output:
        k0     : radiation factor k_{lambda} = lam^3 * gu * Aul / 8 * pi
    """
    lam_m  = 2.998e8 / nu
    k0     = np.power(lam_m, 3) * A / (8.0*np.pi)
    return k0


def opacity_fac_pol_meter(A, nu, jl, ju, Nl, Nu, sig_l, sig_u):
    """
    function   opacity_fac_pol_meter   -- function to compute the opacity
    (radiation factor * population) for both the isotropic and anisotropic part.
    For now, we keep the absorption and stimulated emission contribution separated.

    input:
        A      : Einstein A-coefficent * gu                [t]
        nu     : Frequency in Hz                           [t]
        jl/u   : list of j-upper and lower quantum numbers [t]
        Nl,Nu  : lower and upper populations / gl,gu       [t,ncell]
        sigl/u : lower and upper alignment                 [t,ncell]
    output:
        k0_l   : isotropic contribution to absorption      [t,ncell]
        k2_l   : anisotropic contribution to absorption    [t,ncell]
        k0_u   : isotropic contribution to stim em         [t,ncell]
        k2_u   : anisotropic contribution to stim em       [t,ncell]
    """
    k0 = rad_fac_meter(A, nu)

    k0_l = k0 * Nl
    k0_u = k0 * Nu
    k2_l = k0 * Nl * sig_l
    k2_u = k0 * Nu * sig_u

    wu, wl = w_fac_list(ju, jl)

    k2_l = wl * k2_l
    k2_u = wu * k2_u

    return k0_l, k0_u, k2_l, k2_u


def compute_anisotropy_and_alignment(model, magnetic):
    """
    Computes the anisotropy of the radiation field w.r.t. the given magnetic field.

    Parameters
    ----------
    model : object
        Magritte model object.
    magnetic : array_like
        Array containing the (3D) magnetic field for each point in the model.

    Returns
    -------
    Model in which the anisotropic radiation field and the alignment are computed.
    """

    bb = np.linalg.norm(magnetic, axis=1)
    magnetic = (magnetic.T / bb).T
    bx = magnetic[:,0]
    by = magnetic[:,1]
    bz = magnetic[:,2]

    model.b.set(magnetic)

    model.compute_spectral_discretisation ()
    model.compute_inverse_line_widths     ()
    model.lines.set_emissivity_and_opacity()

    model.compute_radiation_field_feautrier_order_2_anis()

    for l in range(model.parameters.nlspecs()):

        J0       = np.array(model.lines.lineProducingSpecies[l].J)
        J2_0     = np.array(model.lines.lineProducingSpecies[l].J2_0)
        J2_1_Re  = np.array(model.lines.lineProducingSpecies[l].J2_1_Re)
        J2_1_Im  = np.array(model.lines.lineProducingSpecies[l].J2_1_Im)
        J2_2_Re  = np.array(model.lines.lineProducingSpecies[l].J2_2_Re)
        J2_2_Im  = np.array(model.lines.lineProducingSpecies[l].J2_2_Im)

        J2_1 = J2_1_Re + J2_1_Im * 1j
        J2_2 = J2_2_Re + J2_2_Im * 1j

        C_20 = 0.5 * (3.0*bz**2 - 1.0)
        C_21 = -np.sqrt(3.0/2.0) * (bx*bz + 1j*by*bz)
        C_22 = 0.5 * np.sqrt(3.0/2.0) * (bx**2 - by**2 + 2.0j*bx*by)

        J2 = (np.conj(C_20) * J2_0.T + 2.0 * np.real(np.conj(C_21) * J2_1.T + np.conj(C_22) * J2_2.T)).T

        A     = np.array(model.lines.lineProducingSpecies[l].linedata.A,         dtype=np.float64)
        nu    = np.array(model.lines.lineProducingSpecies[l].linedata.frequency, dtype=np.float64)
        j_lev = np.array(model.lines.lineProducingSpecies[l].linedata.weight,    dtype=np.int64  ) - 1
        lau   = np.array(model.lines.lineProducingSpecies[l].linedata.irad,      dtype=np.int64  )
        lal   = np.array(model.lines.lineProducingSpecies[l].linedata.jrad,      dtype=np.int64  )

        # Setup RT matrix elements
        rp, rm, tp = setup_rt(j_lev[lau], j_lev[lal], A, nu)

        # Get bookkeeping variables
        a_lev, nalign = align_book(j_lev)

        model.lines.lineProducingSpecies[l].nalign = nalign
        model.lines.lineProducingSpecies[l].j_lev  = np.array(j_lev)
        model.lines.lineProducingSpecies[l].a_lev  = np.array(a_lev)
        model.lines.lineProducingSpecies[l].J0     = np.array(J0   )
        model.lines.lineProducingSpecies[l].J2     = np.array(J2   )
        model.lines.lineProducingSpecies[l].rp     = np.array(rp   )
        model.lines.lineProducingSpecies[l].rm     = np.array(rm   )
        model.lines.lineProducingSpecies[l].tp     = np.array(tp   )

    model.PORTAL_solve_statistical_equilibrium()

    return model


def compute_image(model, vpix, l=0, line=0, ray_nr=0):
    """
    Imager for the polarized radiation field.

    Parameters
    ----------
    model : object
        Magritte model object.
    vpix : float
        Width (in velocity [m/s]) of the velocity bins.
    l : int
        Index of the line producing species.
    line : int
        Index of the line.
    ray_nr : int
        Index of the ray along which the image is made.

    Returns
    -------
    Model in which the polarized images is computed.
    """

    npoints =          model.parameters.npoints()
    nlev    =          model.lines.lineProducingSpecies[l].linedata.nlev
    A       = np.array(model.lines.lineProducingSpecies[l].linedata.A,         dtype=np.float64)
    rho     = np.array(model.lines.lineProducingSpecies[l].rho,                dtype=np.float64)
    weight  = np.array(model.lines.lineProducingSpecies[l].linedata.weight,    dtype=np.float64)
    nu      = np.array(model.lines.lineProducingSpecies[l].linedata.frequency, dtype=np.float64)
    pops    = np.array(model.lines.lineProducingSpecies[l].population,         dtype=np.float64).reshape((npoints,nlev))
    j_lev   = np.array(model.lines.lineProducingSpecies[l].linedata.weight,    dtype=np.int64  ) - 1
    lau     = np.array(model.lines.lineProducingSpecies[l].linedata.irad,      dtype=np.int64  )
    lal     = np.array(model.lines.lineProducingSpecies[l].linedata.jrad,      dtype=np.int64  )

    # Get bookkeeping variables
    a_lev, nalign = align_book(j_lev)

    N  = pops / weight
    Nl = N[:, lal]
    Nu = N[:, lau]

    sigma = rho[:,a_lev[:,1]] / rho[:,a_lev[:,0]]
    sig_l = sigma[:,lal]
    sig_u = sigma[:,lau]

    jl = j_lev[lal]
    ju = j_lev[lau]

    k0_l, k0_u, k2_l, k2_u = opacity_fac_pol_meter(A * weight[lau], nu, jl, ju, Nl, Nu, sig_l, sig_u)

    model.k_abs_0 = k0_l[:,line]
    model.k_stm_0 = k0_u[:,line]
    model.k_abs_2 = k2_l[:,line]
    model.k_stm_2 = k2_u[:,line]

    nfreqs = model.parameters.nfreqs()
    fcen   = model.lines.lineProducingSpecies[0].linedata.frequency[line]
    dd     = vpix * (nfreqs-1)/2 / constants.c.si.value
    fmin   = fcen - fcen*dd
    fmax   = fcen + fcen*dd

    model.compute_spectral_discretisation (fmin, fmax)
    model.compute_inverse_line_widths     ()
    model.PORTAL_image                    (ray_nr, line)

    return model
