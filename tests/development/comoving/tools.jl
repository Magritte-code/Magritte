#A partial port of the tools.py file to julia
module Tools

#constants
const kb                         = 1.38064852E-23   #///< [J/K] Boltzmann's constant
const h                          = 6.62607004E-34   #///< [J*s] Planck's constant
const c                          = 2.99792458E+8    #///< [m/s] speed of light
const amu                        = 1.66053904E-27   #///< [kg] proton mass
const tcmb                       = 2.7254800        #///< [K] CMB temperature

function relative_error(a,b)
    """
    Returns the relative error between a and b.
    """
    return 2.0.*(a.-b)./(a.+b)
end


function LTEpop(linedata, temperature)
    """
    Returns the LTE level populations give the temperature.

    Parameters
    ----------
    linedata : Magritte Linedata object
        Magritte linedata object of the of the relevant species.
    temperature : float
        Temperature for which to evaluate the LTE level populations.

    Returns
    -------
    out : array_like
        Array containing the LTE level populations for the given temperature.
    """
    pop = zeros(linedata.nlev)
    # Calculate the LTE populations
    for i in 1:linedata.nlev
        pop[i] = linedata.weight[i] * exp(-linedata.energy[i]/(kb*temperature))
    end
    # Normalize to relative populations
    pop = pop / sum(pop)
    # Done
    return pop
end


function lineEmissivity(linedata, pop)
    """
    Returns the line emissivity for each radiative transition.

    Parameters
    ----------
    linedata : Magritte Linedata object
        Magritte linedata object of the of the relevant species.
    pop : array_like
        Populations of the levels.

    Returns
    -------
    out : array_like
        Array containing the line emissivity function for each radiative transition.
    """
    eta = zeros(linedata.nrad)
    for k in 1:linedata.nrad
        i = linedata.irad[k]
        j = linedata.jrad[k]
        eta[k] = h*linedata.frequency[k]/(4.0*pi) * linedata.A[k]*pop[i]
    end
    # Done
    return eta
end


function lineOpacity(linedata, pop)
    """
    Returns the line opacity for each radiative transition.

    Parameters
    ----------
    linedata : Magritte Linedata object
        Magritte linedata object of the of the relevant species.
    pop : array_like
        Populations of the levels.

    Returns
    -------
    out : array_like
        Array containing the line opacity function for each radiative transition.
    """
    chi = zeros(linedata.nrad)
    for k in 1:linedata.nrad
        i = linedata.irad[k]
        j = linedata.jrad[k]
        chi[k] = h*linedata.frequency[k]/(4.0*pi) * (linedata.Ba[k]*pop[j] - linedata.Bs[k]*pop[i])
    end
    # Done
    return chi
end


function lineSource(linedata, pop)
    """
    Returns the line source function for each radiative transition.

    Parameters
    ----------
    linedata : Magritte Linedata object
        Magritte linedata object of the of the relevant species.
    pop : array_like
        Populations of the levels.

    Returns
    -------
    out : array_like
        Array containing the line source function for each radiative transition.
    """
    S = lineEmissivity(linedata, pop)./lineOpacity(linedata, pop)
    # Done
    return S
end


function planck(temperature, frequency)
    """
    Planck function for thermal radiation.

    Parameters
    ----------
    temperature : float
        Temperature at which to evaluate the intensity.
    frequency : float
        Frequency at which to evaluate the intensity.

    Returns
    -------
    out : float
        Planck function evaluated at the frequency for the given temperature.
    """
    return 2.0*h/c^2 * frequency^3 / expm1(h*frequency/(kb*temperature))
end


function I_CMB(frequency)
    """
    Intensity of the cosmic microwave background.

    Parameters
    ----------
    frequency : float
        Frequency at which to evaluate the intensity.

    Returns
    -------
    out : float
        Intensity of the cosmic microwave background evaluated at the frequency.
    """
    return planck(tcmb, frequency)
end


function dnu(linedata, k, temp, vturb2)
    """
    Spectral line width.

    Parameters
    ----------
    linedata : Magritte Linedata object
        Magritte linedata object of the of the relevant species.
    k : int
        Transition number of the line.
    temp : float
        Local temperature [K].
    vturb2 : float
        Square of the turbulent velocity as fraction of the speed of light.

    Returns
    -------
    out : float
        Line width of the line profile function.
    """
    return linedata.frequency[k] * sqrt(2.0*kb*temp/(amu*c^2)*linedata.inverse_mass + vturb2)
end


function profile(linedata, k, temp, vturb2, nu)
    """
    Gaussian line profile function.

    Parameters
    ----------
    linedata : Magritte Linedata object
        Magritte linedata object of the of the relevant species.
    k : int
        Transition number of the line.
    temp : float
        Local temperature [K].
    vturb2 : float
        Square of the turbulent velocity as fraction of the speed of light.
    nu : float
        Frequency at which to evaluate the line profile function.

    Returns
    -------
    out : float
        Gaussian profile function evaluated at frequency nu.
    """
    x = (nu - linedata.frequency[k]) / dnu(linedata, k, temp, vturb2)
    return exp(-x^2) / (sqrt(np.pi) * dnu(linedata, k, temp, vturb2))
end

end
