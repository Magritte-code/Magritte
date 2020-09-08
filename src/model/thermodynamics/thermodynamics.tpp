#include "tools/constants.hpp"


///  profile: line profile function
///    @param[in] inverse_width: inverse profile width
///    @param[in] freq_diff: frequency at which we want evaluate the profile
///    @return profile function evaluated at frequency freq
////////////////////////////////////////////////////////////////////////////
inline Real Thermodynamics :: profile (
    const Real width,
    const Real freq_diff ) const
{
    const Real inverse_width = Real (1.0) / width;
    const Real sqrtExponent  = inverse_width * freq_diff;
    const Real exponent      = sqrtExponent * sqrtExponent;

    return inverse_width * INVERSE_SQRT_PI * expf (-exponent);
}


inline Real Thermodynamics :: profile (
    const Real inverse_mass,
    const Size p,
    const Real freq_line,
    const Real freq ) const
{
    return profile (profile_width (inverse_mass, p, freq_line), freq - freq_line);
}


///  profile_width: line profile width due to thermal and turbulent Doppler shifts
///    @param[in] temperature_gas: temperature of the gas at this cell
///    @param[in] freq_line: frequency of the line under consideration
///    @return width of the correpsonding line profile
//////////////////////////////////////////////////////////////////////////////////
inline Real Thermodynamics :: profile_width (
    const Real inverse_mass,
    const Size p,
    const Real freq_line ) const
{
    return freq_line * profile_width (inverse_mass, p);
}


inline Real Thermodynamics :: profile_width (
    const Real inverse_mass,
    const Size p ) const
{
    return sqrt (TWO_KB_OVER_AMU_CC_SQUARED * inverse_mass * temperature.gas[p]
                 + turbulence.vturb2[p]);
}
