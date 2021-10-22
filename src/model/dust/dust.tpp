#include <cmath>


///  Getter for the dust opacity
///    @param[in] p: point index
///    @param[in] freq: frequency
/////////////////////////////////
inline Real Dust :: get_opacity (const Size p, const Real freq) const
{
    Real opacity = density[p];

    if      (freq <= freqs[0])
    {
        opacity *= kappa[0];
    }
    else if (freq >= freqs[N])
    {
        opacity *= kappa[N];
    }
    else
    {
        // Interpolate, assuming constant spacing in freqs
        const Real index = (freq - freqs[0]) * inverse_delta_f;
        const Size f     = trunc(index);
        const Real step  = index - f;

        opacity *= kappa[f] + (kappa[f+1] - kappa[f]) * step;
    }

    return opacity;
}
