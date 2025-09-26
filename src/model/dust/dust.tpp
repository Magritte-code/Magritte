#include "tools/constants.hpp"
#include "tools/types.hpp"

///  Planck function: copied from solver.cpp
///    @param[in] temp : temperature of the corresponding
///    black body
///    @param[in] freq : frequency at which to evaluate the
///    function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
inline Real Dust::planck(Real temp, Real freq) const {
    return TWO_HH_OVER_CC_SQUARED * (freq * freq * freq) / expm1f(HH_OVER_KB * freq / temp);
}

///  Computes the dust opacity by interpolating the precalculated dust opacities
///    @param[in] p : index of the cell
///    @param[in] freq : Comoving frame frequency at which to evaluate the dust opacity
///    @return dust opacity and emissivity at this point, at this frequency
inline void Dust::compute_dust_opacity_emissivity(
    Size p, Real freq, Real& dust_opacity, Real& dust_emissivity) const {
    std::cout << "DEBUG: Computing dust opacity and emissivity at point " << p << " and frequency "
              << freq << std::endl;
    // Return 0 if no dust present
    if (this->n_dust_frequencies == 0) {
        std::cout
            << "DEBUG: No dust present in the model, returning zero dust opacity and emissivity."
            << std::endl;
        dust_opacity    = 0.0;
        dust_emissivity = 0.0;
        return;
    }
    // If frequency is outside the range of the dust opacities, return the boundary value
    if (freq < dust_frequencies[0]) {
        std::cout << "DEBUG: Frequency too low." << std::endl;
        dust_opacity    = dust_opacities(p, 0);
        dust_emissivity = dust_opacity * planck(dust_temperature[p], freq);
        return;
    }
    if (freq > dust_frequencies[dust_frequencies.size() - 1]) {
        std::cout << "DEBUG: Frequency too high." << std::endl;
        dust_opacity    = dust_opacities(p, dust_frequencies.size() - 1);
        dust_emissivity = dust_opacity * planck(dust_temperature[p], freq);
        return;
    }

    // find the right bin using binary search
    Size f_low  = 0;
    Size f_high = dust_frequencies.size() - 1;
    while (f_high - f_low > 1) {
        Size f_mid = (f_high + f_low) / 2;
        if (dust_frequencies[f_mid] > freq) {
            f_high = f_mid;
        } else {
            f_low = f_mid;
        }
    }

    // and do linear interpolation
    Real slope = (dust_opacities(p, f_high) - dust_opacities(p, f_low))
               / (dust_frequencies[f_high] - dust_frequencies[f_low]);
    dust_opacity    = dust_opacities(p, f_low) + slope * (freq - dust_frequencies[f_low]);
    dust_emissivity = dust_opacity * planck(dust_temperature[p], freq);
    return;
}