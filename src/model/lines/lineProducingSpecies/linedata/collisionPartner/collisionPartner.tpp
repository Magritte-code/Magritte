#include "tools/interpolation.hpp"


///  interpolate_collision_coefficients:
///    @param[in] temperature_gas: local gas temperature
////////////////////////////////////////////////////////
inline void CollisionPartner :: interpolate_collision_coefficients (const Real temperature_gas)
{
    const Size t = search (tmp, temperature_gas);

    if (t == 0)
    {
        Ce_intpld = Ce[0];
        Cd_intpld = Cd[0];
    }
    else if (t == ntmp)
    {
        Ce_intpld = Ce[ntmp-1];
        Cd_intpld = Cd[ntmp-1];
    }
    else
    {
        const Real step = (temperature_gas - tmp[t-1]) / (tmp[t] - tmp[t-1]);

        for (Size k = 0; k < ncol; k++)
        {
            Ce_intpld[k] = Ce[t-1][k] + (Ce[t][k] - Ce[t-1][k]) * step;
            Cd_intpld[k] = Cd[t-1][k] + (Cd[t][k] - Cd[t-1][k]) * step;
        }
    }
}


///  Makes correction for ortho and para in H2 abundance
///    @param[in] temperature_gas : gas temperature
///    @param[in/out] abundance       : H2 abundance
////////////////////////////////////////////////////////
inline void CollisionPartner :: adjust_abundance_for_ortho_or_para (
    const Real temperature_gas, Real& abundance) const
{
    if (orth_or_para_H2 != "n")
    {
        const Real frac_H2_para = 1.0 / (1.0 + 9.0*expf (-170.5/temperature_gas));

        if (orth_or_para_H2 == "o")
        {
            abundance *= (1.0 - frac_H2_para);
        }

        if (orth_or_para_H2 == "p")
        {
            abundance *= frac_H2_para;
        }
    }
}
