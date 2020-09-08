#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "temperature/temperature.hpp"
#include "turbulence/turbulence.hpp"


struct Thermodynamics
{
    Parameters  parameters;
    Temperature temperature;
    Turbulence  turbulence;

    void read  (const Io& io);
    void write (const Io& io) const;

    inline Real profile (
        const Real width,
        const Real freq_diff ) const;

    inline Real profile (
        const Real inverse_mass,
        const Size p,
        const Real freq_line,
        const Real freq ) const;

    inline Real profile_width (
        const Real inverse_mass,
        const Size p,
        const Real freq_line ) const;

    inline Real profile_width (
        const Real inverse_mass,
        const Size p ) const;
};


#include "thermodynamics.tpp"
