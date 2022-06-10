#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "temperature/temperature.hpp"
#include "turbulence/turbulence.hpp"


struct Thermodynamics
{
    std::shared_ptr<Parameters>  parameters;   ///< data structure containing model

    Temperature temperature;                   ///< data structure containing temperature
    Turbulence  turbulence;                    ///< data structure containing turbulence


    Thermodynamics (std::shared_ptr<Parameters> params)
    : temperature (params)
    , turbulence  (params) {};

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
