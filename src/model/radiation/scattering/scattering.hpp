#pragma once

#include "tools/types.hpp"

struct Scattering {
    Real1 opacity_scat; ///< scattering opacity (p,f)
    // Precalculate phase function for all frequencies
    Real3 phase; ///< scattering phase function (r1,r2,f)

    void read(const Io& io);
    void write(const Io& io) const;

    Size nfreqs_scat; ///< number of frequencies in scattering data
    Size nfreqs_red;  ///< number of frequencies
};
