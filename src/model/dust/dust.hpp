#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "model/thermodynamics/thermodynamics.hpp"


struct Dust
{
    Parameters parameters;

    Vector<Real> freqs;     ///< [Hz] Frequencies in opacity table (assume ordered AND equidistantly spaced!)
    Vector<Real> kappa;     ///< opacity (f)
    Vector<Real> density;   ///< dust density (p)

    Size N              ;   ///< Number of elements in opacity table minus one
    Real         delta_f;   ///< frequency spacing in opacity table
    Real inverse_delta_f;   ///< inverse of frequency spacing in opacity table

    void read  (const Io& io);
    void write (const Io& io) const;

    inline Real get_opacity   (const Size p, const Real freq) const;
};


#include "dust.tpp"

