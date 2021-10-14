#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "tools/interpolation.hpp"
#include "model/thermodynamics/thermodynamics.hpp"


struct Dust
{
    Parameters parameters;

    Vector<Real> freqs;     ///< [Hz] Frequencies in opacity table (assume ordered!)
    Vector<Real> kappa;     ///< opacity (f)
    Vector<Real> density;   ///< dust density (p)

    void read  (const Io& io);
    void write (const Io& io) const;

    inline Real get_opacity   (const Size p, const Real freq) const;
};


#include "dust.tpp"

