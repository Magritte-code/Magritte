#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"

/// Data structure for dust or continuum radiation in general; Note: these opacity sources will be
/// treated in LTE, using a single temperature for all at once
struct Dust {
    Size n_dust_frequencies; ///< number of frequencies at which dust opacities are defined
    std::shared_ptr<Parameters> parameters; ///< data structure containing model parameters
    Vector<Real> dust_temperature;          ///< dust temperature at each point
    Vector<Real> dust_frequencies;          ///< frequency grid at which dust opacities are defined
    Matrix<Real> dust_opacities;            ///< dust opacities (point, dust frequency index)
    // Matrix<Real> dust_emissivities; ///< = dust opacities * B_nu(T_dust)
    //  Note: we cannot precalculate dust emissivities, as it depends on the exact frequency

    Dust(std::shared_ptr<Parameters> params) : parameters(params){};

    void read(const Io& io);
    void write(const Io& io) const;

    void compute_dust_opacity_emissivity(
        Size p, Real freq, Real& dust_opacity, Real& dust_emissivity) const;

    inline Real planck(Real temp, Real freq) const;
};

#include "dust.tpp"