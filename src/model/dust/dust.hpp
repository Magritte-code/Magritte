#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"

/// Data structure for dust
struct Dust {
    std::shared_ptr<Parameters> parameters; ///< data structure containing model parameters
    Matrix<Real> dust_opacities;            ///< dust opacities (point, line index)
    Matrix<Real> dust_emissivities;         ///< dust emissivities (point, line index)
    // Note: dust opacities and emissivities are to be precalculated

    Dust(std::shared_ptr<Parameters> params) : parameters(params){};

    void read(const Io& io);
    void write(const Io& io) const;
}