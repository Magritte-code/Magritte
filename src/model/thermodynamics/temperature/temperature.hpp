#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"

struct Temperature {
    std::shared_ptr<Parameters> parameters; ///< data structure containing model

    Vector<Real> gas; ///< [K] gas temperature

    Temperature(std::shared_ptr<Parameters> params) : parameters(params){};

    void read(const Io& io);
    void write(const Io& io) const;
};
