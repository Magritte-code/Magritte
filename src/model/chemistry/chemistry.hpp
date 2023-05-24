#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "species/species.hpp"

///  Data structure for Chemistry
/////////////////////////////////
struct Chemistry {
    std::shared_ptr<Parameters> parameters; ///< data structure containing model parameters

    Species species; ///< data structure containing chemical species

    Chemistry(std::shared_ptr<Parameters> params) : parameters(params), species(params){};

    void read(const Io& io);
    void write(const Io& io) const;
};
