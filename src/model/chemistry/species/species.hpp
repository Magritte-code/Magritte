#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


///  Data structure for Species
///////////////////////////////
struct Species
{
    std::shared_ptr<Parameters> parameters;   ///< data structure containing model parameters

    String1 symbol;                           ///< Chemical symbol of the species

    Double2 abundance_init;                   ///< abundance before chemical evolution
    Double2 abundance;                        ///< (current) abundance in every cell


    Species (std::shared_ptr<Parameters> params)
    : parameters (params) {};

    void read  (const Io& io);
    void write (const Io& io) const;
};
