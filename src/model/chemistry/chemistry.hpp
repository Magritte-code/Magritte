#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "species/species.hpp"


///  Data structure for Chemistry
/////////////////////////////////
struct Chemistry
{
    Parameters parameters;
    Species    species;

    void read  (const Io& io);
    void write (const Io& io) const;
};
