#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


///  Data structure for Species
///////////////////////////////
struct Species
{
    Parameters parameters;

    String1 symbol;

    Double2 abundance_init;   ///< abundance before chemical evolution
    Double2 abundance;        ///< (current) abundance in every cell

    void read  (const Io& io);
    void write (const Io& io) const;
};
