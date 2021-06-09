#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Temperature
{
    Parameters parameters;

    Vector<Real> gas;   ///< [K] gas temperature

    void read  (const Io& io);
    void write (const Io& io) const;
};
