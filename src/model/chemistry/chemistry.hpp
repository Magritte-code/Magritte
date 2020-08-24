#pragma once

#include "io/io.hpp"
#include "species/species.hpp"


///  Data structure for Chemistry
/////////////////////////////////
class Chemistry
{
    public:
        Species  species;

        void read  (const Io& io);
        void write (const Io& io) const;
};
