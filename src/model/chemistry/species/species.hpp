#pragma once

#include "io/io.hpp"
#include "tools/types.hpp"


class Species
{
    public:
        Double2 abundance_init;   ///< abundance before chemical evolution
        Double2 abundance;        ///< (current) abundance in every cell

        void read  (const Io& io);
        void write (const Io& io) const;

    private:
        Size npoints;   ///< number of points
        Size nspecs;    ///< number of chemical species
};
