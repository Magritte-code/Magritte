#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


class Turbulence
{
    public:
        Real1 vturb2;   ///< [.] microturbulence over c all squared

        void read  (const Io &io);
        void write (const Io &io) const;

    private:
        Size npoints;
};
