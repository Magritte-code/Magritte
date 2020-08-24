#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"


class Temperature
{
    public:
        Double1 gas;   ///< [K] gas temperature

        void read  (const Io &io);
        void write (const Io &io) const;

    private:
        Size npoints;
};
