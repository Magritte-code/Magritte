#pragma once


#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"


struct Temperature
{
    Parameters parameters;

    Vector<Real> gas;   ///< [K] gas temperature

    void print()
    {
        cout << "data = " << gas.vec.data() << endl;
        // Real* ptr = gas.dat;
        // cout << "test = " << ptr--          << endl;
        cout << "dat  = " << gas.dat        << endl;

        for (Size i = 0; i < parameters.npoints(); i++)
        {
            cout << "gas[i] = " << gas[i] << "    gas.vec[i] = " << gas.vec[i] << endl;
        }
    }


    void read  (const Io& io);
    void write (const Io& io) const;
};
