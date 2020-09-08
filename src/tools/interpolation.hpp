#pragma once


#include "tools/types.hpp"


///  Binary search for the index of a value right above a value in a list
///  @param[in] x: vector of tabulated argument values
///  @param[in] value: value to search for
///  @return index of x table just above value
/////////////////////////////////////////////////////////////////////////
inline Size search (const Real1 &x, const Real value)
{
    Size start = 0;
    Size stop  = x.size()-1;

    if      (value >= x[stop ]){return stop; }
    else if (value <= x[start]){return start;}

    while (stop > start+1)
    {
        const Size middle = (stop + start) / 2;

        if      (value > x[middle]) {start = middle;}
        else if (value < x[middle]) {stop  = middle;}
        else                        {return  middle;}
    }

    return stop;
}


///  Linear interpolation of f(x) in interval [x1, x2]
///    @param[in] x1: function argument 1
///    @param[in] f1: f(x1)
///    @param[in] x2: function argument 2
///    @param[in] f2: f(x2)
///    @param[in] x: value at which the function has to be interpolated
///    @returns interpolated function value f(x)
///////////////////////////////////////////////////////////////////////
inline Real interpolate_linear (
    const Real x1,
    const Real f1,
    const Real x2,
    const Real f2,
    const Real x )
{
    return (f2-f1)/(x2-x1) * (x-x1) + f1;
}