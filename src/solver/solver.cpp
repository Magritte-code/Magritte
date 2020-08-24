#include "solver.hpp"


Solver :: Solver (const Size l, const Size w)
    : length (l)
    , centre (l/2)
    , width  (w)
{
    dZ   .resize (length);
    nr   .resize (length);
    shift.resize (length);

    first.resize (width);
    last .resize (width);
}
