#include "solver.hpp"


Solver :: Solver (const Size l, const Size w)
    : length (l)
    , centre (l/2)
    , width  (w)
{
    dZ    = (Real*) pc::accelerator::malloc (length*sizeof(Real));
    nr    = (Size*) pc::accelerator::malloc (length*sizeof(Size));
    shift = (Real*) pc::accelerator::malloc (length*sizeof(Real));

    first = (Size*) pc::accelerator::malloc (width *sizeof(Size));
    last  = (Size*) pc::accelerator::malloc (width *sizeof(Size));
}


Solver :: ~Solver ()
{
    pc::accelerator::free (dZ);
    pc::accelerator::free (nr);
    pc::accelerator::free (shift);

    pc::accelerator::free (first);
    pc::accelerator::free (last);
}
