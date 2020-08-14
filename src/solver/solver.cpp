#include "solver.hpp"


Solver :: Solver (const size_t l, const size_t w)
    : length (l)
    , centre (l/2)
    , width  (w)
{
    dZ    = (double*) pc::accelerator::malloc (length*sizeof(double));
    nr    = (size_t*) pc::accelerator::malloc (length*sizeof(size_t));
    shift = (double*) pc::accelerator::malloc (length*sizeof(double));

    first = (size_t*) pc::accelerator::malloc (width *sizeof(size_t));
    last  = (size_t*) pc::accelerator::malloc (width *sizeof(size_t));
}


Solver :: ~Solver ()
{
    pc::accelerator::free (dZ);
    pc::accelerator::free (nr);
    pc::accelerator::free (shift);

    pc::accelerator::free (first);
    pc::accelerator::free (last);
}
