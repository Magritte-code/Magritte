#include "solver.hpp"


Solver :: Solver (const Size l, const Size w)
    : length (l)
    , centre (l/2)
    , width  (w)
{
    dZ     .resize (length);
    nr     .resize (length);
    shift  .resize (length);

    first  .resize (width);
    last   .resize (width);

    I      .resize (width);

    eta_crt.resize (width);
    eta_nxt.resize (width);

    chi_crt.resize (width);
    chi_nxt.resize (width);

    drho   .resize (width);
    dtau   .resize (width);

    tau    .resize (width);
}


void Solver :: initialize (const Size l, const Size w)
{
    dZ     .resize (length);
    nr     .resize (length);
    shift  .resize (length);

    first  .resize (width);
    last   .resize (width);

    I      .resize (width);

    eta_crt.resize (width);
    eta_nxt.resize (width);

    chi_crt.resize (width);
    chi_nxt.resize (width);

    drho   .resize (width);
    dtau   .resize (width);

    tau    .resize (width);
}
