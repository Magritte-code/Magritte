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


    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        eta_crt(i).resize (width);
        eta_nxt(i).resize (width);

        chi_crt(i).resize (width);
        chi_nxt(i).resize (width);

        drho   (i).resize (width);
        dtau   (i).resize (width);

        tau    (i).resize (width);
    }
}


//void Solver :: initialize (const Size l, const Size w)
//{
//    dZ     .resize (length);
//    nr     .resize (length);
//    shift  .resize (length);
//
//    first  .resize (width);
//    last   .resize (width);
//
//
//    eta_crt.resize (width);
//    eta_nxt.resize (width);
//
//    chi_crt.resize (width);
//    chi_nxt.resize (width);
//
//    drho   .resize (width);
//    dtau   .resize (width);
//
//    tau    .resize (width);
//}
