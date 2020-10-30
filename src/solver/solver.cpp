#include "solver.hpp"


Solver :: Solver (const Size l, const Size w, const Size n_o_d)
    : length     (l)
    , centre     (l/2)
    , width      (w)
    , n_off_diag (n_o_d)
{

    cout << "centre = " << centre << endl;
    cout << "length = " << length << endl;
    cout << "width  = " << width  << endl;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        dZ_          (i).resize (length);
        nr_          (i).resize (length);
        shift_       (i).resize (length);

        eta_c_       (i).resize (width);
        eta_n_       (i).resize (width);

        chi_c_       (i).resize (width);
        chi_n_       (i).resize (width);

        inverse_chi_ (i).resize (length);

        tau_         (i).resize (width);

        Su_          (i).resize (length);
        Sv_          (i).resize (length);

        A_           (i).resize (length);
        C_           (i).resize (length);
        inverse_A_   (i).resize (length);
        inverse_C_   (i).resize (length);

        FF_          (i).resize (length);
        FI_          (i).resize (length);
        GG_          (i).resize (length);
        GI_          (i).resize (length);
        GP_          (i).resize (length);

        L_diag_      (i).resize (length);

        L_upper_     (i).resize (n_off_diag, length);
        L_lower_     (i).resize (n_off_diag, length);
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
