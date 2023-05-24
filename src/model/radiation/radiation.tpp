#include "tools/interpolation.hpp"

inline Size Radiation ::index(const Size p, const Size f) const { return f + p * parameters->nfreqs(); }

// inline void Radiation :: rescale_U_and_V (
//     const Real &freq_scaled,
//     const Size  R,
//     const Size  p,
//           Size &notch,
//           Real &U_scaled,
//           Real &V_scaled ) const
//{
//     search_with_notch (frequencies.nu[p], notch, freq_scaled);
//
//     if ( (notch == 0) || (notch == parameters->nfreqs()-1) )
//     {
//         U_scaled = get_U (R, p, notch);
//         V_scaled = get_V (R, p, notch);
//     }
//
//     else
//     {
//         const Size f1 = notch;
//         const Size f2 = notch-1;
//
//         const Real nu1 = frequencies.nu[p][f1];
//         const Real nu2 = frequencies.nu[p][f2];
//
//         const Real U1 = get_U (R, p, f1);
//         const Real U2 = get_U (R, p, f2);
//
//         const Real V1 = get_V (R, p, f1);
//         const Real V2 = get_V (R, p, f2);
//
//         U_scaled = interpolate_linear (nu1, U1, nu2, U2, freq_scaled);
//         V_scaled = interpolate_linear (nu1, V1, nu2, V2, freq_scaled);
//     }
// }
//
//
// inline void Radiation :: rescale_I_bdy (
//     const Real &freq_scaled,
//     const Size  R,
//     const Size  p,
//     const Size  b,
//           Size &notch,
//           Real &I_bdy_scaled) const
//{
//     //cout << "search notch" << endl;
//
//     //cout << "p = " << p << endl;
//     //cout << "freas len " << frequencies.nu.size() << endl;
//     //cout << "freas len " << frequencies.nu[p].size() << endl;
//
//     search_with_notch (frequencies.nu[p], notch, freq_scaled);
//     //cout << "notch = " << notch << endl;
//
//     if ( (notch == 0) || (notch == parameters->nfreqs()-1) )
//     {
//         I_bdy_scaled = get_I_bdy (R, b, notch);
//     }
//
//     else
//     {
//         const Size f1 = notch;
//         const Size f2 = notch-1;
//
//         const Real nu1 = frequencies.nu[p][f1];
//         const Real nu2 = frequencies.nu[p][f2];
//
//         const Real I_bdy1 = get_I_bdy (R, b, f1);
//         const Real I_bdy2 = get_I_bdy (R, b, f2);
//
//         I_bdy_scaled = interpolate_linear (nu1, I_bdy1, nu2, I_bdy2,
//         freq_scaled);
//     }
// }
