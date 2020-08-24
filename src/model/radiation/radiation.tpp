// Magritte: Multidimensional Accelerated General-purpose Radiative Transfer
//
// Developed by: Frederik De Ceuster - University College London & KU Leuven
// _________________________________________________________________________


#include "Functions/interpolation.hpp"
#include "Tools/debug.hpp"


inline long Radiation ::
    index            (
        const long p,
        const long f ) const
{
  return f + p * nfreqs_red;
}


inline vReal Radiation ::
    get_U (
        const long R,
        const long p,
        const long f   ) const
{
  return U[R][index (p,f)];
}

inline vReal Radiation ::
    get_V (
        const long R,
        const long p,
        const long f   ) const
{
  return V[R][index (p,f)];
}


inline vReal Radiation ::
    get_I_bdy (
        const long R,
        const long b,
        const long f   ) const
{
  return I_bdy[R][b][f];
}



#if (GRID_SIMD)


inline double Radiation ::
    get_U (
        const long R,
        const long p,
        const long f,
        const long lane) const
{
  return U[R][index (p,f)].getlane (lane);
}


inline double Radiation ::
    get_V (
        const long R,
        const long p,
        const long f,
        const long lane) const
{
  return V[R][index (p,f)].getlane (lane);
}


inline double Radiation ::
    get_I_bdy (
        const long R,
        const long b,
        const long f,
        const long lane) const
{
  return I_bdy[R][b][f].getlane (lane);
}


#endif




inline double Radiation ::
    get_u (
        const long r,
        const long p,
        const long f ) const

#if (GRID_SIMD)

{
  const long indx = newIndex (f);
  const  int lane = laneNr   (f);

  return u[r][index(p,indx)].getlane(lane);
}

#else

{
  return u[r][index(p,f)];
}

#endif




inline double Radiation ::
    get_v (
        const long r,
        const long p,
        const long f ) const

#if (GRID_SIMD)

{
  const long indx = newIndex (f);
  const  int lane = laneNr   (f);

  return v[r][index(p,indx)].getlane(lane);
}

#else

{
  return v[r][index(p,f)];
}

#endif




inline double Radiation ::
    get_J (
        const long p,
        const long f ) const

#if (GRID_SIMD)

{
  const long indx = newIndex (f);
  const  int lane = laneNr   (f);

  return J[index(p,indx)].getlane(lane);
}

#else

{
  return J[index(p,f)];
}

#endif




inline void Radiation ::
    rescale_U_and_V (
        const vReal &freq_scaled,
        const long   R,
        const long   p,
              long  &notch,
              vReal &U_scaled,
              vReal &V_scaled    ) const

#if (GRID_SIMD)

{

  vReal nu1, nu2, U1, U2, V1, V2;

  GRID_FOR_ALL_LANES (lane)
  {

    const double freq = freq_scaled.getlane (lane);

    search_with_notch (frequencies.nu[p], notch, freq);


    if ( (notch == 0) || (notch == nfreqs-1) )
    {
      const long f1    = newIndex (notch);
      const  int lane1 = laneNr   (notch);

      nu1.putlane (frequencies.nu[p][f1].getlane (lane1), lane);
      nu2.putlane (1+frequencies.nu[p][f1].getlane (lane1), lane);
      // the 1+ is to avoid divide by 0 in interpolation

       U1.putlane (get_U (R, p, f1, lane1), lane);
       U2.putlane (get_U (R, p, f1, lane1), lane);

       V1.putlane (get_V (R, p, f1, lane1), lane);
       V2.putlane (get_V (R, p, f1, lane1), lane);
    }

    else
    {

      const long f1    = newIndex (notch);
      const  int lane1 = laneNr   (notch);

      const long f2    = newIndex (notch-1);
      const  int lane2 = laneNr   (notch-1);

      nu1.putlane (frequencies.nu[p][f1].getlane (lane1), lane);
      nu2.putlane (frequencies.nu[p][f2].getlane (lane2), lane);

       U1.putlane (get_U (R, p, f1, lane1), lane);
       U2.putlane (get_U (R, p, f2, lane2), lane);

       V1.putlane (get_V (R, p, f1, lane1), lane);
       V2.putlane (get_V (R, p, f2, lane2), lane);
    }
  }

  U_scaled = interpolate_linear (nu1, U1, nu2, U2, freq_scaled);
  V_scaled = interpolate_linear (nu1, V1, nu2, V2, freq_scaled);

}

#else

{

  search_with_notch (frequencies.nu[p], notch, freq_scaled);

  if ( (notch == 0) || (notch == nfreqs-1) )
  {
    U_scaled = get_U (R, p, notch);
    V_scaled = get_V (R, p, notch);
  }

  else
  {
    const long f1 = notch;
    const long f2 = notch-1;

    const double nu1 = frequencies.nu[p][f1];
    const double nu2 = frequencies.nu[p][f2];

    const double U1 = get_U (R, p, f1);
    const double U2 = get_U (R, p, f2);

    const double V1 = get_V (R, p, f1);
    const double V2 = get_V (R, p, f2);

    U_scaled = interpolate_linear (nu1, U1, nu2, U2, freq_scaled);
    V_scaled = interpolate_linear (nu1, V1, nu2, V2, freq_scaled);
  }

}

#endif




inline void Radiation ::
    rescale_I_bdy (
        const vReal &freq_scaled,
        const long   R,
        const long   p,
        const long   b,
              long  &notch,
              vReal &I_bdy_scaled) const

#if (GRID_SIMD)

{

  vReal nu1, nu2, I_bdy1, I_bdy2;

  //cout << "freq_scaled = " << freq_scaled << endl;

  GRID_FOR_ALL_LANES (lane)
  {
    double freq = freq_scaled.getlane (lane);

    search_with_notch (frequencies.nu[p], notch, freq);

    // cout << "notch["<<lane<<"] = "<< notch << endl;

    if ( (notch == 0) || (notch == nfreqs-1) )
    {
      const long f1    = newIndex (notch);
      const  int lane1 = laneNr   (notch);

      nu1.putlane (  frequencies.nu[p][f1].getlane (lane1), lane);
      nu2.putlane (1+frequencies.nu[p][f1].getlane (lane1), lane);
      // the 1+ is to avoid divide by 0 in interpolation

      I_bdy1.putlane (get_I_bdy (R, b, f1, lane1), lane);
      I_bdy2.putlane (get_I_bdy (R, b, f1, lane1), lane);

      // cout << "nu1["<<lane<<"] = "<< nu1.getlane(lane) << endl;
      // cout << "nu2["<<lane<<"] = "<< nu2.getlane(lane) << endl;

      // cout << "Ibdy1["<<lane<<"] = "<< I_bdy1.getlane(lane) << endl;
      // cout << "Ibdy2["<<lane<<"] = "<< I_bdy2.getlane(lane) << endl;
    }

    else
    {
      const long f1    = newIndex (notch);
      const  int lane1 = laneNr   (notch);

      const long f2    = newIndex (notch-1);
      const  int lane2 = laneNr   (notch-1);

      nu1.putlane (frequencies.nu[p][f1].getlane (lane1), lane);
      nu2.putlane (frequencies.nu[p][f2].getlane (lane2), lane);

      I_bdy1.putlane (get_I_bdy (R, b, f1, lane1), lane);
      I_bdy2.putlane (get_I_bdy (R, b, f2, lane2), lane);

      // cout << "nu1["<<lane1<<"] = "<< nu1.getlane(lane1) << endl;
      // cout << "nu2["<<lane2<<"] = "<< nu2.getlane(lane2) << endl;

      // cout << "Ibdy1["<<lane1<<"] = "<< I_bdy1.getlane(lane1) << endl;
      // cout << "Ibdy2["<<lane2<<"] = "<< I_bdy2.getlane(lane2) << endl;
    }
  }

  I_bdy_scaled = interpolate_linear (nu1, I_bdy1, nu2, I_bdy2, freq_scaled);

  // cout << "I_bdy_scaled = " << I_bdy_scaled << endl;

}

#else

{

  //cout << "search notch" << endl;

  //cout << "p = " << p << endl;
  //cout << "freas len " << frequencies.nu.size() << endl;
  //cout << "freas len " << frequencies.nu[p].size() << endl;

  search_with_notch (frequencies.nu[p], notch, freq_scaled);
  //cout << "notch = " << notch << endl;

  if ( (notch == 0) || (notch == nfreqs-1) )
  {
    I_bdy_scaled = get_I_bdy (R, b, notch);
  }

  else
  {
    const long f1 = notch;
    const long f2 = notch-1;

    const double nu1 = frequencies.nu[p][f1];
    const double nu2 = frequencies.nu[p][f2];

    const double I_bdy1 = get_I_bdy (R, b, f1);
    const double I_bdy2 = get_I_bdy (R, b, f2);

    I_bdy_scaled = interpolate_linear (nu1, I_bdy1, nu2, I_bdy2, freq_scaled);
  }

}

#endif
