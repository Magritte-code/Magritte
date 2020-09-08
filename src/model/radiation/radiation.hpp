#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "frequencies/frequencies.hpp"
#include "scattering/scattering.hpp"


///  Radiation: data structure for the radiation field
//////////////////////////////////////////////////////
class Radiation
{
    public:
        Frequencies frequencies;
        Scattering  scattering;

        Real2 u;         ///< u intensity             (r, index(p,f))
        Real2 v;         ///< v intensity             (r, index(p,f))

        Real2 U;         ///< U scattered intensity   (r, index(p,f))
        Real2 V;         ///< V scattered intensity   (r, index(p,f))

        Real1 J;         ///< (angular) mean intensity (index(p,f))

        Real3 I_bdy;     ///< intensity at the boundary (r,b,f)

        void read  (const Io& io);
        void write (const Io& io) const;


        inline long index (
            const long p,
            const long f  ) const;

        inline long index (
            const long p,
            const long f,
            const long m  ) const;

        inline vReal get_U (
            const long R,
            const long p,
            const long f   ) const;


        inline vReal get_V (
            const long R,
            const long p,
            const long f   ) const;

        inline vReal get_I_bdy (
            const long R,
            const long p,
            const long f       ) const;


#i  f (GRID_SIMD)

      inline double get_U (
          const long R,
          const long p,
          const long f,
          const long lane ) const;

      inline double get_V (
          const long R,
          const long p,
          const long f,
          const long lane ) const;

      inline double get_I_bdy (
          const long R,
          const long p,
          const long f,
          const long lane     ) const;

#e  ndif


      inline double get_u (
          const long r,
          const long p,
          const long f    ) const;

      inline double get_v (
          const long r,
          const long p,
          const long f    ) const;

      inline double get_J (
          const long p,
          const long f    ) const;


      inline void rescale_U_and_V (
          const vReal &freq_scaled,
          const long   R,
          const long   p,
                long  &notch,
                vReal &U_scaled,
                vReal &V_scaled   ) const;

      inline void rescale_I_bdy (
          const vReal &freq_scaled,
          const long   R,
          const long   p,
          const long   b,
                long  &notch,
                vReal &Ibdy_scaled) const;


      int initialize_J ();
      int MPI_reduce_J ();
      int calc_U_and_V ();


    private:
        size_t npoints;               ///< number of points
        size_t nrays;                 ///< number of rays
        size_t nrays_red;             ///< reduced number of rays
        size_t nfreqs;                ///< number of frequencies
        size_t nfreqs_red;            ///< reduced number of frequencies
        size_t nboundary;             ///< number of boundary cells

        bool use_scattering;          ///< number of boundary cells
};


#include "radiation.tpp"
