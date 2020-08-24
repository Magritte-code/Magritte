#pragma once


#include "io/io.hpp"
#include "tools/types.hpp"
#include "model/thermodynamics/temperature/temperature.hpp"


class Frequencies
{
    public:
        vReal2 nu;                        ///< [Hz] frequencies (ordered in f) (p,f)

        Bool1 appears_in_line_integral;   ///< True if the frequency appears in line integral
        Long1 corresponding_l_for_spec;   ///< number of line species corresponding to frequency
        Long1 corresponding_k_for_tran;   ///< number of transition corresponding to frequency
        Long1 corresponding_z_for_line;   ///< number of line number corresponding to frequency


        inline double get_nu (
            const long p,
            const long f ) const
        {
          return get (nu[p], f);
        }

        void read  (const Io& io);
        void write (const Io& io) const;

        accel inline void set_npoints (const Size n);
        accel inline Size get_npoints () const;

        accel inline void set_nlines (const Size n);
        accel inline Size get_nlines () const;

        accel inline void set_nquads (const Size n);
        accel inline Size get_nquads () const;

    private:
        Size npoints;      ///< number of cells
        Size nlines;       ///< number of lines
        Size nquads;       ///< number frequency quadrature points
//        Size nbins = 0;    ///< number of extra bins per line
//        Size ncont = 0;    ///< number of background bins
        Size nfreqs;       ///< number of frequencies
        Size nfreqs_red;   ///< nfreqs divided by n_simd_lanes
};
