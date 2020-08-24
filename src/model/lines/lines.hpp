#pragma once

#include "io/io.hpp"
#include "tools/types.hpp"
#include "lineProducingSpecies/lineProducingSpecies.hpp"


class Lines
{
    public:
        std::vector <LineProducingSpecies> lineProducingSpecies;

        Double1 line;         ///< [Hz] line center frequencies orderd
        Long1   line_index;   ///< index of the corresponding frequency in line

        Double1 emissivity;   ///< line emissivity (p,l,k)
        Double1 opacity;      ///< line opacity    (p,l,k)

        void read  (const Io &io);
        void write (const Io &io) const;

        accel inline void set_npoints (const Size n);
        accel inline Size get_npoints () const;

        int iteration_using_LTE (
            const Double2 &abundance,
            const Double1 &temperature);

        int iteration_using_statistical_equilibrium (
            const Double2 &abundance,
            const Double1 &temperature,
            const double   pop_prec                 );

        int iteration_using_Ng_acceleration (
            const double   pop_prec         );


        // Inline functions
        inline long index (
            const long p,
            const int  l,
            const int  k) const;

        inline long index (
            const long p,
            const long line_index) const;

        inline void set_emissivity_and_opacity ();


        int gather_emissivities_and_opacities ();


    private:
        Size npoints;
        Size nlines;
        Size nlspecs;

        Long1 nrad_cum;
};


#include "lines.tpp"