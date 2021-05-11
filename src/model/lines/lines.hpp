#pragma once

#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include "lineProducingSpecies/lineProducingSpecies.hpp"
#include "model/thermodynamics/thermodynamics.hpp"


struct Lines
{
    Parameters parameters;

    vector <LineProducingSpecies> lineProducingSpecies;

    // Real1 line;         ///< [Hz] line center frequencies orderd
    Vector<Real> line;     ///< [Hz] line center frequencies (NOT ordered!)
    // Size1 line_index;   ///< index of the corresponding frequency in line

    Vector<Size> nrad_cum;

    // Real1 emissivity;   ///< line emissivity (p,l,k)
    // Real1 opacity;      ///< line opacity    (p,l,k)



    Matrix<Real> emissivity;      ///< line emissivity    (p, lid)
    Matrix<Real> opacity;         ///< line opacity       (p, lid)
    Matrix<Real> inverse_width;   ///< inverse line width (p, lid)

    void read  (const Io& io);
    void write (const Io& io) const;

    void iteration_using_LTE (
        const Double2      &abundance,
        const Vector<Real> &temperature);

    void iteration_using_statistical_equilibrium (
        const Double2      &abundance,
        const Vector<Real> &temperature,
        const Real          pop_prec,
        vector<Size> &points_in_grid            );

    void iteration_using_Ng_acceleration (
        const Real pop_prec,
        vector<Size> &points_in_grid            );

    inline Size      index (const Size p, const Size line_index     ) const;
    inline Size line_index (              const Size l, const Size k) const;
    inline Size      index (const Size p, const Size l, const Size k) const;

    inline void set_emissivity_and_opacity ();
    inline void set_inverse_width (const Thermodynamics& thermodynamics);

    void gather_emissivities_and_opacities ();

    inline void set_all_level_pops(vector<VectorXr> new_population);
    inline vector<VectorXr> get_all_level_pops();

    inline void write_populations_of_iteration(const Io& io, const Size it, const Size lvl) const;
    inline void read_populations_of_iteration(const Io& io, const Size it, const Size lvl);
};


#include "lines.tpp"
