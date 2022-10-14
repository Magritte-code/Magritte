#include "paracabs.hpp"
#include "tools/types.hpp"


///  Indexer for point and line indices
///    @param[in] p          : index of the point
///    @param[in] line_index : index of the line
/////////////////////////////////////////////////
inline Size Lines :: index (const Size p, const Size line_index) const
{
    return line_index + p*parameters->nlines();
}


///  Indexer for cell, line producing species and transition indices
///    @param[in] l : index of the line producing species
///    @param[in] k : index of the line transition
////////////////////////////////////////////////////////////////////
inline Size Lines :: line_index (const Size l, const Size k) const
{
    return k + nrad_cum[l];
}


///  Indexer for point, line producing species and transition indices
///    @param[in] p : index of the point
///    @param[in] l : index of the line producing species
///    @param[in] k : index of the line transition
/////////////////////////////////////// /////////////////////////////
inline Size Lines :: index (const Size p, const Size l, const Size k) const
{
    return index (p, line_index (l, k));
}


///  Setter for line emissivity and opacity
///////////////////////////////////////////
inline void Lines :: set_emissivity_and_opacity ()
{
    threaded_for (p, parameters->npoints(),
    {
        for (Size l = 0; l < parameters->nlspecs(); l++)
        {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Size lid = line_index (l, k);

                emissivity (p, lid) = lineProducingSpecies[l].get_emissivity (p, k);
                   opacity (p, lid) = lineProducingSpecies[l].get_opacity    (p, k);
            }
        }
    })
}


///  Setter for line widths
///    @param[in] thermodynamics : reference to thermodynamics module
/////////////////////////////////////////////////////////////////////
inline void Lines :: set_inverse_width (const Thermodynamics& thermodynamics)
{
    threaded_for (p, parameters->npoints(),
    {
        for (Size l = 0; l < parameters->nlspecs(); l++)
        {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Real invr_mass = lineProducingSpecies[l].linedata.inverse_mass;
                const Real frequency = lineProducingSpecies[l].linedata.frequency[k];

                const Size lid = line_index (l, k);

                inverse_width (p, lid) = (Real) 1.0 / thermodynamics.profile_width (invr_mass, p, frequency);
            }
        }
    })
}



inline Size Lines :: convert_to_table_index (double x) const
{
    //convert to int/size, then clamp?
    // const Real clamped_x = std::clamp(x, -parameters->max_distance_opacity_contribution, parameters->max_distance_opacity_contribution);
    // return std::clamp((Size) std::round((x + parameters->max_distance_opacity_contribution) * (parameters->n_tabulated_profile_funs - 1)/(2.0*parameters->max_distance_opacity_contribution))
    //                         , (Size) 0, parameters->n_tabulated_profile_funs - 1);

    const Real clamped_x = std::clamp(x + MAX_DISTANCE_INTERVAL, 0.0, 2.0*MAX_DISTANCE_INTERVAL);
    // std::cout<<"x: "<<x<<" index: "<<((Size) std::round(clamped_x*MULTIPLICATION_FACTOR))<<std::endl;
    // return ((Size) std::round(clamped_x*MULTIPLICATION_FACTOR));
    return (Size) (clamped_x*MULTIPLICATION_FACTOR + 0.5);
}


///   Precomputes gaussians for removing evaluation
inline void Lines :: set_tabulated_gaussians ()
{
    threaded_for (i, N_TABULATED_PROFILE_FUNS,
    {
        const Real x = -MAX_DISTANCE_INTERVAL
        + (Real)i /((Real)(N_TABULATED_PROFILE_FUNS - 1)) * 2.0 * MAX_DISTANCE_INTERVAL;
        // std::cout<<"x: "<<x<<std::endl;
        const Real minxsqaured = -std::pow(x,2);
        tabulated_gaussians[i] = exp(minxsqaured);
    })
}

inline Real Lines :: compute_tabulated_gaussian (double x) const
{
    return tabulated_gaussians[convert_to_table_index(x)];
}
