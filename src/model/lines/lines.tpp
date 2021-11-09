#include "paracabs.hpp"
#include "tools/types.hpp"


///  Indexer for point and line indices
///    @param[in] p          : index of the point
///    @param[in] line_index : index of the line
/////////////////////////////////////////////////
inline Size Lines :: index (const Size p, const Size line_index) const
{
    return line_index + p*parameters.nlines();
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
    threaded_for (p, parameters.npoints(),
    {
        for (Size l = 0; l < parameters.nlspecs(); l++)
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
    threaded_for (p, parameters.npoints(),
    {
        for (Size l = 0; l < parameters.nlspecs(); l++)
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


///  Calculates the cooling rates
///  NOTE: call this only after computing the radiation field
///////////////////////////////////////////
inline void Lines :: calculate_cooling_rates ()
{
    // Matrix<Real> partial_cooling_rates; //for each point and line species
    // partial_cooling_rates.resize (parameters.npoints(),  parameters.nlspecs());
    threaded_for (p, parameters.npoints(),
    {
        Real cooling_of_point=0;
        // TODO: if certain species may not participate, exclude them here
        for (Size l = 0; l < parameters.nlspecs(); l++)
        {
            Real cooling_of_lspec=0;
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Size lid = line_index (l, k);

                cooling_of_lspec+=emissivity(p, lid)*line[lid]*FOUR_PI;//emissivity*freq(->energy)*4pi(->integrate out direction)
                cooling_of_lspec-=opacity   (p, lid)*line[lid]*FOUR_PI*lineProducingSpecies[l].Jlin[p][k];//also times averaged intensity
            }
            cooling_of_point+=cooling_of_lspec;
        }
        // //TODO: add fancy reduction operation instead of this
        cooling_rates[p]=cooling_of_point;
    })

}
