#include "paracabs.hpp"
#include "tools/types.hpp"


///  Indexer for cell, line producing species and transition indices
///    @param[in] p : index of the cell
///    @param[in] l : index of the line producing species
///    @param[in] k : index of the line transition
////////////////////////////////////////////////////////////////////
inline Size Lines :: index (const Size p, const Size l, const Size k) const
{
    return k + nrad_cum[l] + p*parameters.nlines();
}


///  Indexer for cell and line indices
///    @param[in] p          : index of the cell
///    @param[in] line_index : index of the line
////////////////////////////////////////////////
inline Size Lines :: index (const Size p, const Size line_index) const
{
    return line_index + p*parameters.nlines();
}


///  Setter for line emissivity and opacity
///    @param[in] p : index of the cell
///    @param[in] l : index of the line producing species
/////////////////////////////////////////////////////////
inline void Lines :: set_emissivity_and_opacity ()
{
    threaded_for (p, parameters.npoints(),
    {
        for (Size l = 0; l < parameters.nlspecs(); l++)
        {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Size ind = index (p, l, k);

                emissivity[ind] = lineProducingSpecies[l].get_emissivity (p, k);
                   opacity[ind] = lineProducingSpecies[l].get_opacity    (p, k);
            }
        }
    })
}
