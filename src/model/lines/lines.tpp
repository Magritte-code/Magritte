#include "paracabs.hpp"
#include "tools/types.hpp"

///  Indexer for point and line indices
///    @param[in] p          : index of the point
///    @param[in] line_index : index of the line
/////////////////////////////////////////////////
inline Size Lines ::index(const Size p, const Size line_index) const {
    return line_index + p * parameters->nlines();
}

///  Indexer for cell, line producing species and transition indices
///    @param[in] l : index of the line producing species
///    @param[in] k : index of the line transition
////////////////////////////////////////////////////////////////////
inline Size Lines ::line_index(const Size l, const Size k) const { return k + nrad_cum[l]; }

///  Indexer for point, line producing species and transition indices
///    @param[in] p : index of the point
///    @param[in] l : index of the line producing species
///    @param[in] k : index of the line transition
/////////////////////////////////////// /////////////////////////////
inline Size Lines ::index(const Size p, const Size l, const Size k) const {
    return index(p, line_index(l, k));
}

///  Setter for line emissivity and opacity
///////////////////////////////////////////
inline void Lines ::set_emissivity_and_opacity() {
    threaded_for(p, parameters->npoints(), {
        for (Size l = 0; l < parameters->nlspecs(); l++) {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++) {
                const Size lid = line_index(l, k);

                emissivity(p, lid) = lineProducingSpecies[l].get_emissivity(p, k);
                opacity(p, lid)    = lineProducingSpecies[l].get_opacity(p, k);
            }
        }
    })

    Real min_opacity = opacity.min();
    if (min_opacity < parameters->line_opacity_warning_threshold) {
        cout << "Warning: Significant negative line opacity detected. Magritte does not handle "
                "masers and will ignore this. Minimum line opacity in the model: "
             << min_opacity << " [m^-1 Hz]" << endl;
    } else if (min_opacity < 0) {
        cout << "Minor warning: Negative opacity detected. This can be a numerical error, and will "
                "be ignored by Magritte."
             << endl;
    }
    // Now also setting the opacity to the minimum allowed value
    // threaded_for(p, parameters->npoints(), {
    //     for (Size l = 0; l < parameters->nlspecs(); l++) {
    //         for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++) {
    //             const Size lid = line_index(l, k);

    //             if (opacity(p, lid) < parameters->min_line_opacity) {
    //                 opacity(p, lid) = parameters->min_line_opacity;
    //             }
    //         }
    //     }
    // })
}

///  Setter for line widths
///    @param[in] thermodynamics : reference to thermodynamics module
/////////////////////////////////////////////////////////////////////
inline void Lines ::set_inverse_width(const Thermodynamics& thermodynamics) {
    threaded_for(p, parameters->npoints(), {
        for (Size l = 0; l < parameters->nlspecs(); l++) {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++) {
                const Real invr_mass = lineProducingSpecies[l].linedata.inverse_mass;
                const Real frequency = lineProducingSpecies[l].linedata.frequency[k];

                const Size lid = line_index(l, k);

                inverse_width(p, lid) =
                    (Real)1.0 / thermodynamics.profile_width(invr_mass, p, frequency);
            }
        }
    })
}
