#pragma once
#include "model/model.hpp"
#include "tools/types.hpp"

// Contains a helper class for allowing interpolation when the variation in some quantity is too
// large.
class InterpHelper {

  public:
    Real max_source_diff;
    std::vector<Real> interpolation_criterion; // Note: contains the density normalized line opacity
    ///  Constructor
    ///    @param[in] dshift_max : maximum shift between two points
    /////////////////////////////////////////////////
    InterpHelper(Model& model) : max_source_diff(model.parameters->max_interpolation_diff) {
        // current interpolation criterion based on the line opacity (normalized by the total
        // population)
        interpolation_criterion.resize(model.lines.emissivity.vec.size());
        threaded_for(p, model.parameters->npoints(), {
            for (Size l = 0; l < model.parameters->nlspecs(); l++) {
                for (Size k = 0; k < model.lines.lineProducingSpecies[l].linedata.nrad; k++) {
                    const Size lid = model.lines.line_index(l, k);

                    interpolation_criterion[model.lines.index(p, lid)] =
                        model.lines.opacity(p, lid)
                        / model.lines.lineProducingSpecies[l].population_tot[p];
                }
            }
        })
    }

    /// DO NOT USE THIS CONSTRUCTOR
    // TODO: in solver.hpp/.tpp, replace the object with an std::optional<...> version.
    InterpHelper() : max_source_diff(-1){};

    // NOTE: interpolation is done in log space, as linear space would take far too many points
    inline Size get_n_interp(const Model& model, const Size curr_idx, const Size next_idx) const;
    inline Size get_n_interp_for_line(
        const Model& model, const Size l, const Size curr_idx, const Size next_idx) const;

    // Interpolation functions themselves
    inline Real interpolate_linear(const Real f_start, const Real f_end, const Real factor) const;
    inline Real interpolate_log(const Real f_start, const Real f_end, const Real factor) const;
};

#include "solver_interp_helper.tpp"