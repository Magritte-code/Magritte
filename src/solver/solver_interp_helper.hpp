#pragma once
#include "model/model.hpp"
#include "tools/types.hpp"

//Contains a helper class for allowing interpolation when the variation in some quantity is too large.
class InterpHelper {

public:
    Real max_source_diff;
    // Model model;
    std::vector<Real> line_sources;
    ///  Constructor
    ///    @param[in] dshift_max : maximum shift between two points
    /////////////////////////////////////////////////
    InterpHelper(Model& model) : max_source_diff(model.parameters->max_source_diff) {
        // std::vector<Real> line_sources;
        // model_ptr = std::make_shared(model);
        line_sources.reserve(model.lines.emissivity.vec.size());
        std::transform(model.lines.emissivity.vec.begin(), model.lines.emissivity.vec.end(), model.lines.opacity.vec.begin(), line_sources.begin(), [](Real x, Real y) { return x/y; });

    }

    inline Size get_n_interp(Model& model, const Size curr_idx, const Size next_idx) const;
    inline std::vector<Real> get_subvector_sources(Model& model, const Size curr_idx) const;
    inline Real interpolate_linear(const Real f_start, const Real f_end, const Real factor) const;
    inline Real interpolate_log(const Real f_start, const Real f_end, const Real factor) const;

    // ///  Get the number of points to interpolate between two points
    // ///    @param[in] curr_idx : index of the current point
    // ///    @param[in] next_idx : index of the next point
    // ///    @returns number of points to interpolate
    // /////////////////////////////////////////////////
    // /// TODO: currently, this will result in probably too much interpolation, as the resulting grid is uniform for all lines/frequencies.
    // inline Size get_n_interpl(const Size curr_idx, const Size next_idx) const {
    //     std::vector<Real> curr_sources = get_subvector_sources(curr_idx);
    //     std::vector<Real> next_sources = get_subvector_sources(next_idx);
    //     // get relative difference between sources
    //     std::vector<Real> sources_diff = next_sources / curr_sources;
    //     // convert to ln, and take absolute value
    //     std::transform(sources_diff.begin(), sources_diff.end(), sources_diff.begin(), [](Real x) { return fabs(log(x)); });
    //     const Size source_diff_n_interp = std::ceil(*std::max_element(sources_diff.begin(), sources_diff.end())/std::log(max_source_diff));

    //     if (source_diff_n_interp > 1) {
    //         return source_diff_n_interp;
    //     } else {
    //         return 1;
    //     }
    // }

    // inline std::vector<Real> get_subvector_sources(const Size curr_idx) const {
    //     Size curr_line_start_idx = model.lines.get_index(curr_idx, 0);//assumes all data from any point lies next to eachother in the vector
    //     Size curr_line_end_idx = model.lines.get_index(curr_idx + 1, 0);

    //     return std::vector<Real>(self.line_sources.begin() + curr_line_start_idx, self.line_sources.begin() + curr_line_end_idx);
    // }
};

#include "solver_interp_helper.tpp"