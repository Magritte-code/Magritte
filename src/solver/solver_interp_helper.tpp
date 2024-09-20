
///  Get the number of points to interpolate between two points
///    @param[in] curr_idx : index of the current point
///    @param[in] next_idx : index of the next point
///    @returns number of points to interpolate
/////////////////////////////////////////////////
/// TODO: currently, this will probably result in too much interpolation, as the resulting grid is uniform for all lines/frequencies.
inline Size InterpHelper::get_n_interp(Model& model, const Size curr_idx, const Size next_idx) const {
    std::vector<Real> curr_sources = get_subvector_sources(model, curr_idx);
    std::vector<Real> next_sources = get_subvector_sources(model, next_idx);
    // get relative difference between sources
    std::vector<Real> sources_diff;
    sources_diff.reserve(curr_sources.size());
    std::transform(next_sources.begin(), next_sources.end(), curr_sources.begin(), sources_diff.begin(), [](Real x, Real y) { return x/y; });
    // convert to ln, and take absolute value
    std::transform(sources_diff.begin(), sources_diff.end(), sources_diff.begin(), [](Real x) { return fabs(log(x)); });
    const Size source_diff_n_interp = std::ceil(*std::max_element(sources_diff.begin(), sources_diff.end())/std::log(max_source_diff));

    if (source_diff_n_interp > 1) {
        return source_diff_n_interp;
    } else {
        return 1;
    }
}

inline std::vector<Real> InterpHelper::get_subvector_sources(Model& model, const Size curr_idx) const {
    Size curr_line_start_idx = model.lines.index(curr_idx, 0);//assumes all data from any point lies next to eachother in the vector
    Size curr_line_end_idx = model.lines.index(curr_idx + 1, 0);

    return std::vector<Real>(line_sources.begin() + curr_line_start_idx, line_sources.begin() + curr_line_end_idx);
}

///  Linear interpolation of f(x) in interval [0,1]
inline Real InterpHelper::interpolate_linear(const Real f_start, const Real f_end, const Real factor) const {
    return (x_end - x_start) * factor + x_start;
}

///  Logarithmic interpolation of f(x) in interval [0,1]
inline Real InterpHelper::interpolate_log(const Real f_start, const Real f_end, const Real factor) const {
    return exp((log(x_end) - log(x_start)) * factor + log(x_start));
}