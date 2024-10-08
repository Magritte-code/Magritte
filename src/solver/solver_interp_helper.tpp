
///  Get the number of points to interpolate between two points
///    @param[in] curr_idx : index of the current point
///    @param[in] next_idx : index of the next point
///    @returns number of points to interpolate
/////////////////////////////////////////////////
/// Note: this will result in too much interpolation, as not all frequencies around each line need
/// the same interpolation points uniform for all lines/frequencies. Evidently, this is an upper
/// bound for get_n_interp_for_line.
inline Size InterpHelper::get_n_interp(
    const Model& model, const Size curr_idx, const Size next_idx) const {

    Size curr_line_start_idx = model.lines.index(
        curr_idx, 0); // assumes all data from any point lies next to eachother in the vector
    Size curr_line_end_idx = model.lines.index(curr_idx + 1, 0);

    Size next_line_start_idx = model.lines.index(
        next_idx, 0); // assumes all data from any point lies next to eachother in the vector
    Size next_line_end_idx = model.lines.index(next_idx + 1, 0);

    // get relative difference between sources
    std::vector<Real> sources_diff;
    sources_diff.resize(model.parameters->nlines());

    // in linear space
    // std::transform(interpolation_criterion.begin() + next_line_start_idx,
    //     interpolation_criterion.begin() + next_line_end_idx,
    //     interpolation_criterion.begin() + curr_line_start_idx, sources_diff.begin(),
    //     [](Real x, Real y) { return std::abs(x - y) / std::min(x, y); });
    // const Size source_diff_n_interp = std::ceil(
    //     *std::max_element(sources_diff.begin(), sources_diff.end()) / (max_source_diff - 1.0));

    // in log space
    std::transform(interpolation_criterion.begin() + next_line_start_idx,
        interpolation_criterion.begin() + next_line_end_idx,
        interpolation_criterion.begin() + curr_line_start_idx, sources_diff.begin(),
        [](Real x, Real y) { return fabs(logf(x / y)); });
    const Size source_diff_n_interp = std::ceil(
        *std::max_element(sources_diff.begin(), sources_diff.end()) / logf(max_source_diff));

    if (source_diff_n_interp > 1) {
        return source_diff_n_interp;
    } else {
        return 1;
    }
}

///  Get the number of points to interpolate between two points
///    @param[in] curr_idx : index of the current point
///    @param[in] next_idx : index of the next point
///    @returns number of points to interpolate
/////////////////////////////////////////////////
/// Note: current interpolation assumes that every line is fully seperated from the others, which
/// might not be the case for overlapping lines. FIXME: add better interpolation criterion
inline Size InterpHelper::get_n_interp_for_line(
    const Model& model, const Size l, const Size curr_idx, const Size next_idx) const {
    Real curr_source = interpolation_criterion[model.lines.index(curr_idx, l)];
    Real next_source = interpolation_criterion[model.lines.index(next_idx, l)];

    // in linear space
    // Real source_diff                = std::abs(next_source - curr_source);
    // Real min_source                 = std::min(next_source, curr_source);
    // const Size source_diff_n_interp = std::ceil(source_diff / (min_source * (max_source_diff -
    // 1)));

    // in log space
    Real source_diff                = std::abs(logf(next_source / curr_source));
    const Size source_diff_n_interp = std::ceil(source_diff / logf(max_source_diff));

    if (source_diff_n_interp > 1) {
        return source_diff_n_interp;
    } else {
        return 1;
    }
}

///  Linear interpolation of f(x) in interval [0,1]
inline Real InterpHelper::interpolate_linear(
    const Real f_start, const Real f_end, const Real factor) const {
    if (factor == 0.0) {
        return f_start;
    } else if (factor == 1.0) {
        return f_end;
    }
    return (f_end - f_start) * factor + f_start;
}

///  Logarithmic interpolation of f(x) in interval [0,1]
inline Real InterpHelper::interpolate_log(
    const Real f_start, const Real f_end, const Real factor) const {
    // log is way too expensive to use without essentially doing anything
    if (factor == 0.0) { //<- branch prediction should be here most of the time; assumes
                         // interpolation to only be used for a small part of the model
        return f_start;
    } else if (factor == 1.0) {
        return f_end;
    }
    return f_start * powf(f_end / f_start, factor);
}