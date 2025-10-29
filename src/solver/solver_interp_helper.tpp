
///  Get the number of points to interpolate between two points (inefficient upper bound)
///    @param[in] model : reference to the model
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
    Size n_interp_points = std::ceil(
        *std::max_element(sources_diff.begin(), sources_diff.end()) / logf(max_source_diff));

    // if dust is present, also consider dust opacity difference; upper bound given by iterating
    // over all dust frequencies
    if (model.dust.n_dust_frequencies > 0) {
        std::vector<Real> dust_opacity_diff;
        dust_opacity_diff.resize(model.dust.n_dust_frequencies);

        std::transform(model.dust.dust_opacities.dat + next_idx * model.dust.n_dust_frequencies,
            model.dust.dust_opacities.dat + (next_idx + 1) * model.dust.n_dust_frequencies,
            model.dust.dust_opacities.dat + curr_idx * model.dust.n_dust_frequencies,
            dust_opacity_diff.begin(), [min_op = model.parameters->min_opacity](Real x, Real y) {
                return fabs(logf((x + min_op) / (y + min_op)));
            });

        Size dust_diff_n_interp =
            std::ceil(*std::max_element(dust_opacity_diff.begin(), dust_opacity_diff.end())
                      / logf(max_source_diff)); // TODO: check if other value might be needed
        n_interp_points = std::max(dust_diff_n_interp, n_interp_points);
    }

    if (n_interp_points > 1) {
        return n_interp_points;
    } else {
        return 1;
    }
}

///  Get the number of points to interpolate between two points; internal use only, use
///  get_n_interp_around_freqs instead
///    @param[in] model : reference to the model
///    @param[in] l : line index
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

///  Get the number of points to interpolate between two points for dust opacity; internal use only,
///  use get_n_interp_around_freqs instead
///    @param[in] model : reference to the model
///    @param[in] curr_idx : index of the current point
///    @param[in] next_idx : index of the next point
///    @param[in] curr_freq : current frequency
///    @param[in] next_freq : next frequency
///    @returns the number of interpolation points needed
/// NOTE: use only if dust is present
inline Size InterpHelper::get_n_interp_for_dust(const Model& model, const Size curr_idx,
    const Size next_idx, const Real curr_freq, const Real next_freq) const {
    Size n_interpolation_points = 1;
    // evaluate dust opacity, interpolating on the grid
    Real curr_dust_opacity;
    Real next_dust_opacity;
    // If frequency is outside the range of the dust opacities, return the boundary value
    if (curr_freq < model.dust.dust_frequencies[0]) {
        curr_dust_opacity = model.dust.dust_opacities(curr_idx, 0);
    }
    if (curr_freq > model.dust.dust_frequencies[model.dust.n_dust_frequencies - 1]) {
        curr_dust_opacity = model.dust.dust_opacities(curr_idx, model.dust.n_dust_frequencies - 1);
    }
    // find the right bin using binary search
    Size f_low  = 0;
    Size f_high = model.dust.n_dust_frequencies - 1;
    while (f_high - f_low > 1) {
        Size f_mid = (f_high + f_low) / 2;
        if (model.dust.dust_frequencies[f_mid] > curr_freq) {
            f_high = f_mid;
        } else {
            f_low = f_mid;
        }
    }

    // and do linear interpolation
    Real slope =
        (model.dust.dust_opacities(curr_idx, f_high) - model.dust.dust_opacities(curr_idx, f_low))
        / (model.dust.dust_frequencies[f_high] - model.dust.dust_frequencies[f_low]);
    curr_dust_opacity = model.dust.dust_opacities(curr_idx, f_low)
                      + slope * (curr_freq - model.dust.dust_frequencies[f_low]);

    // Exactly the same for the next position
    if (next_freq < model.dust.dust_frequencies[0]) {
        next_dust_opacity = model.dust.dust_opacities(next_idx, 0);
    }
    if (next_freq > model.dust.dust_frequencies[model.dust.n_dust_frequencies - 1]) {
        next_dust_opacity = model.dust.dust_opacities(next_idx, model.dust.n_dust_frequencies - 1);
    }

    // find the right bin using binary search
    f_low  = 0;
    f_high = model.dust.n_dust_frequencies - 1;
    while (f_high - f_low > 1) {
        Size f_mid = (f_high + f_low) / 2;
        if (model.dust.dust_frequencies[f_mid] > next_freq) {
            f_high = f_mid;
        } else {
            f_low = f_mid;
        }
    }

    // and do linear interpolation
    slope =
        (model.dust.dust_opacities(next_idx, f_high) - model.dust.dust_opacities(next_idx, f_low))
        / (model.dust.dust_frequencies[f_high] - model.dust.dust_frequencies[f_low]);
    next_dust_opacity = model.dust.dust_opacities(next_idx, f_low)
                      + slope * (next_freq - model.dust.dust_frequencies[f_low]);

    Real dust_opacity_diff =
        std::abs(logf((next_dust_opacity + model.parameters->min_opacity)
                      / (curr_dust_opacity + model.parameters->min_opacity))); // in log space
    const Size n_interp = std::ceil(
        dust_opacity_diff / logf(max_source_diff)); // TODO: check if other value might be needed

    if (n_interp > 1) {
        return n_interp;
    } else {
        return 1;
    }
}

///  Get the number of points to interpolate between two points around certain frequencies
///    @param[in] model : reference to the model
///    @param[in] curr_idx : index of the current point
///    @param[in] next_idx : index of the next point
///    @param[in] curr_freq : current frequency
///    @param[in] next_freq : next frequency
///    @returns the number of interpolation points needed
inline Size InterpHelper::get_n_interp_around_freqs(const Model& model, const Size curr_idx,
    const Size next_idx, const Real curr_freq, const Real next_freq) const {

    Size n_interpolation_points = 1;

    // first try to find which lines are close enough to any of the two frequencies
    // and use the corresponding interpolation criteria
    // copied from Solver::compute_S_dtau_line_integrated<CloseLines> in solver.tpp
    // determine frequency bounds for searching nearby
    Real left_freq;
    Real right_freq;

    left_freq  = std::min(curr_freq, next_freq);
    right_freq = std::max(curr_freq, next_freq);

    // using maximum of bounds on the two points to get an
    // upper bound for the line width
    const Real curr_bound_line_width = model.parameters->max_distance_opacity_contribution
                                     * model.thermodynamics.profile_width_upper_bound_with_linefreq(
                                         curr_idx, right_freq, model.lines.max_inverse_mass);
    const Real next_bound_line_width = model.parameters->max_distance_opacity_contribution
                                     * model.thermodynamics.profile_width_upper_bound_with_linefreq(
                                         next_idx, right_freq, model.lines.max_inverse_mass);
    const Real upper_bound_line_width = std::max(curr_bound_line_width, next_bound_line_width);

    const Real left_freq_bound  = left_freq - upper_bound_line_width;
    const Real right_freq_bound = right_freq + upper_bound_line_width;

    // apply default search algorithms on the bounds,
    // obtaining iterators
    auto left_line_bound = std::lower_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound = std::upper_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);

    for (auto freq_sort_l = left_line_bound; freq_sort_l != right_line_bound; freq_sort_l++) {
        const Size sort_l = freq_sort_l - model.lines.sorted_line.begin();
        // Map sorted line index to original line index
        const Size l = model.lines.sorted_line_map[sort_l];

        // evaluate interpolation criterion for this line, take max
        n_interpolation_points =
            std::max(get_n_interp_for_line(model, l, curr_idx, next_idx), n_interpolation_points);
    }

    // and also include the dust opacity in the interpolation if present
    if (model.dust.n_dust_frequencies > 0) {
        n_interpolation_points =
            std::max(get_n_interp_for_dust(model, curr_idx, next_idx, curr_freq, next_freq),
                n_interpolation_points);
    }

    return n_interpolation_points;
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