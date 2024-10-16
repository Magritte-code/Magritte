template <Frame frame, bool use_adaptive_directions> inline void Solver ::setup(Model& model) {
    interp_helper = InterpHelper(model);
    // ray_lengths.resize(model.parameters->hnrays(), model.parameters->npoints());
    const Size length = 2 * get_ray_lengths_max<frame, use_adaptive_directions>(model) + 1;
    const Size width  = model.parameters->nfreqs();
    const Size n_o_d  = model.parameters->n_off_diag;
    std::cout << "length: " << length << std::endl;

    model.set_dshift_max(); // err, probably belongs somewhere else, but we need to compute the max
                            // shift for each point

    setup(length, width, n_o_d);
}

// should be called inside the imager, as these things might
// not be fixed (see adaptive imager)
//  template <Frame frame>
inline void Solver ::setup_new_imager(Model& model, Image& image, const Vector3D& ray_dir) {
    interp_helper = InterpHelper(model);
    // ray_lengths.resize(model.parameters->hnrays(), model.parameters->npoints());
    const Size length = 2 * get_ray_lengths_max_new_imager(model, image, ray_dir) + 1;
    const Size width  = model.parameters->nfreqs();
    const Size n_o_d  = model.parameters->n_off_diag;
    std::cout << "length: " << length << std::endl;

    model.set_dshift_max(); // err, probably belongs somewhere else, but we need to compute the max
                            // shift for each point

    setup(length, width, n_o_d);
}

inline void Solver ::setup(const Size l, const Size w, const Size n_o_d) {
    length     = l;
    centre     = l / 2;
    width      = w;
    n_off_diag = n_o_d;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++) {
        dZ_(i).resize(length);
        nr_(i).resize(length);
        nr_interp_(i).resize(length);
        interp_factor_(i).resize(length);
        shift_(i).resize(length);

        eta_c_(i).resize(width);
        eta_n_(i).resize(width);

        chi_c_(i).resize(width);
        chi_n_(i).resize(width);

        source_c_(i).resize(width);
        source_n_(i).resize(width);

        inverse_chi_(i).resize(length);

        tau_(i).resize(width);

        Su_(i).resize(length);
        Sv_(i).resize(length);

        A_(i).resize(length);
        C_(i).resize(length);
        inverse_A_(i).resize(length);
        inverse_C_(i).resize(length);

        FF_(i).resize(length);
        FI_(i).resize(length);
        GG_(i).resize(length);
        GI_(i).resize(length);
        GP_(i).resize(length);

        L_diag_(i).resize(length);

        L_upper_(i).resize(n_off_diag, length);
        L_lower_(i).resize(n_off_diag, length);
    }
}

template <Frame frame, bool use_adaptive_directions>
accel inline Size Solver ::get_ray_length(Model& model, const Size o, const Size r) const {
    Size l    = 0;   // ray length
    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt)) {
        Size crt = o;
        l += interp_helper.get_n_interp(model, crt, nxt);

        while (model.geometry.not_on_boundary(nxt)) {
            crt = nxt;
            nxt = model.geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);
            l += interp_helper.get_n_interp(model, crt, nxt);

            if (!model.geometry.valid_point(nxt)) {
                printf("ERROR: no valid neighbor o=%u, r=%u, crt=%u\n", o, r, crt);
            }
        }
    }

    return l;
}

template <Frame frame, bool use_adaptive_directions>
inline Size Solver ::get_ray_lengths_max(Model& model) {
    Matrix<Size> ray_lengths;
    ray_lengths.resize(model.parameters->hnrays(), model.parameters->npoints());
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        accelerated_for(o, model.parameters->npoints(), {
            const Size ar = model.geometry.rays.get_antipod_index(rr);

            ray_lengths(rr, o) = get_ray_length<frame, use_adaptive_directions>(model, o, rr)
                               + get_ray_length<frame, use_adaptive_directions>(model, o, ar);
        })

        pc::accelerator::synchronize();
    }

    ray_lengths.copy_ptr_to_vec();

    Size lengths_max = *std::max_element(ray_lengths.vec.begin(), ray_lengths.vec.end());

    return lengths_max;
}

inline Size Solver ::get_ray_lengths_max_new_imager(
    Model& model, Image& image, const Vector3D& ray_dir) {
    const Size npixels = image.ImX.size(); // is ImY.size() (number of pixels in image)

    std::vector<Size> ray_lengths;
    ray_lengths.resize(npixels); // TODO: later on, define which pixels need to
                                 // be used for raytracing (adaptive imager)/or
                                 // define new partial images (might be more
                                 // elegant)

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        ray_lengths[pixidx] = get_ray_length_new_imager(model, origin, start_bdy_point, ray_dir);
    });

    const Size max_ray_length = (*std::max_element(ray_lengths.begin(), ray_lengths.end()));
    return max_ray_length;
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_shortchar_order_0(Model& model) {
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities) {
        model.radiation.I.resize(
            model.parameters->nrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.u.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize(model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();
    }

    // Initialise mean intensity
    model.radiation.initialize_J();

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        for (Size lspec_idx = 0; lspec_idx < model.lines.lineProducingSpecies.size(); lspec_idx++) {
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[lspec_idx];
            accelerated_for(o, model.parameters->npoints(), {
                const Real wt =
                    model.geometry.rays.get_weight<use_adaptive_directions>(o, rr) * two;
                const Size ar = model.geometry.rays.get_antipod_index(rr);

                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                for (Size k = 0; k < lspec.linedata.nrad; k++) {
                    const Size l = model.lines.line_index(lspec_idx, k);
                    // For every individual line, the required interpolation can differ
                    first_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                   model, l, o, rr, -1, centre - 1, centre - 1)
                             + 1;
                    last_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                  model, l, o, ar, +1, centre + 1, centre)
                            - 1;
                    n_tot_() = (last_() + 1) - first_();
                    // std::cout << "o: " << o << " n_tot_ = " << n_tot_() << std::endl;

                    if (n_tot_() > 1) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            const Size f   = lspec.nr_line[o][k][z];
                            const Real I_f = solve_shortchar_order_0_ray_forward<approx,
                                use_adaptive_directions>(model, o, rr, f);
                            const Real I_b = solve_shortchar_order_0_ray_backward<approx,
                                use_adaptive_directions>(model, o, rr, f);

                            model.radiation.I(rr, o, f) = I_f;
                            model.radiation.I(ar, o, f) = I_b;
                            model.radiation.u(rr, o, f) = 0.5 * (I_f + I_b);
                            model.radiation.v(rr, o, f) = 0.5 * (I_f - I_b);
                            model.radiation.J(o, f) += wt * model.radiation.u(rr, o, f);
                        }
                    } else {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            const Size f = lspec.nr_line[o][k][z];
                            const Real bdy_I =
                                boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                            model.radiation.I(rr, o, f) = bdy_I;
                            model.radiation.I(ar, o, f) = bdy_I;
                            model.radiation.u(rr, o, f) = bdy_I;
                            model.radiation.v(rr, o, f) = 0.0;
                            model.radiation.J(o, f) += wt * bdy_I;
                        }
                    }
                }
            })
        }

        pc::accelerator::synchronize();
    }

    model.radiation.I.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2_uv(Model& model) {
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities) {
        model.radiation.u.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
    }

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        const Size ar = model.geometry.rays.get_antipod_index(rr);
        for (Size lspec_idx = 0; lspec_idx < model.parameters->nlspecs(); lspec_idx++) {
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[lspec_idx];
            accelerated_for(o, model.parameters->npoints(), {
                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                for (Size k = 0; k < lspec.linedata.nrad; k++) {
                    const Size l = model.lines.line_index(lspec_idx, k);
                    // For every individual line, the required interpolation can differ
                    first_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                   model, l, o, rr, -1, centre - 1, centre - 1)
                             + 1;
                    last_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                  model, l, o, ar, +1, centre + 1, centre)
                            - 1;
                    n_tot_() = (last_() + 1) - first_();

                    if (n_tot_() > 1) {
                        for (Size q = 0; q < model.parameters->nquads(); q++) {
                            const Size f = lspec.nr_line[o][k][q];
                            solve_feautrier_order_2_uv<approx>(model, o, f);

                            model.radiation.u(rr, o, f) = Su_()[centre];
                            model.radiation.v(rr, o, f) = Sv_()[centre];
                        }
                    } else {
                        for (Size q = 0; q < model.parameters->nquads(); q++) {
                            const Size f = lspec.nr_line[o][k][q];
                            model.radiation.u(rr, o, f) =
                                boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                            model.radiation.v(rr, o, f) = 0.0;
                        }
                    }
                }
            })
        }

        pc::accelerator::synchronize();
    }

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.v.copy_ptr_to_vec();
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2_sparse(Model& model) {
    // Initialise variables
    for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for(o, model.parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                lspec.J(o, k) = 0.0;
            }
        })
    }

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        for (Size lspec_idx = 0; lspec_idx < model.lines.lineProducingSpecies.size(); lspec_idx++) {
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[lspec_idx];
            threaded_for(o, model.parameters->npoints(), {
                const Vector3D nn =
                    model.geometry.rays.get_direction<use_adaptive_directions>(o, rr);
                const Size ar = model.geometry.rays.get_antipod_index(rr);
                const Real wt =
                    model.geometry.rays.get_weight<use_adaptive_directions>(o, rr) * two;

                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                for (Size k = 0; k < lspec.linedata.nrad; k++) {
                    const Size l = model.lines.line_index(lspec_idx, k);
                    // For every individual line, the required interpolation can differ
                    first_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                   model, l, o, rr, -1, centre - 1, centre - 1)
                             + 1;
                    last_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                  model, l, o, ar, +1, centre + 1, centre)
                            - 1;
                    n_tot_() = (last_() + 1) - first_();

                    if (n_tot_() > 1) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            solve_feautrier_order_2<approx>(model, o, lspec.nr_line[o][k][z]);

                            lspec.J(o, k) += lspec.quadrature.weights[z] * wt * Su_()[centre];

                            update_Lambda<approx, use_adaptive_directions>(
                                model, rr, lspec.nr_line[o][k][z]);
                        }
                    } else {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            lspec.J(o, k) +=
                                lspec.quadrature.weights[z] * wt
                                * boundary_intensity(model, o,
                                    model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                        }
                    }
                }
            })
        }
    }
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2_anis(Model& model) {
    // Initialise variables
    for (LineProducingSpecies& lspec : model.lines.lineProducingSpecies) {
        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_0.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for(o, model.parameters->npoints(), {
            for (Size k = 0; k < lspec.linedata.nrad; k++) {
                lspec.J(o, k)       = 0.0;
                lspec.J2_0(o, k)    = 0.0;
                lspec.J2_1_Re(o, k) = 0.0;
                lspec.J2_1_Im(o, k) = 0.0;
                lspec.J2_2_Re(o, k) = 0.0;
                lspec.J2_2_Im(o, k) = 0.0;
            }
        })
    }

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        const Size ar = model.geometry.rays.get_antipod_index(rr);
        for (Size lspec_idx = 0; lspec_idx < model.parameters->nlspecs(); lspec_idx++) {
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[lspec_idx];

            threaded_for(o, model.parameters->npoints(), {
                const Real wt = model.geometry.rays.get_weight<use_adaptive_directions>(o, rr);
                const Vector3D nn =
                    model.geometry.rays.get_direction<use_adaptive_directions>(o, rr);

                const Real wt_0    = inv_sqrt2 * (three * nn.z() * nn.z() - one);
                const Real wt_1_Re = -sqrt3 * nn.x() * nn.z();
                const Real wt_1_Im = -sqrt3 * nn.y() * nn.z();
                const Real wt_2_Re = half * sqrt3 * (nn.x() * nn.x() - nn.y() * nn.y());
                const Real wt_2_Im = sqrt3 * nn.x() * nn.y();

                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                for (Size k = 0; k < lspec.linedata.nrad; k++) {
                    const Size l = model.lines.line_index(lspec_idx, k);

                    first_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                   model, l, o, rr, -1, centre - 1, centre - 1)
                             + 1;
                    last_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                  model, l, o, ar, +1, centre + 1, centre)
                            - 1;
                    n_tot_() = (last_() + 1) - first_();

                    if (n_tot_() > 1) {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            solve_feautrier_order_2<approx>(model, o, lspec.nr_line[o][k][z]);

                            const Real du = lspec.quadrature.weights[z] * wt * Su_()[centre];

                            lspec.J(o, k) += two * du;
                            lspec.J2_0(o, k) += wt_0 * du;
                            lspec.J2_1_Re(o, k) += wt_1_Re * du;
                            lspec.J2_1_Im(o, k) += wt_1_Im * du;
                            lspec.J2_2_Re(o, k) += wt_2_Re * du;
                            lspec.J2_2_Im(o, k) += wt_2_Im * du;
                        }
                    } else {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++) {
                            const Real du =
                                lspec.quadrature.weights[z] * wt
                                * boundary_intensity(model, o,
                                    model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));

                            lspec.J(o, k) += two * du;
                            lspec.J2_0(o, k) += wt_0 * du;
                            lspec.J2_1_Re(o, k) += wt_1_Re * du;
                            lspec.J2_1_Im(o, k) += wt_1_Im * du;
                            lspec.J2_2_Re(o, k) += wt_2_Re * du;
                            lspec.J2_2_Im(o, k) += wt_2_Im * du;
                        }
                    }
                }
            })
        }
    }
}

template <ApproximationType approx, bool use_adaptive_directions>
inline void Solver ::solve_feautrier_order_2(Model& model) {
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities) {
        model.radiation.u.resize(
            model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize(model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.clear();
    }

    // Initialise mean intensity
    model.radiation.initialize_J();

    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {
        cout << "--- rr = " << rr << endl;

        const Size ar = model.geometry.rays.get_antipod_index(rr);
        for (Size lspec_idx = 0; lspec_idx < model.parameters->nlspecs(); lspec_idx++) {
            LineProducingSpecies& lspec = model.lines.lineProducingSpecies[lspec_idx];

            accelerated_for(o, model.parameters->npoints(), {
                nr_()[centre]    = o;
                shift_()[centre] = 1.0;

                for (Size k = 0; k < lspec.linedata.nrad; k++) {
                    const Size l = model.lines.line_index(lspec_idx, k);
                    // For every individual line, the required interpolation can differ
                    first_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                   model, l, o, rr, -1, centre - 1, centre - 1)
                             + 1;
                    last_() = trace_ray_for_line<CoMoving, use_adaptive_directions>(
                                  model, l, o, ar, +1, centre + 1, centre)
                            - 1;
                    n_tot_() = (last_() + 1) - first_();

                    if (n_tot_() > 1) {
                        for (Size q = 0; q < model.parameters->nquads(); q++) {
                            const Size f = lspec.nr_line[o][k][q];
                            solve_feautrier_order_2<approx>(model, o, f);

                            model.radiation.u(rr, o, f) = Su_()[centre];
                            model.radiation.J(o, f) +=
                                Su_()[centre] * two
                                * model.geometry.rays.get_weight<use_adaptive_directions>(o, rr);

                            update_Lambda<approx, use_adaptive_directions>(model, rr, f);
                        }
                    } else {
                        for (Size q = 0; q < model.parameters->nquads(); q++) {
                            const Size f = lspec.nr_line[o][k][q];
                            model.radiation.u(rr, o, f) =
                                boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                            model.radiation.J(o, f) +=
                                two * model.geometry.rays.get_weight<use_adaptive_directions>(o, rr)
                                * model.radiation.u(rr, o, f);
                        }
                    }
                }
            })
        }

        pc::accelerator::synchronize();
    }

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();

    cout << "MPI gathering J..." << endl;
    model.radiation.MPI_reduce_J();
    cout << "Done MPI gathering J." << endl;

    cout << "MPI gathering Lambda..." << endl;
    // Gather contributions to the ALO
    for (auto& lspec : model.lines.lineProducingSpecies) {
        lspec.lambda.MPI_gather();
    }
    cout << "Done MPI gathering Lambda." << endl;
}

template <ApproximationType approx>
inline void Solver ::image_feautrier_order_2(Model& model, const Size rr) {
    Image image = Image(model.geometry, model.radiation.frequencies, Intensity, rr);

    accelerated_for(o, model.parameters->npoints(), {
        const Size ar = model.geometry.rays.get_antipod_index(rr);

        nr_()[centre]    = o;
        shift_()[centre] = model.geometry.get_shift<Rest, false>(o, rr, o, 0.0);

        first_() = trace_ray<Rest, false>(model, o, rr, -1, centre - 1, centre - 1) + 1;
        last_()  = trace_ray<Rest, false>(model, o, ar, +1, centre + 1, centre) - 1;
        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_feautrier_order_2<approx>(model, o, f);

                image.I(o, f) = two * Su_()[last_()]
                              - boundary_intensity(
                                  model, nr_()[last_()], model.radiation.frequencies.nu(o, f));
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(o, f) = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

template <ApproximationType approx>
inline void Solver ::image_feautrier_order_2_new_imager(
    Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix) {
    Image image =
        Image(model.geometry, model.radiation.frequencies, Intensity, ray_dir, nxpix, nypix);
    setup_new_imager(model, image, ray_dir);

    // Note: number of pixels is constant for now, but may be
    // adaptive in the future; (but then while loop will be
    // required anyway)
    const Size npixels             = image.ImX.size(); // is ImY.size(), is I.size()
    const Vector3D origin_velocity = Vector3D(0.0);

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        Real Z = 0.0;
        const Size closest_bdy_point =
            trace_ray_imaging_get_start(model.geometry, origin, start_bdy_point, ray_dir, Z);

        nr_()[centre]    = closest_bdy_point;
        shift_()[centre] = model.geometry.get_shift<Rest>(
            origin, origin_velocity, ray_dir, closest_bdy_point, Z, false);
        first_() = trace_ray_imaging<Rest>(
                       model, origin, closest_bdy_point, ray_dir, -1, Z, centre - 1, centre - 1)
                 + 1;
        last_() = centre; // by definition, only boundary points
                          // can lie in the backward direction

        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_feautrier_order_2<approx>(model, closest_bdy_point, f);
                image.I(pixidx, f) = two * Su_()[last_()]
                                   - boundary_intensity(model, nr_()[last_()],
                                       model.radiation.frequencies.nu(closest_bdy_point, f));
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(pixidx, f) = boundary_intensity(
                    model, closest_bdy_point, model.radiation.frequencies.nu(closest_bdy_point, f));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

template <ApproximationType approx>
inline void Solver ::image_feautrier_order_2_for_point(Model& model, const Size rr, const Size p) {
    // Redefine p to keep everything as similar to
    // image_feautrier_order_2 as possible
    const Size o = p;

    const Size ar = model.geometry.rays.get_antipod_index(rr);

    nr_()[centre]    = o;
    shift_()[centre] = model.geometry.get_shift<Rest, false>(o, rr, o, 0.0);

    first_() = trace_ray<Rest, false>(model, o, rr, -1, centre - 1, centre - 1) + 1;
    last_()  = trace_ray<Rest, false>(model, o, ar, +1, centre + 1, centre) - 1;
    n_tot_() = (last_() + 1) - first_();

    model.S_ray.resize(n_tot_(), model.parameters->nfreqs());
    model.dtau_ray.resize(n_tot_(), model.parameters->nfreqs());
    model.u_ray.resize(n_tot_(), model.parameters->nfreqs());

    if (n_tot_() > 1) {
        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            image_feautrier_order_2_for_point_loc<approx>(model, o, f);
        }
    }
}

template <ApproximationType approx>
inline void Solver ::image_optical_depth(Model& model, const Size rr) {
    Image image = Image(model.geometry, model.radiation.frequencies, Intensity, rr);

    accelerated_for(o, model.parameters->npoints(), {
        const Size ar = model.geometry.rays.get_antipod_index(rr);

        nr_()[centre]    = o;
        shift_()[centre] = model.geometry.get_shift<Rest, false>(o, rr, o, 0.0);

        first_() = trace_ray<Rest, false>(model, o, rr, -1, centre - 1, centre - 1) + 1;
        last_()  = trace_ray<Rest, false>(model, o, ar, +1, centre + 1, centre) - 1;
        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_optical_depth<approx>(model, o, f);

                image.I(o, f) = optical_depth_();
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(o, f) = 0.0;
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

template <ApproximationType approx>
inline void Solver ::image_optical_depth_new_imager(
    Model& model, const Vector3D& ray_dir, const Size nxpix, const Size nypix) {
    Image image =
        Image(model.geometry, model.radiation.frequencies, OpticalDepth, ray_dir, nxpix, nypix);
    setup_new_imager(model, image, ray_dir);

    // Note: number of pixels is constant for now, but may be
    // adaptive in the future; (but then while loop will be
    // required anyway)
    const Size npixels             = image.ImX.size(); // is ImY.size(), is I.size()
    const Vector3D origin_velocity = Vector3D(0.0);

    const Size start_bdy_point = model.geometry.get_closest_bdy_point_in_custom_raydir(ray_dir);

    accelerated_for(pixidx, npixels, {
        const Vector3D origin =
            image.surface_coords_to_3D_coordinates(image.ImX[pixidx], image.ImY[pixidx]);
        Real Z = 0.0;
        const Size closest_bdy_point =
            trace_ray_imaging_get_start(model.geometry, origin, start_bdy_point, ray_dir, Z);

        nr_()[centre]    = closest_bdy_point;
        shift_()[centre] = model.geometry.get_shift<Rest>(
            origin, origin_velocity, ray_dir, closest_bdy_point, Z, false);
        first_() = trace_ray_imaging<Rest>(
                       model, origin, closest_bdy_point, ray_dir, -1, Z, centre - 1, centre - 1)
                 + 1;
        last_() = centre; // by definition, only boundary points
                          // can lie in the backward direction

        n_tot_() = (last_() + 1) - first_();

        if (n_tot_() > 1) {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image_optical_depth<approx>(model, closest_bdy_point, f);

                image.I(pixidx, f) = optical_depth_();
            }
        } else {
            for (Size f = 0; f < model.parameters->nfreqs(); f++) {
                image.I(pixidx, f) = 0.0;
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back(image);
}

// Traces ray and sets data required for the solvers
// Does this for all lines at the same time, so the required interpolation might be
// excessive
template <Frame frame, bool use_adaptive_directions>
accel inline Size Solver ::trace_ray(
    const Model& model, const Size o, const Size r, const int increment, Size id1, Size id2) {
    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt)) {
        Size crt         = o;
        double shift_crt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, crt, 0.0);
        double shift_nxt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

        set_data(model, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);

        while (model.geometry.not_on_boundary(nxt)) {
            crt       = nxt;
            shift_crt = shift_nxt;

            nxt       = model.geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);
            shift_nxt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

            set_data(model, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);
        }
    }

    return id1;
}

// Traces ray and sets data required for the solvers
// Note: only does this for a single line at a time, as the interpolation might be different
// for the different lines
template <Frame frame, bool use_adaptive_directions>
accel inline Size Solver ::trace_ray_for_line(const Model& model, const Size l, const Size o,
    const Size r, const int increment, Size id1, Size id2) {
    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt)) {
        Size crt         = o;
        double shift_crt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, crt, 0.0);
        double shift_nxt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

        set_data_for_line(model, l, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);

        while (model.geometry.not_on_boundary(nxt)) {
            crt       = nxt;
            shift_crt = shift_nxt;

            nxt       = model.geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);
            shift_nxt = model.geometry.get_shift<frame, use_adaptive_directions>(o, r, nxt, Z);

            set_data_for_line(model, l, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);
        }
    }

    return id1;
}

/// Obtain the starting point for the imaging raytracer
/// @param[in] geometry: The geometry of the model
/// @param[in] origin: The origin of the ray (outside of the model)
/// @param[in] start_bdy: The starting boundary point for the search
/// @param[in] raydir: The direction of the ray
/// @param[out] Z: The distance from the origin to the boundary point
accel inline Size Solver ::trace_ray_imaging_get_start(const Geometry& geometry,
    const Vector3D& origin, const Size start_bdy, const Vector3D& raydir, Real& Z) {

    Size initial_point = start_bdy;
    // first figure out which boundary point lies closest to
    // the custom ray
    // TODO: is slightly inefficient implementation, can be
    // improved by only checking bdy point neighbors. We use
    // here the assumption that the outer boundary is convex,
    // to obtain the closest point on the boundary
    while (true) {
        Size next_attempt =
            geometry.get_boundary_point_closer_to_custom_ray(origin, raydir, initial_point);
        if (next_attempt == initial_point) {
            break;
        }
        initial_point = next_attempt;
    }

    Z = geometry.get_distance_origin_to_boundary(origin, raydir, initial_point);

    return initial_point;
}

/// Sets the data required for the imaging raytracer
/// @param[in] model: The model
/// @param[in] origin: The origin of the ray
/// @param[in] start_bdy: The starting boundary point
/// @param[in] raydir: The direction of the ray
/// @param[in] increment: The increment of the data index (-1 or +1)
/// @param[in] Z: The distance from the origin to the boundary point
/// @param[in/out] id1: The first data index (used for all data)
/// @param[in/out] id2: The second data index (only used for dZ)
template <Frame frame>
accel inline Size Solver ::trace_ray_imaging(const Model& model, const Vector3D& origin,
    const Size start_bdy, const Vector3D& raydir, const int increment,
    Real& Z, // distance from origin can be non-zero to start,
             // as this measures the distance from the
             // projection plane
    Size id1, Size id2) {

    double dZ                      = 0.0; // last increment in Z
    const Vector3D origin_velocity = Vector3D(0.0, 0.0, 0.0);
    Size crt                       = start_bdy;

    Size nxt = model.geometry.get_next<Imagetracer>(origin, raydir, start_bdy, Z, dZ);

    bool accessed_non_boundary = false;
    Size original_id1 =
        id1; // initial ray index to return in case of only boundary points on the ray

    if (model.geometry.valid_point(nxt)) {
        // Size crt = o;
        double shift_crt =
            model.geometry.get_shift<frame>(origin, origin_velocity, raydir, crt, Z - dZ, false);
        double shift_nxt =
            model.geometry.get_shift<frame>(origin, origin_velocity, raydir, nxt, Z, false);

        set_data(model, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);

        while (model.geometry.not_on_boundary(nxt)
               || (model.geometry.points.position[nxt] - model.geometry.points.center).dot(raydir)
                      < 0) {
            crt       = nxt;
            shift_crt = shift_nxt;

            if (model.geometry.not_on_boundary(nxt)) {
                accessed_non_boundary = true;
            }

            nxt = model.geometry.get_next<Imagetracer>(origin, raydir, nxt, Z, dZ);

            if (!model.geometry.valid_point(nxt)) {
                break;
            }

            shift_nxt =
                model.geometry.get_shift<frame>(origin, origin_velocity, raydir, nxt, Z, false);

            set_data(model, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);
        }
    }

    // Spherically symmetric models have a correct boundary treatment already, thus no need to
    // reduce their accuracy in a region between the outer boundary and the nearest point.
    if (!accessed_non_boundary && !model.parameters->spherical_symmetry()) {
        return original_id1;
    }
    return id1;
}

// Traces rays for the imager, and sets the data required
// Note: only does this for a single line at a time, as the interpolation might be different
// for the different lines
/// TODO: implement in the imagers, to save computation time
template <Frame frame>
accel inline Size Solver ::trace_ray_imaging_for_line(const Model& model, const Size l,
    const Vector3D& origin, const Size start_bdy, const Vector3D& raydir, const int increment,
    Real& Z, // distance from origin can be non-zero to start,
             // as this measures the distance from the
             // projection plane
    Size id1, Size id2) {

    double dZ                      = 0.0; // last increment in Z
    const Vector3D origin_velocity = Vector3D(0.0, 0.0, 0.0);
    Size crt                       = start_bdy;

    Size nxt = model.geometry.get_next<Imagetracer>(origin, raydir, start_bdy, Z, dZ);

    bool accessed_non_boundary = false;
    Size original_id1 =
        id1; // initial ray index to return in case of only boundary points on the rays

    if (model.geometry.valid_point(nxt)) {
        double shift_crt =
            model.geometry.get_shift<frame>(origin, origin_velocity, raydir, crt, Z - dZ, false);
        double shift_nxt = model.geometry.get_shift<frame>(origin, origin_velocity, raydir, nxt, Z);

        set_data_for_line(model, l, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);

        while (model.geometry.not_on_boundary(nxt)
               || (model.geometry.points.position[nxt] - model.geometry.points.center).dot(raydir)
                      < 0) {
            crt       = nxt;
            shift_crt = shift_nxt;

            if (model.geometry.not_on_boundary(nxt)) {
                accessed_non_boundary = true;
            }

            nxt = model.geometry.get_next<Imagetracer>(origin, raydir, nxt, Z, dZ);

            if (!model.geometry.valid_point(nxt)) {
                break;
            }

            shift_nxt = model.geometry.get_shift<frame>(origin, origin_velocity, raydir, nxt, Z);

            set_data_for_line(model, l, crt, nxt, shift_crt, shift_nxt, dZ, increment, id1, id2);
        }
    }

    if (!accessed_non_boundary) {
        return original_id1;
    }

    return id1;
}

/// Returns the length of a ray for the new imager
/// @param[in] model : model object
/// @param[in] origin : origin of the ray (outside of the model)
/// @param[in] start_bdy : starting boundary point
/// @param[in] raydir : direction of the ray
/// @return : length of the ray
/// Note: does not efficiently use the interpolation, as all lines will assumed to be relevant
accel inline Size Solver ::get_ray_length_new_imager(
    const Model& model, const Vector3D& origin, const Size start_bdy, const Vector3D& raydir) {
    Size l             = 0; // ray length, which we need to compute
    double Z           = 0.0;
    double dZ          = 0.0; // last increment in Z
    Size initial_point = trace_ray_imaging_get_start(model.geometry, origin, start_bdy, raydir, Z);
    Size crt           = initial_point;
    Size nxt           = model.geometry.get_next<Imagetracer>(origin, raydir, crt, Z, dZ);

    if (model.geometry.valid_point(nxt)) {
        l += interp_helper.get_n_interp(model, crt, nxt);

        while (model.geometry.not_on_boundary(nxt)
               || (model.geometry.points.position[nxt] - model.geometry.points.center).dot(raydir)
                      < 0) {
            crt = nxt;
            nxt = model.geometry.get_next<Imagetracer>(origin, raydir, nxt, Z, dZ);
            if (!model.geometry.valid_point(nxt)) {
                break;
            }
            l += interp_helper.get_n_interp(model, crt, nxt);
        }
    }

    return l;
}

/// Sets data for a ray, including required interpolation for the computations.
/// @param[in] model : model object
/// @param[in] crt : current point
/// @param[in] nxt : next point
/// @param[in] shift_crt : shift at current point
/// @param[in] shift_nxt : shift at next point
/// @param[in] dZ_loc : local distance increment
/// @param[in] increment : increment for the data vectors (either 1 or -1, depending on direction)
/// @param[in/out] id1 : index for the first data vector (everything except dZ)
/// @param[in/out] id2 : index for the second data vector (dZ); is stored from first_() to last_()-1
/// @note The interpolation is done assuming all lines are relevant for all frequencies. This
/// results in increased computation time
accel inline void Solver ::set_data(const Model& model, const Size crt, const Size nxt,
    const double shift_crt, const double shift_nxt, const double dZ_loc, const int increment,
    Size& id1, Size& id2) {
    Vector<double>& dZ            = dZ_();
    Vector<Size>& nr              = nr_();
    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();

    Vector<double>& shift = shift_();

    const double dshift     = shift_nxt - shift_crt;
    const double dshift_abs = fabs(dshift);

    // get number of interpolation points
    Size n_interp_points = interp_helper.get_n_interp(model, crt, nxt);

    if (n_interp_points > 1) {
        const double dZ_interpl     = dZ_loc / n_interp_points;
        const double dshift_interpl = dshift / n_interp_points;
        for (Size i = 1; i <= n_interp_points; i++) {
            nr[id1]            = nxt;
            shift[id1]         = shift_crt + i * dshift_interpl;
            nr_interp[id1]     = crt;
            interp_factor[id1] = 1.0 - ((double)i) / n_interp_points;

            dZ[id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }
    } else {
        nr[id1]            = nxt;
        shift[id1]         = shift_nxt;
        nr_interp[id1]     = crt;
        interp_factor[id1] = 0.0;
        dZ[id2]            = dZ_loc;

        id1 += increment;
        id2 += increment;
    }
}

/// Sets data for a ray, including required interpolation, assuming only a single line is important
/// for the computations.
/// @param[in] model : model object
/// @param[in] l : line index
/// @param[in] crt : current point
/// @param[in] nxt : next point
/// @param[in] shift_crt : shift at current point
/// @param[in] shift_nxt : shift at next point
/// @param[in] dZ_loc : local distance increment
/// @param[in] increment : increment for the data vectors (either 1 or -1, depending on direction)
/// @param[in/out] id1 : index for the first data vector (everything except dZ)
/// @param[in/out] id2 : index for the second data vector (dZ); is stored from first_() to last_()-1
accel inline void Solver ::set_data_for_line(const Model& model, const Size l, const Size crt,
    const Size nxt, const double shift_crt, const double shift_nxt, const double dZ_loc,
    const int increment, Size& id1, Size& id2) {
    Vector<double>& dZ            = dZ_();
    Vector<Size>& nr              = nr_();
    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();

    Vector<double>& shift = shift_();

    const double dshift     = shift_nxt - shift_crt;
    const double dshift_abs = fabs(dshift);

    // get number of interpolation points
    Size n_interp_points = interp_helper.get_n_interp_for_line(model, l, crt, nxt);

    if (n_interp_points > 1) {
        const double dZ_interpl     = dZ_loc / n_interp_points;
        const double dshift_interpl = dshift / n_interp_points;
        for (Size i = 1; i <= n_interp_points; i++) {
            nr[id1]            = nxt;
            shift[id1]         = shift_crt + i * dshift_interpl;
            nr_interp[id1]     = crt;
            interp_factor[id1] = 1.0 - ((double)i) / n_interp_points;

            dZ[id2] = dZ_interpl;

            id1 += increment;
            id2 += increment;
        }
    } else {
        nr[id1]            = nxt;
        shift[id1]         = shift_nxt;
        nr_interp[id1]     = crt;
        interp_factor[id1] = 0.0;
        dZ[id2]            = dZ_loc;

        id1 += increment;
        id2 += increment;
    }
}

///  Gaussian line profile function
///    @param[in] width : profile width
///    @param[in] diff  : frequency difference with line
///    centre
///    @return profile function evaluated with this
///    frequency difference
////////////////////////////////////////////////////////////////////////
accel inline Real Solver ::gaussian(const Real inverse_width, const Real diff) const {
    const Real sqrt_exp = inverse_width * diff;

    return inverse_width * INVERSE_SQRT_PI * expf(-sqrt_exp * sqrt_exp);
}

///  Planck function
///    @param[in] temp : temperature of the corresponding
///    black body
///    @param[in] freq : frequency at which to evaluate the
///    function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
accel inline Real Solver ::planck(const Real temp, const Real freq) const {
    return TWO_HH_OVER_CC_SQUARED * (freq * freq * freq) / expm1f(HH_OVER_KB * freq / temp);
}

///  Getter for the boundary conditions
///    @param[in] model  : reference to model object
///    @param[in] p      : point index of the boundary point
///    @param[in] freq   : frequency at which to evaluate
///    boundary condition
///    @returns incoming radiation intensity at the boundary
////////////////////////////////////////////////////////////////////////////
accel inline Real Solver ::boundary_intensity(
    const Model& model, const Size p, const Real freq) const {
    const Size bdy_id = model.geometry.boundary.point2boundary[p];

    switch (model.geometry.boundary.boundary_condition[bdy_id]) {
    case Zero:
        return 0.0;
    case Thermal:
        return planck(model.geometry.boundary.boundary_temperature[bdy_id], freq);
    default:
        return planck(T_CMB, freq);
    }
}

///  Getter for the emissivity (eta) and the opacity (chi)
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
template <>
accel inline void Solver ::get_eta_and_chi<None>(const Model& model, const Size p,
    const Size ll, // dummy variable
    const Real freq, Real& eta, Real& chi) const {
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters->nlines(); l++) {
        const Real diff = freq - model.lines.line[l];
        const Real prof = gaussian(model.lines.inverse_width(p, l), diff);

        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity(p, l);
    }
}

/// Getter for the emissivity (eta) and the opacity (chi), with interpolation support
/// @param[in] model : reference to model object
/// @param[in] p1 : index of the first point
/// @param[in] p2 : index of the second point
/// @param[in] interp_factor : interpolation factor between both points (0->p1, 1->p2)
/// @param[in] l : line index corresponding to the frequency
/// @param[in] freq : frequency (in co-moving frame)
/// @param[out] eta : emissivity
/// @param[out] chi : opacity
template <>
accel inline void Solver ::get_eta_and_chi_interpolated<None>(const Model& model, const Size p1,
    const Size p2, const Real interp_factor,
    const Size ll, // dummy variable
    const Real freq, Real& eta, Real& chi) const {
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters->nlines(); l++) {
        const Real diff             = freq - model.lines.line[l];
        const Real interp_inv_width = interp_helper.interpolate_linear(
            model.lines.inverse_width(p1, l), model.lines.inverse_width(p2, l), interp_factor);
        const Real prof = gaussian(interp_inv_width, diff);

        const Real interp_emissivity = interp_helper.interpolate_log(
            model.lines.emissivity(p1, l), model.lines.emissivity(p2, l), interp_factor);
        const Real interp_opacity = interp_helper.interpolate_log(
            model.lines.opacity(p1, l), model.lines.opacity(p2, l), interp_factor);

        eta += prof * interp_emissivity;
        chi += prof * interp_opacity;
    }
}

///  Getter for the emissivity (eta) and the opacity (chi)
///  function uses only the nearby lines to save computation
///  time
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
template <>
accel inline void Solver ::get_eta_and_chi<CloseLines>(const Model& model, const Size p,
    const Size ll, // dummy variable
    const Real freq, Real& eta, Real& chi) const {
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    const Real upper_bound_line_width =
        model.parameters->max_distance_opacity_contribution
        * model.thermodynamics.profile_width_upper_bound_with_linefreq(
            p, freq, model.lines.max_inverse_mass);
    const Real left_freq_bound  = freq - upper_bound_line_width;
    const Real right_freq_bound = freq + upper_bound_line_width;

    // Just using default search algorithms, obtaining
    // iterators
    auto left_line_bound = std::lower_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound = std::upper_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);

    for (auto freq_sort_l = left_line_bound; freq_sort_l != right_line_bound; freq_sort_l++) {
        const Size sort_l = freq_sort_l - model.lines.sorted_line.begin();
        // mapping sorted line index to original line index
        const Size l = model.lines.sorted_line_map[sort_l];
        // const Real diff = freq - model.lines.line[l];
        const Real diff = freq - *freq_sort_l; // should be equal to the
                                               // previous line of code
        const Real inv_width = model.lines.inverse_width(p, l);
        const Real prof      = gaussian(model.lines.inverse_width(p, l), diff);
        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity(p, l);
    }
}

/// Getter for the emissivity (eta) and the opacity (chi), with interpolation support.
/// This function uses only the nearby lines to save computation time
/// @param[in] model : reference to model object
/// @param[in] p1 : index of the first point
/// @param[in] p2 : index of the second point
/// @param[in] interp_factor : interpolation factor between both points (0->p1, 1->p2)
/// @param[in] l : line index corresponding to the frequency
/// @param[in] freq : frequency (in co-moving frame)
/// @param[out] eta : emissivity
/// @param[out] chi : opacity
template <>
accel inline void Solver ::get_eta_and_chi_interpolated<CloseLines>(const Model& model,
    const Size p1, const Size p2, const Real interp_factor,
    const Size ll, // dummy variable
    const Real freq, Real& eta, Real& chi) const {
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    const Real upper_bound_line_width =
        model.parameters->max_distance_opacity_contribution
        * model.thermodynamics.profile_width_upper_bound_with_linefreq(p1, freq,
            model.lines.max_inverse_mass); // note: assumes that the profile width does
                                           // not significantly change between p1 and p2
    const Real left_freq_bound  = freq - upper_bound_line_width;
    const Real right_freq_bound = freq + upper_bound_line_width;

    // Just using default search algorithms, obtaining
    // iterators
    auto left_line_bound = std::lower_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound = std::upper_bound(
        model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);

    for (auto freq_sort_l = left_line_bound; freq_sort_l != right_line_bound; freq_sort_l++) {
        const Size sort_l = freq_sort_l - model.lines.sorted_line.begin();
        // mapping sorted line index to original line index
        const Size l = model.lines.sorted_line_map[sort_l];
        // const Real diff = freq - model.lines.line[l];
        const Real diff = freq - *freq_sort_l; // should be equal to the
                                               // previous line of code
        const Real interp_inv_width = interp_helper.interpolate_linear(
            model.lines.inverse_width(p1, l), model.lines.inverse_width(p2, l), interp_factor);
        const Real prof = gaussian(interp_inv_width, diff);

        const Real interp_emissivity = interp_helper.interpolate_log(
            model.lines.emissivity(p1, l), model.lines.emissivity(p2, l), interp_factor);
        const Real interp_opacity = interp_helper.interpolate_log(
            model.lines.opacity(p1, l), model.lines.opacity(p2, l), interp_factor);

        eta += prof * interp_emissivity;
        chi += prof * interp_opacity;
    }
}

///  Getter for the emissivity (eta) and the opacity (chi)
///  in the "one line" approixmation
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  l     : line index corresponding to the
///    frequency
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////////////////////////////////////
template <>
accel inline void Solver ::get_eta_and_chi<OneLine>(
    const Model& model, const Size p, const Size l, const Real freq, Real& eta, Real& chi) const {
    const Real diff = freq - model.lines.line[l];
    const Real prof = gaussian(model.lines.inverse_width(p, l), diff);

    eta = prof * model.lines.emissivity(p, l);
    chi = prof * model.lines.opacity(p, l) + model.parameters->min_opacity;
}

/// Getter for the emissivity (eta) and the opacity (chi), with interpolation support
/// This functions uses only a single line for the computations
/// @param[in] model : reference to model object
/// @param[in] p1 : index of the first point
/// @param[in] p2 : index of the second point
/// @param[in] interp_factor : interpolation factor between both points (0->p1, 1->p2)
/// @param[in] l : line index corresponding to the frequency
/// @param[in] freq : frequency (in co-moving frame)
/// @param[out] eta : emissivity
/// @param[out] chi : opacity
template <>
accel inline void Solver ::get_eta_and_chi_interpolated<OneLine>(const Model& model, const Size p1,
    const Size p2, const Real interp_factor, const Size l, const Real freq, Real& eta,
    Real& chi) const {

    const Real diff             = freq - model.lines.line[l];
    const Real interp_inv_width = interp_helper.interpolate_linear(
        model.lines.inverse_width(p1, l), model.lines.inverse_width(p2, l), interp_factor);
    const Real prof = gaussian(interp_inv_width, diff);

    const Real interp_emissivity = interp_helper.interpolate_log(
        model.lines.emissivity(p1, l), model.lines.emissivity(p2, l), interp_factor);
    const Real interp_opacity = interp_helper.interpolate_log(
        model.lines.opacity(p1, l), model.lines.opacity(p2, l), interp_factor);

    eta = prof * interp_emissivity;
    chi = prof * interp_opacity + model.parameters->min_opacity;
}

///  Apply trapezium rule to x_crt and x_nxt
///    @param[in] x_crt : current value of x
///    @param[in] x_nxt : next value of x
///    @param[in] dZ    : distance inscrement along ray
///    @returns integral x over dZ
///////////////////////////////////////////////////////
accel inline Real trap(const Real x_crt, const Real x_nxt, const double dZ) {
    return half * (x_crt + x_nxt) * dZ;
}

// Approach for shortchar solver which included raytracing in the solver itself; was incompatible
// with the interpolation implementation
// template <ApproximationType approx, bool use_adaptive_directions>
// accel inline void Solver ::solve_shortchar_order_0(Model& model,
// const Size o, const Size r) {
//     Vector<Real>& eta_c = eta_c_();
//     Vector<Real>& eta_n = eta_n_();
//
//     Vector<Real>& chi_c = chi_c_();
//     Vector<Real>& chi_n = chi_n_();
//
//     Vector<Real>& source_c = source_c_();
//     Vector<Real>& source_n = source_n_();
//
//     Vector<Real>& tau = tau_();
//
//     double Z  = 0.0; // distance along ray
//     double dZ = 0.0; // last distance increment
//
//     Size crt = o;
//     Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);
//     Real term_c, term_n, dtau;
//     bool compute_curr_opacity, prev_compute_curr_opacity;
//     const Real ray_weight = model.geometry.rays.get_weight<use_adaptive_directions>(o,
//     r);
//
//     if (model.geometry.valid_point(nxt)) {
//         double shift_c = 1.0;
//         double shift_n = model.geometry.get_shift<CoMoving, use_adaptive_directions>(o,
//         r, nxt, Z);
//
//         for (Size f = 0; f < model.parameters->nfreqs(); f++) {
//             const Real freq = model.radiation.frequencies.nu(o, f);
//             const Size l    = model.radiation.frequencies.corresponding_line[f]; // line
//             index
//
//             compute_curr_opacity = true; // for the first point, we need to compute
//                                          // both the curr and next opacity (and source)
//
//             // TODO: refactor such that this also uses set_data()
//             Size crt_interp_idx = crt;
//             Size nxt_interp_idx = nxt;
//             Real curr_interp    = 1.0;
//             Real next_interp    = 1.0;
//
//             compute_source_dtau<approx>(model, crt, nxt, crt_interp_idx, nxt_interp_idx,
//             l,
//                 freq * shift_c, freq * shift_n, shift_c, shift_n, dZ, curr_interp,
//                 next_interp, compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f],
//                 source_n[f]);
//             dtau = std::max(model.parameters->min_dtau, dtau);
//
//             // proper implementation of 2nd order shortchar (not
//             // yet times reducing factor of exp(-tau))
//             //  model.radiation.I(r,o,f) = term_c *
//             //  (expm1(-dtau)+dtau) / dtau
//             //                           + term_n *
//             //                           (-expm1(-dtau)-dtau*expf(-dtau))
//             //                           /dtau;
//             // Rewrite, trying to use less exponentials
//             const Real factor = expm1f(-dtau) / dtau;
//
//             model.radiation.I(r, o, f) =
//                 factor * (source_c[f] - source_n[f] * (1.0 + dtau)) + source_c[f] -
//                 source_n[f];
//             tau[f] = dtau;
//
//             // Compute local lambda operator
//             const Size l_spec =
//                 model.radiation.frequencies.corresponding_l_for_spec[f]; // index of
//                 species
//             const Size k = model.radiation.frequencies.corresponding_k_for_tran[f]; //
//             index of
//                                                                                     //
//                                                                                     transition
//             const Size z =
//                 model.radiation.frequencies.corresponding_z_for_line[f]; // index of
//                                                                          // quadrature
//                                                                          point
//             const Real w_ang = ray_weight;
//
//             LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l_spec];
//
//             const Real freq_line = lspec.linedata.frequency[k];
//             const Real invr_mass = lspec.linedata.inverse_mass;
//             const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] *
//             w_ang;
//
//             Real eta, chi; // eta is dummy var
//             // chi is not necessarily computed, so compute it to
//             // be sure
//             get_eta_and_chi<approx>(model, o, k, freq_line, eta, chi);
//             Real inverse_chi = 1.0 / chi;
//             Real phi         = model.thermodynamics.profile(invr_mass, o, freq_line,
//             freq);
//             // const Real lambda_factor =
//             // (dtau+expm1f(-dtau))/dtau;// If one wants to
//             // compute lambda a bit more accurately in case of
//             // dtau≃0. Real L   = constante * freq * phi *
//             // lambda_factor * inverse_chi;
//             Real L = constante * freq * phi * (factor + 1.0)
//                    * inverse_chi; // using factor+1.0, the
//                                   // computed lambda elements can
//                                   // be negative if dtau very
//                                   // small; but then the lambda
//                                   // elements are also negligible
//             lspec.lambda.add_element(o, k, o, L);
//
//             // TODO: possible nonlocal lambda part // FIXME:
//             // probably incorrect chi used
//             //  L   = constante * freq * phi * (-factor *
//             //  (1.0+dtau) - 1.0) * inverse_chi;
//             //  lspec.lambda.add_element(o, k, nxt, L);
//         }
//
//         // For all frequencies, we need to use the same method
//         // for computing the optical depth
//         //  bool
//         //  prev_compute_curr_opacity=compute_curr_opacity;//technically,
//         //  we could also keep this bool individually for every
//         //  frequency
//         prev_compute_curr_opacity = compute_curr_opacity; // technically, we could also
//                                                           // keep this bool individually
//                                                           // for every frequency
//
//         while (model.geometry.not_on_boundary(nxt)) {
//             crt     = nxt;
//             shift_c = shift_n;
//
//             model.geometry.get_next<use_adaptive_directions>(o, r, crt, nxt, Z, dZ,
//             shift_n);
//
//             for (Size f = 0; f < model.parameters->nfreqs(); f++) {
//                 source_c[f]     = source_n[f];
//                 chi_c[f]        = chi_n[f];
//                 const Real freq = model.radiation.frequencies.nu(o, f);
//                 const Size l    = model.radiation.frequencies.corresponding_line[f];
//
//                 compute_curr_opacity = prev_compute_curr_opacity;
//
//                 Size crt_interp_idx = crt;
//                 Size nxt_interp_idx = nxt;
//                 Real curr_interp    = 1.0;
//                 Real next_interp    = 1.0;
//
//                 compute_source_dtau<approx>(model, crt, nxt, crt_interp_idx,
//                 nxt_interp_idx, l,
//                     freq * shift_c, freq * shift_n, shift_c, shift_n, dZ, curr_interp,
//                     next_interp, compute_curr_opacity, dtau, chi_c[f], chi_n[f],
//                     source_c[f], source_n[f]);
//                 dtau = std::max(model.parameters->min_dtau, dtau);

//                 // proper implementation of 2nd order shortchar (not
//                 // yet times reducing factor of exp(-tau))
//                 //  model.radiation.I(r,o,f) += expf(-tau[f]) *
//                 //                           ( term_c *
//                 //                           (expm1(-dtau)+dtau) /
//                 //                           dtau
//                 //                           + term_n *
//                 //                           (-expm1(-dtau)-dtau*expf(-dtau))
//                 //                           /dtau);
//                 // Rewrite, trying to use less exponentials
//                 model.radiation.I(r, o, f) +=
//                     expf(-tau[f])
//                     * (expm1f(-dtau) / dtau * (source_c[f] - source_n[f] * (1.0 + dtau))
//                         + source_c[f] - source_n[f]);
//                 // TODO: check order of addition, as we might be
//                 // starting with the largest contributions, before
//                 // adding the smaller ones...
//                 tau[f] += dtau;
//             }
//
//             // save setting for use for all frequencies for the
//             // next interval
//             prev_compute_curr_opacity = compute_curr_opacity;
//         }
//
//         for (Size f = 0; f < model.parameters->nfreqs(); f++) {
//             const Real freq = model.radiation.frequencies.nu(o, f);
//
//             model.radiation.I(r, o, f) +=
//                 boundary_intensity(model, nxt, freq * shift_n) * expf(-tau[f]);
//             model.radiation.J(o, f) += ray_weight * model.radiation.I(r, o, f);
//         }
//     }
//
//     else {
//         for (Size f = 0; f < model.parameters->nfreqs(); f++) {
//             const Real freq = model.radiation.frequencies.nu(o, f);
//
//             model.radiation.I(r, o, f) = boundary_intensity(model, crt, freq);
//             model.radiation.J(o, f) += ray_weight * model.radiation.I(r, o, f);
//         }
//     }
// }

/// Solve the radiative transfer equation using the short characteristics method
/// @param model : reference to model object
/// @param o : origin point index
/// @param r : ray index (hnrays)
/// @param f : frequency index
template <ApproximationType approx, bool use_adaptive_directions>
accel inline Real Solver ::solve_shortchar_order_0_ray_forward(
    Model& model, const Size o, const Size r, const Size f) {
    Real I = 0.0;
    Real eta_c, eta_n, chi_c, chi_n, source_c, source_n;

    Vector<Real>& tau     = tau_();
    Vector<double>& shift = shift_();
    Vector<Size>& nr      = nr_();
    Vector<double>& dZ    = dZ_();

    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();

    Size crt, nxt, crt_interp_idx, nxt_interp_idx;

    Real term_c, term_n, dtau, crt_interp, nxt_interp;
    bool compute_curr_opacity = true;
    const Real ray_weight     = model.geometry.rays.get_weight<use_adaptive_directions>(o, r);

    const Real freq = model.radiation.frequencies.nu(o, f);

    I = boundary_intensity(model, nr[first_()], freq * shift[first_()]);
    double shift_c, shift_n;
    // interate until we reach the middle point
    for (Size idx = first_() + 1; idx <= centre; idx++) {
        crt            = nr[idx - 1];
        nxt            = nr[idx];
        crt_interp_idx = nr_interp[idx - 1];
        nxt_interp_idx = nr_interp[idx];
        crt_interp     = interp_factor[idx - 1];
        nxt_interp     = interp_factor[idx];
        shift_c        = shift[idx - 1];
        shift_n        = shift[idx];

        chi_c    = chi_n;
        source_c = source_n;

        const Real freq = model.radiation.frequencies.nu(o, f);
        const Size l    = model.radiation.frequencies.corresponding_line[f];

        compute_source_dtau<approx>(model, crt, nxt, crt_interp_idx, nxt_interp_idx, l,
            freq * shift_c, freq * shift_n, shift_c, shift_n, dZ[idx - 1], crt_interp, nxt_interp,
            compute_curr_opacity, dtau, chi_c, chi_n, source_c, source_n);
        dtau = std::max(model.parameters->min_dtau, dtau);

        // model.radiation.I(r, o, f) = exp(-dtau) * model.radiation.I(r, o, f)
        //                            + term_c * (-expm1(-dtau) - dtau * expf(-dtau)) / dtau
        //                            + term_n * (expm1(-dtau) + dtau) / dtau;
        // slight rewrite, should be more stable
        // I = expf(-dtau) * I
        //   + (source_c * (-expm1f(-dtau) - dtau * expf(-dtau)) + source_n * (expm1f(-dtau) +
        //   dtau)) / dtau;
        I = expf(-dtau) * I + expm1f(-dtau) / dtau * (source_c - source_n * (1.0 + dtau)) + source_c
          - source_n;
    }
    // Compute local lambda operator
    const Real factor = expm1f(-dtau) / dtau;

    const Size l_spec = model.radiation.frequencies.corresponding_l_for_spec[f]; // index of species
    const Size k      = model.radiation.frequencies.corresponding_k_for_tran[f]; // index of
                                                                                 // transition
    const Size z = model.radiation.frequencies.corresponding_z_for_line[f];      // index of
                                                                                 // quadrature point
    const Real w_ang = ray_weight;

    LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l_spec];

    const Real freq_line = lspec.linedata.frequency[k];
    const Real invr_mass = lspec.linedata.inverse_mass;
    const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;

    Real eta, chi; // eta is dummy var
    // chi needs to be recomputed, as we have no guarantee the for loop runs at least once
    get_eta_and_chi<approx>(model, o, k, freq_line, eta, chi);
    Real inverse_chi = 1.0 / chi;
    Real phi         = model.thermodynamics.profile(invr_mass, o, freq_line, freq);
    // const Real lambda_factor =
    // (dtau+expm1f(-dtau))/dtau;// If one wants to
    // compute lambda a bit more accurately in case of
    // dtau≃0. Real L   = constante * freq * phi *
    // lambda_factor * inverse_chi;
    Real L = constante * freq * phi * (factor + 1.0) * inverse_chi; // using factor+1.0, the
                                                                    // computed lambda elements can
                                                                    // be negative if dtau very
                                                                    // small; but then the lambda
                                                                    // elements are also negligible
    lspec.lambda.add_element(o, k, o, L);

    return I;
}

/// Computes the intensity using the short characteristics method, for the backward direction
/// @param model : reference to model object
/// @param o : origin point index
/// @param r : ray index (hnrays)
/// @param f : frequency index
// TODO: reorder computations to start from the centre; in this way, we can avoid the constant I =
// I*exp(-dtau) + ..., which gives more numerical rounding errors; (also do this for the forward
// ray)
template <ApproximationType approx, bool use_adaptive_directions>
accel inline Real Solver ::solve_shortchar_order_0_ray_backward(
    Model& model, const Size o, const Size r, const Size f) {
    Real I = 0.0;
    Real eta_c, eta_n, chi_c, chi_n, source_c, source_n;

    Vector<Real>& tau     = tau_();
    Vector<double>& shift = shift_();
    Vector<Size>& nr      = nr_();
    Vector<double>& dZ    = dZ_();

    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();

    Size crt, nxt, crt_interp_idx, nxt_interp_idx;

    Real term_c, term_n, dtau, crt_interp, nxt_interp;
    bool compute_curr_opacity = true;
    const Real ray_weight     = model.geometry.rays.get_weight<use_adaptive_directions>(o, r);

    // Set boundary condition
    const Real freq = model.radiation.frequencies.nu(o, f);

    I = boundary_intensity(model, nr[last_()], freq * (shift[last_()]));

    double shift_c, shift_n;
    // interate until we reach the middle point
    for (Size indexp1 = last_(); indexp1 >= centre + 1; indexp1--) {
        crt            = nr[indexp1];
        nxt            = nr[indexp1 - 1];
        crt_interp_idx = nr_interp[indexp1];
        nxt_interp_idx = nr_interp[indexp1 - 1];
        crt_interp     = interp_factor[indexp1];
        nxt_interp     = interp_factor[indexp1 - 1];
        shift_c        = shift[indexp1];
        shift_n        = shift[indexp1 - 1];

        chi_c    = chi_n;
        source_c = source_n;

        const Real freq = model.radiation.frequencies.nu(o, f);
        const Size l    = model.radiation.frequencies.corresponding_line[f];

        compute_source_dtau<approx>(model, crt, nxt, crt_interp_idx, nxt_interp_idx, l,
            freq * shift_c, freq * shift_n, shift_c, shift_n, dZ[indexp1 - 1], crt_interp,
            nxt_interp, compute_curr_opacity, dtau, chi_c, chi_n, source_c, source_n);
        dtau = std::max(model.parameters->min_dtau, dtau);

        // model.radiation.I(r,o,f) =
        // exp(-dtau) * model.radiation.I(r, o, f)
        //  + term_c * (-expm1(-dtau) - dtau * expf(-dtau)) / dtau
        //  + term_n * (expm1(-dtau) + dtau) / dtau;
        // slight rewrite, should be more stable
        // I = expf(-dtau) * I
        //   + (source_c * (-expm1f(-dtau) - dtau * expf(-dtau)) + source_n * (expm1f(-dtau) +
        //   dtau))/ dtau;
        I = expf(-dtau) * I + expm1f(-dtau) / dtau * (source_c - source_n * (1.0 + dtau)) + source_c
          - source_n;
    }
    // Compute local lambda operator
    const Real factor = expm1f(-dtau) / dtau;

    const Size l_spec = model.radiation.frequencies.corresponding_l_for_spec[f]; // index of species
    const Size k      = model.radiation.frequencies.corresponding_k_for_tran[f]; // index of
                                                                                 // transition
    const Size z = model.radiation.frequencies.corresponding_z_for_line[f];      // index of
                                                                                 // quadrature point
    const Real w_ang = ray_weight;

    LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l_spec];

    const Real freq_line = lspec.linedata.frequency[k];
    const Real invr_mass = lspec.linedata.inverse_mass;
    const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;

    Real eta, chi; // eta is dummy var
    // chi has to be recomputed, as we have no guarantee the for loop runs at least once
    get_eta_and_chi<approx>(model, o, k, freq_line, eta, chi);
    Real inverse_chi = 1.0 / chi;
    Real phi         = model.thermodynamics.profile(invr_mass, o, freq_line, freq);
    // const Real lambda_factor =
    // (dtau+expm1f(-dtau))/dtau;// If one wants to
    // compute lambda a bit more accurately in case of
    // dtau≃0. Real L   = constante * freq * phi *
    // lambda_factor * inverse_chi;
    Real L = constante * freq * phi * (factor + 1.0) * inverse_chi; // using factor+1.0, the
                                                                    // computed lambda elements can
                                                                    // be negative if dtau very
                                                                    // small; but then the lambda
                                                                    // elements are also negligible
    lspec.lambda.add_element(o, k, o, L);

    return I;
}

// only for the feautrier solver
template <ApproximationType approx, bool use_adaptive_directions>
accel inline void Solver ::update_Lambda(Model& model, const Size rr, const Size f) {
    const Frequencies& freqs        = model.radiation.frequencies;
    const Thermodynamics& thermodyn = model.thermodynamics;

    if (freqs.appears_in_line_integral[f]) {
        const Size first = first_();
        const Size last  = last_();
        const Size n_tot = n_tot_();

        Vector<Size>& nr      = nr_();
        Vector<double>& shift = shift_();
        Vector<Real>& L_diag  = L_diag_();
        Matrix<Real>& L_upper = L_upper_();
        Matrix<Real>& L_lower = L_lower_();
        // Vector<Real  >& inverse_chi = inverse_chi_();

        const Real w_ang =
            two * model.geometry.rays.get_weight<use_adaptive_directions>(nr[centre], rr);

        const Size l = freqs.corresponding_l_for_spec[f]; // index of species
        const Size k = freqs.corresponding_k_for_tran[f]; // index of transition
        const Size z = freqs.corresponding_z_for_line[f]; // index of
                                                          // quadrature point

        LineProducingSpecies& lspec = model.lines.lineProducingSpecies[l];

        const Real freq_line = lspec.linedata.frequency[k];
        const Real invr_mass = lspec.linedata.inverse_mass;
        const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;
        // TODO: approximation if we do not have overlapping
        // lines
        //  Real inverse_opacity = 1.0/model.lines.opacity
        //  (nr[centre], k);
        //  //Inverse line opacity; includes 1/HH_OVER_FOUR_PI
        //  Real L   = constante
        //  * L_diag[centre] * inverse_opacity;
        Real eta, chi; // eta is dummy var
        get_eta_and_chi<approx>(model, nr[centre], k, freq_line, eta, chi);
        Real inverse_chi = 1.0 / chi;

        Real frq = freqs.nu(nr[centre], f) * shift[centre];
        Real phi = thermodyn.profile(invr_mass, nr[centre], freq_line, frq);
        Real L   = constante * frq * phi * L_diag[centre] * inverse_chi;

        lspec.lambda.add_element(nr[centre], k, nr[centre], L);

        for (long m = 0; (m < n_off_diag) && (m + 1 < n_tot); m++) {
            if (centre >= first + m + 1) // centre-m-1 >= first
            {
                const long n = centre - m - 1;

                // TODO: approximation if we do not have overlapping
                // lines // also check which opacity we need
                // (centre?)
                //  inverse_opacity = 1.0/model.lines.opacity
                //  (nr[n], k); //Inverse line opacity; no longer
                //  1/HH_OVER_FOUR_PI L   = constante * L_lower(m,n)
                //  * inverse_opacity;
                get_eta_and_chi<approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi = 1.0 / chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile(invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_lower(m, n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }

            if (centre + m + 1 <= last) // centre+m+1 < last
            {
                const long n = centre + m + 1;

                // TODO: approximation if we do not have overlapping
                // lines
                //  inverse_opacity = 1.0/model.lines.opacity
                //  (nr[n], k); //Inverse line opacity; no includes
                //  1/HH_OVER_FOUR_PI L   = constante * L_upper(m,n)
                //  * inverse_opacity;
                get_eta_and_chi<approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi = 1.0 / chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile(invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_upper(m, n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }
        }
    }
}

///  Computer for the optical depth and source function when
///  computing using the formal line integration In case of
///  low velocity differences, almost two times slower For
///  high velocity increments however, this does not need
///  any extra interpolation points (-> way faster) (and
///  some extra optimizations are made in the limit of
///  extremely large doppler shift (as erf goes to its limit
///  value)
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] curr_point_interp_idx : index of interpolation point for curr_point
///    @param[in] next_point_interp_idx : index of interpolation point for next_point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in
///    comoving frame)
///    @param[in] nextfreq : frequency at next point (in
///    comoving frame)
///    @param[in] dZ : position increment
///    @param[in] curr_interp : interpolation factor for curr_point
///    @param[in] next_interp : interpolation factor for next_point
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point
///    to compute
///    @param[out] Snext : source function at next point to
///    compute
/////////////////////////////////////////////////////////////////////
template <>
inline void Solver ::compute_S_dtau_line_integrated<OneLine>(Model& model, Size currpoint,
    Size nextpoint, Size currpoint_interp_idx, Size nextpoint_interp_idx, Size lineidx,
    Real currfreq, Real nextfreq, Real curr_interp, Real next_interp, Real dZ, Real& dtau,
    Real& Scurr, Real& Snext) {
    // FIXME: line idx is wrong
    dtau = compute_dtau_single_line(model, currpoint, nextpoint, currpoint_interp_idx,
        nextpoint_interp_idx, lineidx, currfreq, nextfreq, curr_interp, next_interp, dZ);
    // const Real curr_opacity = interp_helper.interpolate_log(model.lines.opacity(currpoint,
    // lineidx),
    //     model.lines.opacity(currpoint_interp_idx, lineidx), curr_interp);
    // const Real next_opacity = interp_helper.interpolate_log(model.lines.opacity(nextpoint,
    // lineidx),
    //     model.lines.opacity(nextpoint_interp_idx, lineidx), next_interp);
    // const Real curr_emissivity =
    //     interp_helper.interpolate_log(model.lines.emissivity(currpoint, lineidx),
    //         model.lines.emissivity(currpoint_interp_idx, lineidx), curr_interp);
    // const Real next_emissivity =
    //     interp_helper.interpolate_log(model.lines.emissivity(nextpoint, lineidx),
    //         model.lines.emissivity(nextpoint_interp_idx, lineidx), next_interp);
    // Scurr = curr_emissivity / curr_opacity;
    // Snext = next_emissivity / next_opacity;

    Real Scurr_p =
        model.lines.emissivity(currpoint, lineidx) / model.lines.opacity(currpoint, lineidx);
    Real Snext_p =
        model.lines.emissivity(nextpoint, lineidx) / model.lines.opacity(nextpoint, lineidx);
    Real Scurr_interp_p = model.lines.emissivity(currpoint_interp_idx, lineidx)
                        / model.lines.opacity(currpoint_interp_idx, lineidx);
    Real Snext_interp_p = model.lines.emissivity(nextpoint_interp_idx, lineidx)
                        / model.lines.opacity(nextpoint_interp_idx, lineidx);

    Scurr = interp_helper.interpolate_log(Scurr_p, Scurr_interp_p, curr_interp);
    Snext = interp_helper.interpolate_log(Snext_p, Snext_interp_p, next_interp);
}

///  Computer for the optical depth and source function when
///  computing using the formal line integration In case of
///  low velocity differences, almost two times slower For
///  high velocity increments however, this does not need
///  any extra interpolation points (-> way faster) (and
///  some extra optimizations are made in the limit of
///  extremely large doppler shift (as erf goes to its limit
///  value)
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] curr_point_interp_idx : index of interpolation point for curr_point
///    @param[in] next_point_interp_idx : index of interpolation point for next_point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in
///    comoving frame)
///    @param[in] nextfreq : frequency at next point (in
///    comoving frame)
///    @param[in] dZ : position increment
///    @param[in] curr_interp : interpolation factor for curr_point
///    @param[in] next_interp : interpolation factor for next_point
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point
///    to compute
///    @param[out] Snext : source function at next point to
///    compute
/////////////////////////////////////////////////////////////////////
template <>
inline void Solver ::compute_S_dtau_line_integrated<None>(Model& model, Size currpoint,
    Size nextpoint, Size currpoint_interp_idx, Size nextpoint_interp_idx, Size lineidx,
    Real currfreq, Real nextfreq, Real curr_interp, Real next_interp, Real dZ, Real& dtau,
    Real& Scurr, Real& Snext) {
    Real sum_dtau             = 0.0;
    Real sum_dtau_times_Scurr = 0.0;
    Real sum_dtau_times_Snext = 0.0;
    for (Size l = 0; l < model.parameters->nlines(); l++) {
        Real line_dtau = compute_dtau_single_line(model, currpoint, nextpoint, currpoint_interp_idx,
            nextpoint_interp_idx, l, currfreq, nextfreq, curr_interp, next_interp, dZ);
        // const Real curr_opacity =
        //     interp_helper.interpolate_log(model.lines.opacity(currpoint, lineidx),
        //         model.lines.opacity(currpoint_interp_idx, lineidx), curr_interp);
        // const Real next_opacity =
        //     interp_helper.interpolate_log(model.lines.opacity(nextpoint, lineidx),
        //         model.lines.opacity(nextpoint_interp_idx, lineidx), next_interp);
        // const Real curr_emissivity =
        //     interp_helper.interpolate_log(model.lines.emissivity(currpoint, lineidx),
        //         model.lines.emissivity(currpoint_interp_idx, lineidx), curr_interp);
        // const Real next_emissivity =
        //     interp_helper.interpolate_log(model.lines.emissivity(nextpoint, lineidx),
        //         model.lines.emissivity(nextpoint_interp_idx, lineidx), next_interp);
        // Real line_Scurr = curr_emissivity / curr_opacity;
        // Real line_Snext = next_emissivity / next_opacity;

        Real Scurr_p = model.lines.emissivity(currpoint, l) / model.lines.opacity(currpoint, l);
        Real Snext_p = model.lines.emissivity(nextpoint, l) / model.lines.opacity(nextpoint, l);
        Real Scurr_interp_p = model.lines.emissivity(currpoint_interp_idx, l)
                            / model.lines.opacity(currpoint_interp_idx, l);
        Real Snext_interp_p = model.lines.emissivity(nextpoint_interp_idx, l)
                            / model.lines.opacity(nextpoint_interp_idx, l);

        Real line_Scurr = interp_helper.interpolate_log(Scurr_p, Scurr_interp_p, curr_interp);
        Real line_Snext = interp_helper.interpolate_log(Snext_p, Snext_interp_p, next_interp);

        sum_dtau += line_dtau;
        sum_dtau_times_Scurr += line_dtau * line_Scurr;
        sum_dtau_times_Snext += line_dtau * line_Snext;
    }
    dtau  = sum_dtau;
    Scurr = sum_dtau_times_Scurr / sum_dtau;
    Snext = sum_dtau_times_Snext / sum_dtau;
}

///  Computer for the optical depth and source function when
///  computing using the formal line integration. In case of
///  low velocity differences, almost two times slower. For
///  high velocity increments however, this does not need
///  any extra interpolation points (-> way faster) (and
///  some extra optimizations are made in the limit of
///  extremely large doppler shift (as erf goes to its limit
///  value). This function only takes into account the nearby
///  lines, saving some computation time.
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] curr_point_interp_idx : index of interpolation point for curr_point
///    @param[in] next_point_interp_idx : index of interpolation point for next_point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in
///    comoving frame)
///    @param[in] nextfreq : frequency at next point (in
///    comoving frame)
///    @param[in] curr_interp : interpolation factor for curr_point
///    @param[in] next_interp : interpolation factor for next_point
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point
///    to compute
///    @param[out] Snext : source function at next point to
///    compute
/////////////////////////////////////////////////////////////////////
template <>
inline void Solver ::compute_S_dtau_line_integrated<CloseLines>(Model& model, Size currpoint,
    Size nextpoint, Size currpoint_interp_idx, Size nextpoint_interp_idx, Size lineidx,
    Real currfreq, Real nextfreq, Real curr_interp, Real next_interp, Real dZ, Real& dtau,
    Real& Scurr, Real& Snext) {
    Real sum_dtau = 0.0; // division by zero might occur otherwise
    // Real sum_dtau=model.parameters->min_dtau; //division by
    // zero might occur otherwise
    Real sum_dtau_times_Scurr = 0.0;
    Real sum_dtau_times_Snext = 0.0;

    Real left_freq;
    Real right_freq;

    // err, compiler will probably figure out that I just want
    // these two values ordered
    if (currfreq < nextfreq) {
        left_freq  = currfreq;
        right_freq = nextfreq;
    } else {
        right_freq = currfreq;
        left_freq  = nextfreq;
    }

    // using maximum of bounds on the two points to get an
    // upper bound for the line width
    const Real curr_bound_line_width = model.parameters->max_distance_opacity_contribution
                                     * model.thermodynamics.profile_width_upper_bound_with_linefreq(
                                         currpoint, right_freq, model.lines.max_inverse_mass);
    const Real next_bound_line_width = model.parameters->max_distance_opacity_contribution
                                     * model.thermodynamics.profile_width_upper_bound_with_linefreq(
                                         nextpoint, right_freq, model.lines.max_inverse_mass);
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

        Real line_dtau = compute_dtau_single_line(model, currpoint, nextpoint, currpoint_interp_idx,
            nextpoint_interp_idx, l, currfreq, nextfreq, curr_interp, next_interp, dZ);
        // const Real curr_opacity =
        //     interp_helper.interpolate_log(model.lines.opacity(currpoint, lineidx),
        //         model.lines.opacity(currpoint_interp_idx, lineidx), curr_interp);
        // const Real next_opacity =
        //     interp_helper.interpolate_log(model.lines.opacity(nextpoint, lineidx),
        //         model.lines.opacity(nextpoint_interp_idx, lineidx), next_interp);
        // const Real curr_emissivity =
        //     interp_helper.interpolate_log(model.lines.emissivity(currpoint, lineidx),
        //         model.lines.emissivity(currpoint_interp_idx, lineidx), curr_interp);
        // const Real next_emissivity =
        //     interp_helper.interpolate_log(model.lines.emissivity(nextpoint, lineidx),
        //         model.lines.emissivity(nextpoint_interp_idx, lineidx), next_interp);
        // Real line_Scurr = curr_emissivity / curr_opacity;
        // Real line_Snext = next_emissivity / next_opacity;

        Real Scurr_p = model.lines.emissivity(currpoint, l) / model.lines.opacity(currpoint, l);
        Real Snext_p = model.lines.emissivity(nextpoint, l) / model.lines.opacity(nextpoint, l);
        Real Scurr_interp_p = model.lines.emissivity(currpoint_interp_idx, l)
                            / model.lines.opacity(currpoint_interp_idx, l);
        Real Snext_interp_p = model.lines.emissivity(nextpoint_interp_idx, l)
                            / model.lines.opacity(nextpoint_interp_idx, l);

        Real line_Scurr = interp_helper.interpolate_log(Scurr_p, Scurr_interp_p, curr_interp);
        Real line_Snext = interp_helper.interpolate_log(Snext_p, Snext_interp_p, next_interp);

        sum_dtau += line_dtau;
        sum_dtau_times_Scurr += line_dtau * line_Scurr;
        sum_dtau_times_Snext += line_dtau * line_Snext;
    }
    dtau = sum_dtau;
    // needs extra bounding, as nothing may be added in the
    // first place (above for loop may have looped over 0
    // elements)
    const Real bound_min_dtau = model.parameters->min_opacity * dZ;
    // Correct way of bounding from below; should be able to
    // deal with very minor computation errors around 0.
    if (-bound_min_dtau < dtau) {
        dtau = std::max(bound_min_dtau, dtau);
    }
    // Note: 0 source functions can be returned if no lines
    // are nearby; but then the negligible lower bound gets
    // returned
    Scurr = sum_dtau_times_Scurr / dtau;
    Snext = sum_dtau_times_Snext / dtau;

    // note: due to interaction with dtau when computing all
    // sources individually, we do need to recompute Scurr and
    // Snext for all position increments
}

/// Computes the source function and optical depth in a
/// hybrid manner
///    @param[in/out] compute_curr_opacity: for deciding
///    whether we need to compute the current opacity when
///    using the trapezoidal rule
///    @param[in] currpoint : index of current point
///    @param[in] nextpoint : index of next point
///    @param[in] currpoint_interp_idx : index of interpolation point for current point
///    @param[in] nextpoint_interp_idx : index of interpolation point for next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in
///    comoving frame)
///    @param[in] nextfreq : frequency at next point (in
///    comoving frame)
///    @param[in] currshift : shift at curr point
///    @param[in] nextshift : shift at next point
///    @param[in] dZ : position increment
///    @param[in] curr_interp : interpolation fraction at current point
///    @param[in] next_interp : interpolation fraction at next point
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point
///    to compute
///    @param[out] Snext : source function at next point to
///    compute
/// Warning: depending on how large the doppler shift is,
/// the opacity is NOT computed.
template <ApproximationType approx>
accel inline void Solver ::compute_source_dtau(Model& model, Size currpoint, Size nextpoint,
    Size currpoint_interp_idx, Size nextpoint_interp_idx, Size line, Real curr_freq, Real next_freq,
    double curr_shift, double next_shift, Real dZ, double curr_interp, double next_interp,
    bool& compute_curr_opacity, Real& dtaunext, Real& chicurr, Real& chinext, Real& Scurr,
    Real& Snext) {
    // deciding which optical depth computation to use,
    // depending on the doppler shift
    const double dshift     = next_shift - curr_shift;
    const double dshift_abs = fabs(dshift);
    const double dshift_max = std::min(model.dshift_max[currpoint], model.dshift_max[nextpoint]);
    const bool using_large_shift = (dshift_abs > dshift_max);

    // fancy computation for large doppler shifts
    if (using_large_shift) {
        compute_curr_opacity = true;
        compute_S_dtau_line_integrated<approx>(model, currpoint, nextpoint, currpoint_interp_idx,
            nextpoint_interp_idx, line, curr_freq, next_freq, curr_interp, next_interp, dZ,
            dtaunext, Scurr, Snext);
        // OPACITY IS NOT COMPUTED IN THIS BRANCH!
    } else {
        // default computation using trapezoidal rule
        if (compute_curr_opacity) // fancy computation does not
                                  // compute the current
                                  // opacity, so we might need
                                  // to recompute it here
        {
            compute_curr_opacity = false;
            Real eta_c           = 0.0; // current emissivity
            // also get previous opacity (emissivity does not
            // matter)
            // get_eta_and_chi<approx>(model, currpoint, line, curr_freq, eta_c, chicurr);
            get_eta_and_chi_interpolated<approx>(model, currpoint, currpoint_interp_idx,
                curr_interp, line, curr_freq, eta_c, chicurr);
            Scurr = eta_c / chicurr; // might as well compute the
                                     // source function too
        }

        Real eta_n = 0.0;
        // Get new radiative properties
        get_eta_and_chi_interpolated<approx>(
            model, nextpoint, nextpoint_interp_idx, next_interp, line, next_freq, eta_n, chinext);

        Snext = eta_n / chinext;

        dtaunext = half * (chicurr + chinext) * dZ;
    }
}

///  Solver for Feautrier equation along ray pairs using the
///  (ordinary) 2nd-order solver, without adaptive optical
///  depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::solve_feautrier_order_2(Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ            = dZ_();
    Vector<Size>& nr              = nr_();
    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();
    Vector<double>& shift         = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], nr_interp[first],
        nr_interp[first + 1], l, freq * shift[first], freq * shift[first + 1], shift[first],
        shift[first + 1], dZ[first], interp_factor[first], interp_factor[first + 1],
        compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first]         = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = 1.0 / C[first]; // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) /
    /// C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], nr_interp[n], nr_interp[n + 1], l,
            freq * shift[n], freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], interp_factor[n],
            interp_factor[n + 1], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n]        = dtau_avg * dtau_c;
        inverse_C[n]        = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source
        /// function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;

    if (n_off_diag == 0) {
        if (centre < last) {
            /// Write economically: G[last] = (B[last] - A[last])
            /// / A[last];
            GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
            GP[last] = GG[last] / (one + GG[last]);

            for (long n = last - 1; n > centre; n--) // use long in reverse loops!
            {
                Su[n] += Su[n + 1] * FI[n];

                GG[n] = (C[n] * GP[n + 1] + one) * inverse_A[n];
                GP[n] = GG[n] / (one + GG[n]);
            }

            Su[centre] += Su[centre + 1] * FI[centre];
            L_diag[centre] = inverse_C[centre] / (FF[centre] + GP[centre + 1]);
        } else {
            L_diag[centre] = (one + FF[centre - 1]) / (Bl_min_Al + Bl * FF[centre - 1]);
        }
    } else {
        /// Write economically: G[last] = (B[last] - A[last]) /
        /// A[last];
        GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
        GI[last] = one / (one + GG[last]);
        GP[last] = GG[last] * GI[last];

        L_diag[last] = (one + FF[last - 1]) / (Bl_min_Al + Bl * FF[last - 1]);

        for (long n = last - 1; n > first; n--) // use long in reverse loops!
        {
            Su[n] += Su[n + 1] * FI[n];

            GG[n] = (C[n] * GP[n + 1] + one) * inverse_A[n];
            GI[n] = one / (one + GG[n]);
            GP[n] = GG[n] * GI[n];

            L_diag[n] = inverse_C[n] / (FF[n] + GP[n + 1]);
        }

        Su[first] += Su[first + 1] * FI[first];
        L_diag[first] = (one + GG[first + 1]) / (Bf_min_Cf + Bf * GG[first + 1]);

        for (long n = last - 1; n >= first; n--) // use long in reverse loops!
        {
            L_upper(0, n + 1) = L_diag[n + 1] * FI[n];
            L_lower(0, n)     = L_diag[n] * GI[n + 1];
        }

        for (Size m = 1; (m < n_off_diag) && (m < n_tot - 1); m++) {
            for (long n = last - 1 - m; n >= first; n--) // use long in reverse loops!
            {
                L_upper(m, n + m + 1) = L_upper(m - 1, n + m + 1) * FI[n];
                L_lower(m, n)         = L_lower(m - 1, n) * GI[n + m + 1];
            }
        }
    }
}

///   Computes the optical depth assuming only a single line
///   exists.
///    @param[in] curridx: index of the current point
///    @param[in] nextidx: index of the next point
///    @param[in] currpoint_interp_idx: index of the interpolation point for the current
///    point
///    @param[in] nextpoint_interp_idx: index of the interpolation point for the next point
///    @param[in] lineidx: index of the line for which to
///    compute the optical depth
///    @param[in] curr_freq: current frequency (in comoving
///    frame)
///    @param[in] next_freq: next frequency (in comoving
///    frame)
///    @param[in] curr_interp: interpolation factor for the current point
///    @param[in] next_interp: interpolation factor for the next point
///    @param[in] dz: distance increment
inline Real Solver ::compute_dtau_single_line(Model& model, Size curridx, Size nextidx,
    Size currpoint_interp_idx, Size nextpoint_interp_idx, Size lineidx, Real curr_freq,
    Real next_freq, Real curr_interp, Real next_interp, Real dz) {
    const Real linefreq = model.lines.line[lineidx];
    const Real inv_width_curr_point =
        interp_helper.interpolate_linear(model.lines.inverse_width(curridx, lineidx),
            model.lines.inverse_width(currpoint_interp_idx, lineidx), curr_interp);
    const Real inv_width_next_point =
        interp_helper.interpolate_linear(model.lines.inverse_width(nextidx, lineidx),
            model.lines.inverse_width(nextpoint_interp_idx, lineidx), next_interp);
    const Real average_inverse_line_width = (inv_width_curr_point + inv_width_next_point) / 2.0;

    // opacity is stored divided by the linefreq, so multiply
    // by it
    const Real curr_line_opacity =
        interp_helper.interpolate_log(model.lines.opacity(curridx, lineidx),
            model.lines.opacity(currpoint_interp_idx, lineidx), curr_interp);
    const Real next_line_opacity =
        interp_helper.interpolate_log(model.lines.opacity(nextidx, lineidx),
            model.lines.opacity(nextpoint_interp_idx, lineidx), next_interp);

    // if frequencies are equal, division by zero (due to the
    // optical depth formula) happens if we were not to use
    // this branch
    if (curr_freq == next_freq) {
        std::cout << "Warning: Frequencies are equal, using default computation for optical depth"
                  << std::endl;
        // doing the default computation instead (no shifting)
        const Real diff = curr_freq - model.lines.line[lineidx]; // curr_freq==next_freq,
                                                                 // so choice is arbitrary
        const Real prof            = gaussian(average_inverse_line_width, diff);
        const Real average_opacity = (curr_line_opacity + next_line_opacity) / 2.0;

        return dz * (prof * average_opacity + model.parameters->min_opacity);
    }

    // We assume a linear interpolation of these dimensionless
    // frequency positions We will also assume the line width
    // to be somewhat constant, replacing the values with the
    // averages
    const Real next_pos = (linefreq - next_freq) * average_inverse_line_width;
    const Real curr_pos = (linefreq - curr_freq) * average_inverse_line_width;

    // In this way, the diff_pos can be computed quite simple,
    // and we do not have a discrepancy between the
    // interpolation and the bounds
    const Real diff_pos = next_pos - curr_pos;

    /// the more correct approach, taking into account also
    /// the line opacity change; however, it does not make too
    /// much of a difference in the actual result and is quite
    /// a bit slower
    // const Real
    // delta_opacity=(next_line_opacity-curr_line_opacity);
    // const Real deltanu=-next_freq+curr_freq;//differences
    // in curr freqs; +-1 due to shift being defined in the
    // other direction
    //
    // //note: opacity can also be extrapolated; however the
    // correction term (expterm) accounts for that const Real
    // interp_opacity=curr_line_opacity+delta_opacity*(curr_freq-linefreq)/deltanu;
    //
    // //This term is a constant term, giving the usual ... as
    // if the opacity were static const Real
    // erfterm=interp_opacity/diff_pos/2.0*(std::erf(next_pos)-std::erf(curr_pos));
    // //This term corrects for the fact that the opacity
    // between points changes const Real
    // expterm=delta_opacity/2.0*INVERSE_SQRT_PI/diff_pos/diff_pos*(std::exp(-curr_pos*curr_pos)-std::exp(-next_pos*next_pos));
    // return
    // dz*std::max(average_inverse_line_width*(erfterm+expterm),
    // model.parameters->min_opacity);

    // If we instead use an average opacity, the computation
    // is quite a bit faster
    const Real average_opacity = (next_line_opacity + curr_line_opacity) / 2.0;
    const Real erfterm =
        average_opacity / diff_pos / 2.0 * (std::erff(next_pos) - std::erff(curr_pos));
    // correcting to bound opacity from below to the minimum
    // opacity (assumes positive opacities occuring in the
    // model)
    return dz * (average_inverse_line_width * erfterm + model.parameters->min_opacity);
}

///  Solver for Feautrier equation along ray pairs using the
///  (ordinary) 2nd-order solver, without adaptive optical
///  depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::image_feautrier_order_2(Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ            = dZ_();
    Vector<Size>& nr              = nr_();
    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();
    Vector<double>& shift         = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], nr_interp[first],
        nr_interp[first + 1], l, freq * shift[first], freq * shift[first + 1], shift[first],
        shift[first + 1], dZ[first], interp_factor[first], interp_factor[first + 1],
        compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) /
    /// C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], nr_interp[n], nr_interp[n + 1], l,
            freq * shift[n], freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], interp_factor[n],
            interp_factor[n + 1], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n]        = dtau_avg * dtau_c;
        inverse_C[n]        = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source
        /// function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;

    // for (long n = last-1; n > first; n--) // use long in
    // reverse loops!
    // {
    //     Su[n] += Su[n+1] * FI[n];
    // }

    // Su[first] += Su[first+1] * FI[first];
}

///  Solver for Feautrier equation along ray pairs using the
///  (ordinary) 2nd-order solver, without adaptive optical
///  depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::image_feautrier_order_2_for_point_loc(
    Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ            = dZ_();
    Vector<Size>& nr              = nr_();
    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();
    Vector<double>& shift         = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], nr_interp[first],
        nr_interp[first + 1], l, freq * shift[first], freq * shift[first + 1], shift[first],
        shift[first + 1], dZ[first], interp_factor[first], interp_factor[first + 1],
        compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    // err, source function might be slightly different when
    // looking at it from curr and next point
    //  this is due to the weighting by the line optical
    //  depths; this might not be saved TODO think whether
    //  this is correct
    model.S_ray(0, f) = term_c;
    model.S_ray(1, f) = term_n;

    model.dtau_ray(0, f) = 0.0;
    model.dtau_ray(1, f) = dtau_n;

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) /
    /// C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], nr_interp[n], nr_interp[n + 1], l,
            freq * shift[n], freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], interp_factor[n],
            interp_factor[n + 1], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n]        = dtau_avg * dtau_c;
        inverse_C[n]        = dtau_avg * dtau_n;

        model.S_ray(n + 1 - first, f)    = term_n;
        model.dtau_ray(n + 1 - first, f) = dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source
        /// function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;

    model.u_ray(last - first, f) = Su[last];

    for (long n = last - 1; n > first; n--) // use long in reverse loops!
    {
        Su[n] += Su[n + 1] * FI[n];

        model.u_ray(n - first, f) = Su[n];
    }

    Su[first] += Su[first + 1] * FI[first];

    model.u_ray(0, f) = Su[first];
}

/// Image the optical depth
///////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::image_optical_depth(Model& model, const Size o, const Size f) {
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ            = dZ_();
    Vector<Size>& nr              = nr_();
    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();
    Vector<double>& shift         = shift_();

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    Real tau = 0.0;

    bool compute_curr_opacity = true; // for the first point, we need to compute both
                                      // the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], nr_interp[first],
        nr_interp[first + 1], l, freq * shift[first], freq * shift[first + 1], shift[first],
        shift[first + 1], dZ[first], interp_factor[first], interp_factor[first + 1],
        compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);
    tau += dtau_n;

    /// iterate over all other points on the ray
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], nr_interp[n], nr_interp[n + 1], l,
            freq * shift[n], freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], interp_factor[n],
            interp_factor[n + 1], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);
        tau += dtau_n;
    }

    optical_depth_() = tau;
}

///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///  Computes both the mean intensity u and the flux v
//////////////////////////////////////////////////////////////////////////
template <ApproximationType approx>
accel inline void Solver ::solve_feautrier_order_2_uv(Model& model, const Size o, const Size f) {
    // Note: in comments, we have another option for computing v (Sv at every location), by
    // mirroring the solution steps for u with slightly different data. This might be
    // slower, and a bit more work to make second-order accurate. Currently, only sv at the
    // center position gets computed using du/dtau = -v
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_();
    const Size n_tot = n_tot_();

    Vector<double>& dZ            = dZ_();
    Vector<Size>& nr              = nr_();
    Vector<Size>& nr_interp       = nr_interp_();
    Vector<double>& interp_factor = interp_factor_();
    Vector<double>& shift         = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_();
    Vector<Real>& C         = C_();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both the
                                      // curr and next opacity(and source)

    compute_source_dtau<approx>(model, nr[first], nr[first + 1], nr_interp[first],
        nr_interp[first + 1], l, freq * shift[first], freq * shift[first + 1], shift[first],
        shift[first + 1], dZ[first], interp_factor[first], interp_factor[first + 1],
        compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first]         = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = one / C[first]; // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity(model, nr[first], freq * shift[first]);
    // TODO: if re-implementing other option, make this dSdtau second order accurate
    // Real dSdtau          = (term_n - term_c) / dtau_n;
    // only first order accurate dSdtau, as otherwise the solver needs to be rewritten
    // Current algorithm only allows us to access the data at current and next point

    Su[first] = term_c + two * I_bdy_f * inverse_dtau_f;
    // Sv[first] = -dSdtau + two * inverse_dtau_f * (I_bdy_f - term_c);

    Su[first] /= Bf;
    // Sv[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first + 1; n < last; n++) {
        term_c = term_n;
        dtau_c = dtau_n;
        eta_c  = eta_n;
        chi_c  = chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n + 1], nr_interp[n], nr_interp[n + 1], l,
            freq * shift[n], freq * shift[n + 1], shift[n], shift[n + 1], dZ[n], interp_factor[n],
            interp_factor[n + 1], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        // TODO: if re-implementing other option, make this dSdtau second order accurate
        // dSdtau              = (term_n - term_c) / dtau_n;
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n - 1] * FI[n - 1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n - 1] + Su[n]) * FI[n] * inverse_C[n];
        // Sv[n] = (A[n] * Sv[n - 1] - dSdtau) * FI[n] * inverse_C[n];
    }

    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last - 1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity(model, nr[last], freq * shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    // Sv[last] = -dSdtau - two * inverse_dtau_l * (I_bdy_l - term_n);
    // Different sign for Sv last boundary condition extra term (2/dtau(I-S))! (should be
    // assymetric)

    Su[last] = (A[last] * Su[last - 1] + Su[last]) * (one + FF[last - 1]) * denominator;
    // Sv[last] = (A[last] * Sv[last - 1] + Sv[last]) * (one + FF[last - 1]) * denominator;

    if (centre < last) {
        for (long n = last - 1; n >= centre; n--) // use long in reverse loops !
        {
            Su[n] += Su[n + 1] * FI[n];
            // Sv[n] += Sv[n + 1] * FI[n];
        }
    } else {
        // Compute v using boundary at the end
        Sv[last] = Su[last] - I_bdy_l;
        return;
    }
    // Do one extra step for computing the derivative of the mean intensity
    if (centre > first) {
        Su[centre - 1] += Su[centre] * FI[centre - 1];
        // Recompute the optical depth increments, as we forgot to save them; FIXME: add
        // threadprivate dtau_() to solver.hpp
        compute_curr_opacity = true;
        compute_source_dtau<approx>(model, nr[centre - 1], nr[centre], nr_interp[centre - 1],
            nr_interp[centre], l, freq * shift[centre - 1], freq * shift[centre], shift[centre - 1],
            shift[centre], dZ[centre - 1], interp_factor[centre - 1], interp_factor[centre],
            compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);
        const Real dtaumin = dtau_n;
        compute_source_dtau<approx>(model, nr[centre], nr[centre + 1], nr_interp[centre],
            nr_interp[centre + 1], l, freq * shift[centre], freq * shift[centre + 1], shift[centre],
            shift[centre + 1], dZ[centre], interp_factor[centre], interp_factor[centre + 1],
            compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);
        const Real dtauplus = dtau_n;
        // TODO: optimize this calculation
        const Real coeffmin  = -dtauplus / (dtaumin * dtaumin + dtaumin * dtauplus);
        const Real coeffplus = dtaumin / (dtaumin * dtauplus + dtauplus * dtauplus);
        const Real coeffzero = -coeffmin - coeffplus;

        Sv[centre] =
            -(coeffmin * Su[centre - 1] + coeffplus * Su[centre + 1] + coeffzero * Su[centre]);
    } else {
        // Compute v using boundary at the start
        Sv[first] = I_bdy_f - Su[first];
        return;
    }
}

accel inline void Solver ::set_eta_and_chi(Model& model, const Size rr) const {
    model.eta.resize(model.parameters->npoints(), model.parameters->nfreqs());
    model.chi.resize(model.parameters->npoints(), model.parameters->nfreqs());

    for (Size p = 0; p < model.parameters->npoints(); p++) {
        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            // Extract the Doppler shift
            const double shift = model.geometry.get_shift<Rest, false>(0, rr, p, 0.0);
            const Real freq    = model.radiation.frequencies.nu(0, f);
            const Size l       = model.radiation.frequencies.corresponding_line[f];

            get_eta_and_chi<None>(model, p, l, freq * shift, model.eta(p, f), model.chi(p, f));
        }
    }
}

accel inline void Solver ::set_boundary_condition(Model& model) const {
    model.boundary_condition.resize(model.parameters->nboundary(), model.parameters->nfreqs());

    for (Size b = 0; b < model.parameters->nboundary(); b++) {
        const Size p = model.geometry.boundary.boundary2point[b];

        for (Size f = 0; f < model.parameters->nfreqs(); f++) {
            const Real freq = model.radiation.frequencies.nu(0, f);

            model.boundary_condition(b, f) = boundary_intensity(model, p, freq);
        }
    }
}

template <bool use_adaptive_directions> inline void Solver ::set_column(Model& model) const {
    model.column.resize(model.parameters->nrays(), model.parameters->npoints());

    for (Size rr = 0; rr < model.parameters->hnrays(); rr++) {

        cout << "--- rr = " << rr << endl;

        accelerated_for(o, model.parameters->npoints(), {
            const Size ar       = model.geometry.rays.get_antipod_index(rr);
            model.column(rr, o) = get_column<use_adaptive_directions>(model, o, rr);
            model.column(ar, o) = get_column<use_adaptive_directions>(model, o, ar);
        })
    }
}

template <bool use_adaptive_directions>
accel inline Real Solver ::get_column(const Model& model, const Size o, const Size r) const {
    Real column = 0.0;

    double Z  = 0.0; // distance from origin (o)
    double dZ = 0.0; // last increment in Z

    Size nxt = model.geometry.get_next<use_adaptive_directions>(o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt)) {
        Size crt = o;

        column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;

        while (model.geometry.not_on_boundary(nxt)) {
            crt = nxt;
            nxt = model.geometry.get_next<use_adaptive_directions>(o, r, nxt, Z, dZ);

            column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;
        }
    }

    return column;
}
