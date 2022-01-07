#include <cmath>
#include <Eigen/Core>

inline void Collocation :: set_interacting_bases(Model& model)
{
  //Set the kd_tree
  Eigen::Matrix<double, Eigen::Dynamic, 3> temp_mat;
  // Size1 index_conversion;
  auto tuple=model.create_mat_for_kd_tree_of_lvl(0);
  temp_mat=std::get<0>(tuple);
  index_conversion=std::get<1>(tuple);
  // The reverse map of index_conversion
  // reverse_index_conversion
  Size temp_counter=0;
  for (Size index: index_conversion)
  {
      reverse_index_conversion.insert(std::pair<Size,Size>(index,temp_counter));
      temp_counter++;
  }

  kd_tree mat_index(3, temp_mat, 10 /* max leaf */);
  mat_index.index->buildIndex();

    for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    {
        vector<Size> neighbors_coarser_grid=model.get_coarser_neighbors_kd_tree(pointid, mat_index, index_conversion);
        Real max_dist=0.0;
        Size farthest_point_idx=parameters.npoints();
        // Size furthest_point_idx=parameters.npoints();
        // TODO: prune furthest point from this
        for (Size coarser_neighbor: neighbors_coarser_grid)
        {
            Vector3D temp_vector=(model.geometry.points.position[coarser_neighbor]-model.geometry.points.position[pointid]);
            Real temp_dist=std::sqrt(temp_vector.squaredNorm());
            // if (temp_dist>max_dist)
            // {
            max_dist=temp_dist;
            farthest_point_idx=coarser_neighbor;
            // }
        }
        // Vector3D temp_vector=(model.geometry.points.position[neighbors_coarser_grid.back()]-model.geometry.points.position[pointid]);
        rbf_radii[pointid]=max_dist;
        std::cout<<"radius: "<<rbf_radii[pointid]<<std::endl;
        for (Size coarser_neighbor: neighbors_coarser_grid)
        {
            if (coarser_neighbor!=farthest_point_idx)
            {//set for every neighbors that they have a non-zero value in the rbf of pointid
            other_radial_basis_interacting_with[coarser_neighbor].insert(pointid);//this also includes itself
            }
        }

    }

    for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    {

        //currently, i am choosing the radius such that the last point lies on the boundary of the domain of the radial basis function
        //Therefore, we do not need to include it in the vector of interacting basis functions
        // vector<Size> subvector={neighbors_coarser_grid.begin(),neighbors_coarser_grid.end()-1};


        //calculate the max doppler shift possible
        Real max_doppler_shift=1.0;
        // vector<Real> min_doppler_shift(parameters.nrays(),1.0);
        Vector3D velocity_curr_point=model.geometry.points.velocity[pointid];
        // for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
        // {
        // for (Size neighbor: subvector)
        //TODO: possibly change this for the comoving frame...
        for (Size neighbor: other_radial_basis_interacting_with[pointid])
        {
            max_doppler_shift=std::max(max_doppler_shift,Real(
        1+std::sqrt((model.geometry.points.velocity[neighbor]-velocity_curr_point).squaredNorm())));
            // min_doppler_shift[rayidx]=std::min(max_doppler_shift,
        // 1+std::sqrt((model.geometry.points.velocity[neighbor]-velocity_curr_point).dot(model.geometry.rays.direction[rayidx])/CC);
        }


        // For every frequency, calculate which frequencies lie nearby
        //now using all the quadrature frequencies
        // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
        //TODO: use more efficient estimates?
        for (Size freqid=0; freqid<parameters.nfreqs(); freqid++)
        {
            // const Size corresp_lineid=?; for more efficient estimate, just try whether the line centers lie close enough
            //The line frequencies are not ordered, so just search the whole list...
            // for (Size other_lineid=0; other_lineid<parameters.nlines(); other_lineid++)
            for (Size other_freqid=0; other_freqid<parameters.nfreqs(); other_freqid++)
            {
                Size l=model.radiation.frequencies.corresponding_l_for_spec[other_freqid];   ///< number of line species corresponding to frequency
                Size k=model.radiation.frequencies.corresponding_k_for_tran[other_freqid];   ///< number of transition corresponding to frequency
                Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];   ///< number of line number corresponding to frequency
                // std::cout<<"l, k, z: "<<l<<", "<<k<<", "<<z<<std::endl;
                Size corresp_other_lineid=model.lines.line_index(l, k);
                // Somewhat loose uniform (over direction and neighbors) estimate
                // Start by bounding the apparent frequency nu' for the point (pointidx, lineid): nu/max_dpplr_shift<nu'<nu*max_dpplr_shift
                // if ((model.lines.line[lineid]/max_doppler_shift<model.lines.line[other_lineid]+1/model.lines.inverse_width(pointid,other_lineid)*TRUNCATION_SIGMA*max_doppler_shift)
                // &&  (model.lines.line[lineid]*max_doppler_shift>model.lines.line[other_lineid]-1/model.lines.inverse_width(pointid,other_lineid)*TRUNCATION_SIGMA*max_doppler_shift))
                //TODO: other possibility is to use non_doppler_shifted_inverse_widths instead of lines.inverse_width
                if ((model.radiation.frequencies.nu(pointid,freqid)/max_doppler_shift<model.radiation.frequencies.nu(pointid,other_freqid)+1/model.lines.inverse_width(pointid,corresp_other_lineid)*TRUNCATION_SIGMA*max_doppler_shift)
                &&  (model.radiation.frequencies.nu(pointid,freqid)*max_doppler_shift>model.radiation.frequencies.nu(pointid,other_freqid)-1/model.lines.inverse_width(pointid,corresp_other_lineid)*TRUNCATION_SIGMA*max_doppler_shift))
                {
                    other_frequency_basis_interacting_with[pointid][freqid].insert(other_freqid);
                }
            }
        }


        // for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        // {
        //     other_frequency_basis_interacting_with[rayid][pointid].resize(subvector.size());
        //     for (Size temp=0; temp<subvector.size(); temp++)
        //     {
        //         other_frequency_basis_interacting_with[rayid][pointid][temp].resize(parameters.nfreqs());
        //     }
        // }
        // also set what interactions there are between the frequency basis functions, which unfortunately also depend on the interacting points
        // Size reduced_index=0;
        // for (Size other_point: subvector)
        // {
        //     for (Size lineid=0; lineid<parameters.nfreqs(); lineid++)
        //     {
        //         for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        //         {
        //             Real curr_freq=frequencies[rayid][pointid][lineid];
        //             vector<Real> searchfreqs=frequencies[rayid][other_point];
        //             for (Size search_freq_idx=0; search_freq_idx<parameters.nfreqs(); search_freq_idx++)
        //             {   // The frequency basis function at the other point and search frequency is only nonzero for the frequency frequencies[pointid][lineid] if the main frequency of that gaussian does not lie to far from it (compared to the width)
        //                 Real abs_freq_diff=std::abs(curr_freq-searchfreqs[search_freq_idx]);
        //                 if (!freq_cutoff_condition(abs_freq_diff,frequencies_widths[rayid][other_point][search_freq_idx]))
        //                 {
        //                     other_frequency_basis_interacting_with[rayid][pointid][reduced_index][lineid].insert(search_freq_idx);
        //                 }
        //             }
        //         }
        //     }
        //     reduced_index++;
        // }
    }

}


/// Dev note: in geomeetry.tpp, similar functions are defined; however, these functions did not state exactly how the ray direction and the velocity difference were specified
///  @param[in] velocity_of_distant_object: the velocity of the distant object minus your local velocity
///  @param[in] ray_direction: the direction of the ray coming towards you
inline Real calculate_doppler_shift_observer_frame(Vector3D velocity_of_distant_object, Vector3D ray_direction)
{
    return 1+(velocity_of_distant_object.dot(ray_direction));//velocity is already divided by the light speed
}



// inline void Collocation :: setup_eval_matrix(Model& model)
// {
//     vector<Triplet<Real, Size>> eval_triplets; //triplets for the evaluation matrix
//     for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
//     {
//         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
//         {
//             for (Size pointidx=0; pointidx<point_locations.size(); pointidx++)
//             {
//                 //Now creating all the triplets, thus we need all nonzero basis triplets
//                 std::set<std::vector<Size>> nonzero_triplets;
//                 get_nonzero_basis_triplets(nonzero_triplets, rayidx, lineidx, pointidx);
//
//                 const Size curr_mat_idx=get_mat_index(rayidx, lineidx, pointidx);
//
//
//             }
//         }
//     }
// }

inline void Collocation :: setup_basis(Model& model)
{

    //then we also set all the radii and locations for the radial basis functions
    //currently just uses all the points // might need to be changed later to all the points in a certain grid // note: currently also adding the boundary points
    point_locations.resize(parameters.npoints());
    rbf_radii.resize(parameters.npoints());
    other_radial_basis_interacting_with.resize(parameters.npoints());
    //We need all frequencies at ray and each point (can be different at each point due to doppler shift)
    frequencies.resize(parameters.nrays());
    frequencies_inverse_widths.resize(parameters.nrays());
    other_frequency_basis_interacting_with.resize(parameters.npoints());  //?with a bit more computational effort, this could be reduced a bit in size
    // for (Size i:model.geometry.points.multiscale.get_current_points_in_grid())
    for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    {
        point_locations[pointid]=model.geometry.points.position[pointid];
        // other_frequency_basis_interacting_with[pointid].resize(parameters.nlines());
        other_frequency_basis_interacting_with[pointid].resize(parameters.nfreqs());
    }

    // We need the inverse line widths, so we might as well make sure that they are calculated now//TODO: remove this here and just say that we need them before
    model.compute_inverse_line_widths();

    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        frequencies[rayid].resize(parameters.npoints());
        frequencies_inverse_widths[rayid].resize(parameters.npoints());

        for (Size pointid=0; pointid<parameters.npoints(); pointid++)
        {
            // vector<Real> non_shifted_frequencies=model.lines.line;
            // vector<Real> non_shifted_inverse_widths=model.lines.inverse_width[pointid];
            Real doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointid],model.geometry.rays.direction[rayid]);
            // frequencies[rayid][pointid].resize(parameters.nlines());
            // frequencies_inverse_widths[rayid][pointid].resize(parameters.nlines());
            frequencies[rayid][pointid].resize(parameters.nfreqs());
            frequencies_inverse_widths[rayid][pointid].resize(parameters.nfreqs());

            // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
            for (Size freqid=0; freqid<parameters.nfreqs(); freqid++)
            {
                // frequencies[rayid][pointid][lineid]=doppler_shift*model.lines.line[lineid];
                frequencies[rayid][pointid][freqid]=doppler_shift*model.radiation.frequencies.nu[freqid];
                Size l=model.radiation.frequencies.corresponding_l_for_spec[freqid];   ///< number of line species corresponding to frequency
                Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];   ///< number of transition corresponding to frequency
                Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];   ///< number of line number corresponding to frequency
                // std::cout<<"l, k, z: "<<l<<", "<<k<<", "<<z<<std::endl;
                Size corresp_lineid=model.lines.line_index(l, k);
                // model.lines.lineProducingSpecies[l].nr_line[pointid][k][z];
                // std::cout<<"corresp lineid: "<<corresp_lineid<<std::endl;
                frequencies_inverse_widths[rayid][pointid][freqid]=model.lines.inverse_width(pointid,corresp_lineid)/doppler_shift;
            }

        }
    }

    // non_doppler_shifted_frequencies.resize(parameters.nlines());
    // non_doppler_shifted_inverse_widths.resize(parameters.nlines());
    // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
    // {
    //     non_doppler_shifted_frequencies[lineid]=model.lines.line[lineid];
    //     non_doppler_shifted_inverse_widths[lineid].resize(parameters.npoints());
    //     for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    //     {
    //         non_doppler_shifted_inverse_widths[lineid][pointid]=model.lines.inverse_width(pointid,lineid);
    //     }
    // }

    non_doppler_shifted_frequencies.resize(parameters.nfreqs());
    non_doppler_shifted_inverse_widths.resize(parameters.nfreqs());
    for (Size freqid=0; freqid<parameters.nfreqs(); freqid++)
    {
        non_doppler_shifted_frequencies[freqid]=model.radiation.frequencies.nu[freqid];
        std::cout<<"non doppler shifted frequencies: "<<non_doppler_shifted_frequencies[freqid]<<std::endl;
        non_doppler_shifted_inverse_widths[freqid].resize(parameters.npoints());
        for (Size pointid=0; pointid<parameters.npoints(); pointid++)
        {
            Size l=model.radiation.frequencies.corresponding_l_for_spec[freqid];   ///< number of line species corresponding to frequency
            Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];   ///< number of transition corresponding to frequency
            Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];   ///< number of line number corresponding to frequency
            // std::cout<<"l, k, z: "<<l<<", "<<k<<", "<<z<<std::endl;
              Size corresp_lineid=model.lines.line_index(l, k);
            // std::cout<<"corresp lineid: "<<corresp_lineid<<std::endl;
            non_doppler_shifted_inverse_widths[freqid][pointid]=model.lines.inverse_width(pointid,corresp_lineid);
            // std::cout<<"inverse_width: "<<non_doppler_shifted_inverse_widths[freqid][pointid]<<std::endl;
        }
    }


    // Also do not forget to initialize the basis coefficients // TODO: figure out which exact format (std::vector, Eigen, ...) to use
    // basis_coefficients=VectorXr::Zero(parameters.nrays()*parameters.nlines()*parameters.npoints());
    basis_coefficients=VectorXr::Zero(parameters.nrays()*parameters.nfreqs()*parameters.npoints());

    set_interacting_bases(model);
    //Now determine which basis's interact with eachother

    //err, it turns out I first need to determine the typical rbf radius in order to determine whether we need a slope...
    // So that is why this seems out of place.
    //TODO: when finally deciding to remove the observer frame formulation, please try to remove the ray dependence
    rbf_using_slope.resize(parameters.nrays());
    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        rbf_using_slope[rayid].resize(parameters.npoints());

        for (Size pointid=0; pointid<parameters.npoints(); pointid++)
        {
            const Real point_radius=rbf_radii[pointid];
            rbf_using_slope[rayid][pointid].resize(parameters.nfreqs());

            for (Size freqid=0; freqid<parameters.nfreqs(); freqid++)
            {
                const Real point_opacity=get_opacity(model, rayid, freqid, pointid);
                //TODO: maybe do something else than a binary choice between slope and no slope
                if (point_opacity*point_radius>SLOPE_THRESHOLD)
                {
                    rbf_using_slope[rayid][pointid][freqid]=false;
                }
                else
                {
                    rbf_using_slope[rayid][pointid][freqid]=true;
                }
                std::cout<<"rid, pid, fid: "<<rayid<<", "<<pointid<<", "<<freqid<<std::endl;
                std::cout<<"using slope: "<<rbf_using_slope[rayid][pointid][freqid]<<std::endl;
            }
        }
    }

    // set_eval_matrix(model);
}


inline bool Collocation :: freq_cutoff_condition(Real abs_freq_diff, Real freq_basis_inverse_width)
{
    return false;//debug
    return (abs_freq_diff*freq_basis_inverse_width>TRUNCATION_SIGMA);
}

/// The basis function associated with each direction
inline Real Collocation :: basis_direction(Size rayindex)
{
    //as we do not need directional derivatives, piecewise constant basis functions are sufficient;
    return 1;
}

/// The basis function associated with each direction integrated over the solid angle
inline Real Collocation :: basis_direction_int(Geometry& geometry, Size rayindex)
{
    return FOUR_PI*geometry.rays.weight[rayindex];
}

///  The basis function associated with each frequency, evaluated at frequency currnu
inline Real Collocation :: basis_freq(Size rayidx, Size freqidx, Size pointidx, Real currfreq)
{
    //TODO: figure out what to use exactly
    //also implement some logical data structure (possibly first read out all the sorted frequencies, then construct the basis functions)
    Real abs_freq_diff=0;
    Real freq_inv_width=0;

    // Real background=boundary_intensity(Model& model, Size rayidx, Size freqidx, Size pointidx);
    if (USING_COMOVING_FRAME)
    {
        freq_inv_width=non_doppler_shifted_inverse_widths[freqidx][pointidx];//TEST:/2 for testing wider function
        abs_freq_diff=std::abs(currfreq-non_doppler_shifted_frequencies[freqidx]-FREQ_OFFSET/freq_inv_width);
        // std::cout<<"freq inverse width: "<<freq_inv_width<<std::endl;
        // std::cout<<"shift: "<<FREQ_OFFSET/freq_inv_width<<std::endl;
        // std::cout<<"base freq: "<<non_doppler_shifted_frequencies[freqidx]<<std::endl;
        // std::cout<<"abs_freq_diff: "<<abs_freq_diff<<std::endl;
        // std::cout<<"freqidx: "<<freqidx<<std::endl;
    }
    else
    {
        abs_freq_diff=std::abs(currfreq-frequencies[rayidx][pointidx][freqidx]);
        freq_inv_width=frequencies_inverse_widths[rayidx][pointidx][freqidx];
    }

    //FIXME: freq offset somewhat makes the cutoff condition a bit nonsymmetrical
    if (freq_cutoff_condition(abs_freq_diff,freq_inv_width))
    {
        return 0;
    }
    else
    {
        // std::cout<<"basis_freq: "<<INVERSE_SQRT_PI*frequencies_inverse_widths[rayidx][pointidx][lineidx]*std::exp(-std::pow(abs_freq_diff*frequencies_inverse_widths[rayidx][pointidx][lineidx],2))<<std::endl;
        return INVERSE_SQRT_PI*freq_inv_width*std::exp(-std::pow(abs_freq_diff*freq_inv_width,2));
    }
}

inline Real Collocation :: basis_freq_lp_int(Size rayidx, Size freqidx, Size pointidx, Real lp_freq, Real lp_freq_inv_width)
{
    Real point_freq=0;//frequencies[rayidx][pointidx][lineidx];
    Real point_freq_inv_width=0;//frequencies_inverse_widths[rayidx][pointidx][lineidx];

    // Real abs_freq_diff=0;
    // Real freq_inv_width=0;
    if (USING_COMOVING_FRAME)
    {
        point_freq_inv_width=non_doppler_shifted_inverse_widths[freqidx][pointidx];//TEST:/2 for testing wider function
        point_freq=non_doppler_shifted_frequencies[freqidx]+FREQ_OFFSET/point_freq_inv_width;
    }
    else
    {
        point_freq=frequencies[rayidx][pointidx][freqidx];
        point_freq_inv_width=frequencies_inverse_widths[rayidx][pointidx][freqidx];
    }

    const Real inv_variance=1/(1/std::pow(point_freq_inv_width,2)+1/std::pow(lp_freq_inv_width,2));// note: still includes the factor (1/sqrt(2)) squared (because variance)
    const Real delta_freq=lp_freq-point_freq;

    // std::cout<<"rel_delta_freq: "<<point_freq_inv_width*delta_freq<<std::endl;
    // std::cout<<"freq int: "<<INVERSE_SQRT_PI*std::sqrt(2)*std::sqrt(inv_variance)*std::exp(-(std::pow(delta_freq,2)*inv_variance))<<std::endl;
    //NOTE: factor *std::sqrt(2) is maybe necessary, but i don't know for sure... TODO: recalculate this analytically
    return INVERSE_SQRT_PI*std::sqrt(inv_variance)*std::exp(-(std::pow(delta_freq,2)*inv_variance));

}




///  The derivative of the frequency basis function
inline Real Collocation :: basis_freq_der(Size rayidx, Size freqidx, Size pointidx, Real currfreq)
{
    Real abs_freq_diff=0;
    Real freq_diff=0;
    Real freq_inv_width=0;
    if (USING_COMOVING_FRAME)
    {
        freq_inv_width=non_doppler_shifted_inverse_widths[freqidx][pointidx];//TEST:/2 for testing wider function
        freq_diff=currfreq-non_doppler_shifted_frequencies[freqidx]-FREQ_OFFSET/freq_inv_width;
        // abs_freq_diff=std::abs(currfreq-non_doppler_shifted_frequencies[freqidx]-FREQ_OFFSET/freq_inv_width);
    }
    else
    {
        freq_diff=currfreq-frequencies[rayidx][pointidx][freqidx];
        // abs_freq_diff=std::abs(currfreq-frequencies[rayidx][pointidx][freqidx]);
        freq_inv_width=frequencies_inverse_widths[rayidx][pointidx][freqidx];
    }
    abs_freq_diff=std::abs(freq_diff);//TODO: maybe replace def of freq_cutoff_condition by accepting the freq diff instead of the abs freq diff

    if (freq_cutoff_condition(abs_freq_diff,freq_inv_width))
    {
        return 0;
    }
    else
    {
        return INVERSE_SQRT_PI*freq_inv_width*-2*(freq_diff*freq_inv_width)*freq_inv_width*std::exp(-std::pow(abs_freq_diff*freq_inv_width,2));
    }
}

//DEPRECATED
// ///  Function to manipulate the (already scaled) distance
// inline Real Collocation :: distance_manip(Size pointid, Real original_distance)
// {
//     //test option to check that everything works correctly7
//     // return original_distance;
//     //option to scale by using the n-th root
//     // return std::pow(original_distance, DISTANCE_EXPONENT);
//     //option to scale by using a truncated logarithm TODO: define eps for each point, (look at closest distance)
//     // if (original_distance<DISTANCE_EPS)
//     // {return 0;}
//     // else
//     // {return 1-std::log(original_distance)/std::log(DISTANCE_EPS);}
//     //truncated logarithm fitted exactly into the domain and range [0,1]
//     return (1-std::log(original_distance+DISTANCE_EPS)/std::log(DISTANCE_EPS))/(1-std::log(1+DISTANCE_EPS)/std::log(DISTANCE_EPS));
// }
//
// inline Real Collocation :: distance_manip_der(Size pointid, Real original_distance)
// {
//     //test option to check that everything works correctly
//     // return 1;
//     //option to scale by using the n-th root
//     // return DISTANCE_EXPONENT*std::pow(original_distance, DISTANCE_EXPONENT-1);//TODO: check if nans occur if using pow in this way (negative exponent)
//     //option to scale by using a truncated logarithm TODO: define eps for each point, (look at closest distance)
//     //if (original_distance<eps)
//     //{return 0;}
//     //else
//     //{return -1/(original_distance*std::log(eps);}
//     //truncated logarithm fitted exactly into the domain and range [0,1]
//     return -1/(original_distance+DISTANCE_EPS)/std::log(DISTANCE_EPS)/(1-std::log(1+DISTANCE_EPS)/std::log(DISTANCE_EPS));
//
// }

///Returns the slope
///   @param[in]: relative distance on ray compared to the rbf radius (should lie between -1 and 1, inclusive)
//TODO: also add the other slope candidate(s)
inline Real Collocation :: get_slope(Real x, Real rayidx, Real pointidx, Real freqidx)
{
   if (!rbf_using_slope[rayidx][pointidx][freqidx])
   {
      return 1;
   }
   else
   {
      Real slope=1.0/2.0+1.0/(1.0+std::exp(-SLOPE_STEEPNESS*x));
      return slope;
   }
}

///Returns the slope derivative
///   @param[in]: relative distance on ray compared to the rbf radius (should lie between -1 and 1, inclusive)
///NOTE: the usual factor /radius should still be multiplied manually
//TODO: also add the other slope candidate(s)
inline Real Collocation :: get_slope_derivative(Real x, Real rayidx, Real pointidx, Real freqidx)
{
    if (!rbf_using_slope[rayidx][pointidx][freqidx])
    {
        return 0;
    }
    else
    {
        Real exp_factor=std::exp(-SLOPE_STEEPNESS*x);
        Real slope_derivative=SLOPE_STEEPNESS*exp_factor/std::pow(1+exp_factor,2);
        return slope_derivative;
    }

}

///  The radial basis function for the position
inline Real Collocation :: basis_point(Size centerpoint, Vector3D& location, Size rayindex, Size freqidx, Geometry& geometry)
{
    Vector3D diff_vector=location-point_locations[centerpoint];
    Real distance=std::sqrt(diff_vector.dot(diff_vector));
    Real radius=rbf_radii[centerpoint];
    // std::cout<<"distance: "<<distance<<std::endl;
    // std::cout<<"radius: "<<radius<<std::endl;
    //if the location is too far from the center, the evalutation of the compact radial basis function is 0.
    if (distance>=radius)
    {
        return 0;
    }
    else
    {
        Real rel_dist=distance/radius;
        // Real manip_rel_dist=distance_manip(centerpoint, rel_dist);
        // std::cout<<"rel_dist: "<<rel_dist<<std::endl;
        //This basis function has nice continuous derivatives at r=0 and r=1 (being 0 at both places);
        Real basis=-4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        // Real basis=1-rel_dist;
        // Real basis=std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));
        Vector3D raydirection=geometry.rays.direction[rayindex];

        // Real x=0;
        // if (rel_dist>0)
        // {
        //     Real costheta=raydirection.dot(diff_vector)/std::sqrt(diff_vector.squaredNorm());
        //     Real x=costheta*manip_rel_dist;
        // }
        Real x=raydirection.dot(diff_vector)/radius;
        //the slope 1+x in the direction of the ray;
        // Real slope=1+x/SLOPE_FACTOR;
        //the slope 1/2+1/(1+e^(-kx))
        // Real slope=1.0/2.0+1.0/(1.0+std::exp(-SLOPE_STEEPNESS*x));
        // std::cout<<"slope: "<<slope<<std::endl;
        // Real slope=1;
        Real slope=get_slope(x, rayindex, centerpoint, freqidx);

        return basis*slope;
        // return basis;
        // return 1-rel_dist;
    }
}

///  The directional derivative of the position basis function
inline Real Collocation :: basis_point_der(Size centerpoint, Vector3D& location, Size rayindex, Size freqidx, Geometry& geometry)
{
    Vector3D diff_vector=location-point_locations[centerpoint];
    Real distance=std::sqrt(diff_vector.dot(diff_vector));
    // Real distance=distance_manip(centerpoint, std::sqrt(diff_vector.dot(diff_vector)));
    // Real manip_der_factor=distance_manip_der(centerpoint, std::sqrt(diff_vector.dot(diff_vector)));
    Real radius=rbf_radii[centerpoint];
    //if the location is too far from the center, the evalutation of the compact radial basis function is 0.
    if (distance>=radius)
    {
        return 0;
    }
    else if (distance==0)
    {
        Real rel_dist=distance/radius;
        // std::cout<<"manip der factor: "<<manip_der_factor<<std::endl;
        // Real manip_rel_dist=distance_manip(centerpoint, rel_dist);
        // Real manip_der_factor=distance_manip_der(centerpoint, 0);

        // return -2/radius;//for 1-2x+x^2 as basis function
        // Real basis=1-rel_dist;//for 1-x as basis function
        // return 0; //derivative at r=0 is defined to be zero

        // the slope 1+x/SLOPE_FACTOR
        // Real slope_derivative=1/SLOPE_FACTOR;

        // the slope 1/2+1/(1+e^(-kx))
        // Real exp_factor=std::exp(-SLOPE_STEEPNESS*x);
        // Real slope_derivative=SLOPE_STEEPNESS*std::exp(-SLOPE_STEEPNESS*x)/std::pow(1+std::exp(-SLOPE_STEEPNESS*x),2);
        //evaluated at x=0
        // Real slope_derivative=SLOPE_STEEPNESS/4;
        Real slope_derivative=get_slope_derivative(0, rayindex, centerpoint, freqidx);
        // Real slope_derivative=0;
        Real basis=-4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        // Real basis=std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));
        return slope_derivative*basis/radius;
        // return manip_der_factor*slope_derivative*basis/radius;
    }
    else
    {
        Real rel_dist=distance/radius;
        // Real manip_der_factor=distance_manip_der(centerpoint, rel_dist);
        // std::cout<<"manip der factor: "<<manip_der_factor<<std::endl;
        // Real manip_rel_dist=distance_manip(centerpoint, rel_dist);

        Vector3D raydirection=geometry.rays.direction[rayindex];
        Real costheta=raydirection.dot(diff_vector)/std::sqrt(diff_vector.squaredNorm());
        // Real x=raydirection.dot(diff_vector)/radius;
        Real x=costheta*rel_dist;
        //the slope 1+x in the direction of the ray;
        // Real slope=1+x/SLOPE_FACTOR;
        // Real slope_derivative=1/SLOPE_FACTOR;
        //the slope 1/2+1/(1+e^(-kx)) in direction of the ray
        // Real slope=1.0/2.0+1.0/(1.0+std::exp(-SLOPE_STEEPNESS*x));
        // Real exp_factor=std::exp(-SLOPE_STEEPNESS*x);
        // Real slope_derivative=SLOPE_STEEPNESS*exp_factor/std::pow(1+exp_factor,2);
        Real slope=get_slope(x, rayindex, centerpoint, freqidx);
        Real slope_derivative=get_slope_derivative(x, rayindex, centerpoint, freqidx);

        // Real slope=1;
        // Real slope_derivative=0;
        Real basis=-4.0/(1.0+std::pow(rel_dist,3))+6.0/(1.0+std::pow(rel_dist,2))-1;
        // Real basis=1-rel_dist;
        // Real basis=std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));

        // std::cout<<"raydirection=[x,y,z]: "<<raydirection.x()<<", "<<raydirection.y()<<", "<<raydirection.z()<<std::endl;
        //The derivatives in the perpendicular direction (to the radial direction) are 0, so we can calculate the derivative
        // by taking the radial derivative and multiplying it by cos(theta) in which theta is the angle between the radial and the actual direction
        //Below, we see the simple rule for calculating the cosine
        //Note: The ray directions are normalised, so we can omit the normalization of the ray direction
        // Minus sign because we actually need the gradient (not -1 times the gradient); this is due to the definition of the diff_vector, which is natural for RBF's but does actually not give us the angle we want

        //derivative of -4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        Real radial_derivative = 12*(std::pow(rel_dist,2))/(std::pow(1.0+std::pow(rel_dist,3),2))-12*rel_dist/(std::pow(1.0+std::pow(rel_dist,2),2));
        // Real slope_derivative = 1;
        //derivative of 1-rel_dist
        // Real radial_derivative = -1;
        //derivative of 1-2rel_dist+rel_dist^2
        // Real radial_derivative = -2+2*rel_dist;
        //derivative of e^-(x*trunc)^2
        // Real radial_derivative = -std::pow(TRUNCATION_SIGMA,2)*2*rel_dist*std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));
        // Because we switched the perspective of the basis and point, we need the minus sign
        // Normally you take a basis, and compute the derivative at a given point; in this case, we do the opposite:
        // we take a point and calculate the derivative at that position in some basis functions.
        // std::cout<<"der without slope: "<<costheta*radial_derivative/radius<<std::endl;
        // std::cout<<"der with slope: "<<costheta*radial_derivative*slope/radius+slope_derivative*basis/radius<<std::endl;
        //currently, the slope uses the default distance, not the manipulated one
        return costheta*radial_derivative*slope/radius+slope_derivative*basis/radius;
        // return costheta*radial_derivative/radius;
    }
}

/// Fills the triplets which correspond to the nonzero basis function products at a given triplet location
///   @param[in/out] basis_triplets_to_fill: An empty set which will contain the basis triplets after execution
inline void Collocation :: get_nonzero_basis_triplets(std::set<std::vector<Size>>& basis_triplets_to_fill, Size rayidx, Size freqidx, Size pointidx)
{
    // std::set<std::vector<Size>> basis_triplets;
    // I currently assume that the different ray direction bases do not interact with eachother
    for (Size interacting_points: other_radial_basis_interacting_with[pointidx])
    {
        for (Size interacting_freqs: other_frequency_basis_interacting_with[pointidx][freqidx])
        {   //TODO: we are constantly allocating new memory; this might be slow; better to just compute this once...
            vector<Size> temp_triplet{rayidx, interacting_freqs, interacting_points};
            basis_triplets_to_fill.insert(temp_triplet);
        }
    }
    // return basis_triplets;
}

///  Returns the matrix index corresponding to the general index (rayidx, lineidx, pointidx)
inline Size Collocation :: get_mat_index(Size rayidx, Size freqidx, Size pointidx)
{
    // return rayidx*parameters.nlines()*point_locations.size()+lineidx*point_locations.size()+pointidx;
    // return rayidx*parameters.nfreqs()*point_locations.size()+freqidx*point_locations.size()+pointidx;
    return rayidx*parameters.nfreqs()*point_locations.size()+pointidx*parameters.nfreqs()+freqidx;
}

// inline Real Collocation :: calculate_doppler_shift(Size pointidx, Size rayidx, Geometry& geometry)
// {
//     Vector3D raydirection=model.geometry.rays.direction[rayindex];
//     return 1+std::sqrt(raydirection.dot())
// }

//FIXME: reomve this, as its function is already improved upon by dividing by the diagonal.
//Computes the balancing factor needed in order to make the boundary condition equation the same magnitude as the other nearby equations in the matrix
//maybe TODO: also take the scattering stuff into account (when implemented)
inline Real Collocation :: compute_balancing_factor_boundary(Size rayidx, Size freqidx, Size pointidx, Real local_velocity_gradient, Model& model)
{
    return 1;
    // const Real numerator=basis_point_der(pointidx, point_locations[pointidx], rayidx, freqidx, model.geometry)*basis_freq(rayidx, freqidx, pointidx, non_doppler_shifted_frequencies[freqidx])-local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]*basis_point(pointidx, point_locations[pointidx], rayidx, freqidx, model.geometry)*basis_freq_der(rayidx, freqidx, pointidx, non_doppler_shifted_frequencies[freqidx]);
    // const Real denominator=basis_point(pointidx, point_locations[pointidx], rayidx, freqidx, model.geometry)*get_opacity(model, rayidx, freqidx, pointidx)*basis_freq(rayidx, freqidx, pointidx, non_doppler_shifted_frequencies[freqidx]);
    // // std::cout<<"numerator/denominator: "<<numerator/denominator<<std::endl;
    //
    // return 1 + numerator/denominator;
}


//
inline void Collocation :: setup_basis_matrix_Eigen(Model& model)
{
    //We are currently not using the lambda operator, so set to zero
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    Size mat_size=parameters.nrays()*parameters.nfreqs()*point_locations.size();
    vector<Triplet<Real, Size>> eigen_triplets;
    // eigen_triplets.reserve(mat_size*16*nrays());//matrix size times number of interacting points times the number of rays (scattering)
    eigen_triplets.reserve(mat_size*16);//matrix size times number of interacting points (does not include scattering)

    for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
    {
        std::cout<<"rayidx: "<<rayidx<<std::endl;
        std::cout<<"raydirection: "<<model.geometry.rays.direction[rayidx].x()<<std::endl;
        for (Size freqidx=0; freqidx<parameters.nfreqs(); freqidx++)
        {


            //TODO: either store all doppler shifted frequencies or calculate some relative doppler shifts
            //The doppler shift term n.grad(n.v(x)/c)*nu does not depend on the basis functions. Thus we calculate it here
            //TODO: actually check how it's calculated in Geometry
            // Real doppler_shift_term=0;//temporary value FIXME

            for (Size pointidx=0; pointidx<point_locations.size(); pointidx++)
            {
                //Now creating all the triplets, thus we need all nonzero basis triplets
                std::set<std::vector<Size>> nonzero_triplets;
                get_nonzero_basis_triplets(nonzero_triplets, rayidx, freqidx, pointidx);

                const Size curr_mat_idx=get_mat_index(rayidx, freqidx, pointidx);
                const Real curr_opacity=get_opacity(model, rayidx, freqidx, pointidx);
                // TODO: think about whether the frequency should be evaluated in the frame of the neighbors...//no, as we evaluate it at a single point each time (no interpolation stuff required)

                //TODO: maybe refactor this part
                Real local_velocity_gradient=0;

                double dist=0;
                double ddist=0;
                Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                if (next_point!=parameters.npoints())//if there actually exists a next point
                {
                    // std::cout<<"dist: "<<dist<<std::endl;
                    local_velocity_gradient=model.geometry.rays.direction[rayidx].dot(model.geometry.points.velocity[next_point]-model.geometry.points.velocity[pointidx])/dist;
                    // std::cout<<"local velocity grad: "<<local_velocity_gradient<<std::endl;
                }
                else
                {
                    double dist=0;
                    double ddist=0;
                    Size prev_point=model.geometry.get_next(pointidx, model.geometry.rays.antipod[rayidx], pointidx, dist, ddist);
                    if (prev_point!=parameters.npoints())
                    {
                        local_velocity_gradient=model.geometry.rays.direction[rayidx].dot(model.geometry.points.velocity[pointidx]-model.geometry.points.velocity[prev_point])/dist;
                    }//otherwise we only have a single point on the ray in this direction; so defining a non-zero velocity gradient seems silly
                }
                //TODO: refactor until here


                Real curr_freq=0;
                if (USING_COMOVING_FRAME)
                {
                    curr_freq=non_doppler_shifted_frequencies[freqidx];
                }
                else
                {
                    curr_freq=frequencies[rayidx][pointidx][freqidx];
                }
                Vector3D curr_location=point_locations[pointidx];
                // std::cout<<"this pointidx: "<<pointidx<<std::endl;
                // std::cout<<"curr_opacity: "<<curr_opacity<<std::endl;

                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(pointidx))
                {
                    double dist=0;//these two variables are actually not used, but required in the function definition get_next
                    double ddist=0;

                    Size antipod=model.geometry.rays.antipod[rayidx];

                    const Size point_behind=model.geometry.get_next(pointidx, antipod, pointidx, dist, ddist);
                    if (point_behind==parameters.npoints())
                    {//if there lies no point behind this point, we need to evaluate the boundary condition here
                        boundary_condition_required=true;
                    }
                }
                // boundary_condition_required=false;
                if (boundary_condition_required)
                {
                    //This term makes sure that the magnitude of the matrix terms corresponding to this boundary condition is similar to the other nearby terms
                    const Real balancing_factor=compute_balancing_factor_boundary(rayidx, freqidx, pointidx, local_velocity_gradient, model);
                    //Merely evaluating the intensity at the boundary
                    for (std::vector<Size> triplet: nonzero_triplets)
                    {
                        const Size other_mat_idx=get_mat_index(triplet[0],triplet[1],triplet[2]);

                        if (USING_COMOVING_FRAME)
                        {
                            // const Real doppler_shifted_frequency=curr_freq*calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx]
                            //                                                                                       , model.geometry.rays.direction[rayidx]);
                            // Real doppler_shifted_frequency=curr_freq;//try it out non-doppler shifted
                            // Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point(triplet[2], curr_location, triplet[0], model.geometry);
                            Real basis_eval=balancing_factor*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], triplet[1], model.geometry);
                            eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, other_mat_idx, basis_eval));
                        }
                        else
                        {
                            //In order to make this term of the same size as the others, multiply by the opacity (only in emissivity formulation)
                            // Real basis_eval=curr_opacity*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                            Real basis_eval=balancing_factor*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], triplet[1], model.geometry);
                            // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                            // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                            eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, other_mat_idx, basis_eval));
                        }
                        std::cout<<"bdy condition triplets: "<<triplet[0]<<", "<<triplet[1]<<", "<<triplet[2]<<std::endl;
                        std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                        std::cout<<"basis_point: "<<basis_point(triplet[2], curr_location, triplet[0], triplet[1], model.geometry)<<std::endl;

                    }
                }
                else
                {
                    for (std::vector<Size> triplet: nonzero_triplets)
                    {
                        const Size other_mat_idx=get_mat_index(triplet[0],triplet[1],triplet[2]);

                        if (USING_COMOVING_FRAME)
                        {   //TODO: give a more solid reasoning why the doppler shift should be calculated using the rayidx, instead of triplet[0]
                            // const Real doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx]
                            //                                                                 , model.geometry.rays.direction[rayidx]);
                            // const Real doppler_shifted_frequency=curr_freq*doppler_shift;
                            // const Real doppler_shifted_opacity=get_opacity(model, rayidx, lineidx, pointidx, doppler_shift);
                            // Real doppler_shifted_frequency=curr_freq;//try it out non-doppler shifted
                            // Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry)/curr_opacity;///doppler_shifted_opacity;///curr_opacity;
                            // Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point(triplet[2], curr_location, triplet[0], model.geometry);
                            Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], triplet[1], model.geometry)/curr_opacity;///doppler_shifted_opacity;///curr_opacity;
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], triplet[1], model.geometry);
                            std::cout<<"curr_opacity: "<<curr_opacity<<std::endl;
                            // std::cout<<"basis_eval+dir_der: "<<basis_eval+directional_der<<std::endl;
                            // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"dir_der: "<<directional_der<<std::endl;
                            // std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_point: "<<basis_point(triplet[2], curr_location, triplet[0], model.geometry)<<std::endl;
                            Real doppler_shift_term=0;//TODO: add frequency quadrature;; TODO: think about approximating velocity gradient somehow at the point itself?
                            // if (freqidx!=triplet[1])//doppler shift only applies when not at the same frequency
                            // {//FIXME: freqidx!=triplet[1] and compute the velocity gradient somehow
                                // Vector3D diff_vector=point_locations[triplet[2]]-curr_location;
                                // Real distance=std::sqrt(diff_vector.squaredNorm());
                                // Real velocity_grad=model.geometry.rays.direction[rayidx].dot(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx])/distance;
                                //divided by the light speed
                            doppler_shift_term=-local_velocity_gradient*curr_freq/curr_opacity*
                                  basis_direction(triplet[0])*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], triplet[1], model.geometry);
                            // }
                            // std::cout<<"velocity grad: "<<local_velocity_gradient<<std::endl;
                            // std::cout<<"freq: "<<non_doppler_shifted_frequencies[freqidx]<<std::endl;
                            // std::cout<<"1/opacity: "<<1.0/curr_opacity<<std::endl;
                            // std::cout<<"doppler_shift_spat_der_ratio"<<local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)/basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], model.geometry)/basis_point_der(triplet[2], curr_location, triplet[0], model.geometry)<<std::endl;
                            // std::cout<<"inverse width: "<<non_doppler_shifted_inverse_widths[triplet[1]][triplet[2]]<<std::endl;
                            // std::cout<<"doppler_shift_term: "<<doppler_shift_term<<std::endl;
                            // std::cout<<"doppler_shift_basis_ratio: "<<-local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]/curr_opacity*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)/basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // doppler_shift_term=0;

                            eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, other_mat_idx, directional_der+basis_eval+doppler_shift_term));

                            // //For every scattering direction! add scattering//also divide by curr_opacity
                            // for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                            // {
                            //     const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                            //     //TODO: add toggle to whether we are actually using scattering
                            //     Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                            //     Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                            //     //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                            //     // eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                            // }
                        }
                        else
                        {
                            //Start with the directional derivative (in source function formulation divided by the opacity)
                            // Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry);
                            Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], triplet[1], model.geometry)/curr_opacity;
                            // std::cout<<"directional_der: "<<directional_der<<std::endl;
                            // directional_der=0;
                            //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                            // Real opacity=curr_opacity*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], triplet[1], model.geometry);
                            // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                            // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_point_eval: "<<basis_point(triplet[2], curr_location)<<std::endl;
                            // basis_eval=0;
                            Real doppler_shift_term=0;
                            // doppler_shift_term=-local_velocity_gradient*curr_freq/curr_opacity*
                            //       basis_direction(triplet[0])*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], triplet[1], model.geometry);

                            // if (pointidx!=triplet[2])//doppler shift only applies when not at the same position
                            // {
                            //     Vector3D diff_vector=point_locations[triplet[2]]-curr_location;
                            //     Real distance=std::sqrt(diff_vector.squaredNorm());
                            //     Real velocity_grad=model.geometry.rays.direction[rayidx].dot(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx])/distance;
                            //     //divided by the light speed
                            //     doppler_shift_term=-velocity_grad*frequencies[rayidx][pointidx][lineidx]/curr_opacity*
                            //       basis_direction(triplet[0])*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], model.geometry);
                            // }

                            eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, other_mat_idx, directional_der+basis_eval+doppler_shift_term));

                            // //For every scattering direction! add scattering//also divide by curr_opacity
                            // for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                            // {
                            //     const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                            //     //TODO: add toggle to whether we are actually using scattering
                            //     Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                            //     Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                            //     //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                            //     // eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                            // }

                        }

                    }
                }
            }
        }
    }

    // std::cout<<"resizing collocation mat"<<std::endl;

    //after all this, construct the matrix
    collocation_mat.resize (mat_size, mat_size);

    // std::cout<<"setting collocation mat"<<std::endl;


    collocation_mat.setFromTriplets (eigen_triplets.begin(), eigen_triplets.end());

    //TODO: add toggle for computing svd and condition number
    if (model.COMPUTING_SVD_AND_CONDITION_NUMBER)
    {
        Eigen::JacobiSVD<MatrixXr> svd(collocation_mat);
        Real cond = svd.singularValues()(0)
        / svd.singularValues()(svd.singularValues().size()-1);
        std::cout<<"condition number: "<<cond<<std::endl;
        // collocation_mat.data().squeeze();
        std::cout << MatrixXr(collocation_mat) << std::endl;
    }
}


// Returns the opacity (takes into account the gaussian line profile)
inline Real Collocation :: get_opacity(Model& model, Size rayidx, Size freqidx, Size pointidx)
{
    Real toreturn=0;
    // if using other line profile/frequency basis function than gaussian, please replace the for(...) with the following line:
    // for (Size other_line_idx=0; other_line_idx<parameters.nlines(); other_line_idx++)
    // for (Size other_line_idx: other_frequency_basis_interacting_with[pointidx][lineidx])
    for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
    {
        Real delta_freq=0;
        if (USING_COMOVING_FRAME)
        {
            delta_freq=non_doppler_shifted_frequencies[freqidx]-model.lines.line[lineidx];
        }
        else
        {
            //local opacity should use the local frequencies, not the global ones...; otherwise it goes very quickly to zero if we have a minor doppler shift during the ray
            // I just don't feel like doppler shifting everything -> only very minor correction in the inv_width factor in front
            delta_freq=non_doppler_shifted_frequencies[freqidx]-model.lines.line[lineidx];
            // delta_freq=frequencies[rayidx][pointidx][freqidx]-model.lines.line[lineidx];
        }
        const Real inv_width=model.lines.inverse_width(pointidx,lineidx);
        toreturn+=INVERSE_SQRT_PI*inv_width*std::exp(-std::pow(delta_freq*inv_width,2))
                  *model.lines.line[lineidx]*model.lines.opacity(pointidx, lineidx);
    }
    return toreturn;
}

// Returns the emissivity (takes into account the gaussian line profile)
inline Real Collocation :: get_emissivity(Model& model, Size rayidx, Size freqidx, Size pointidx)
{
    Real toreturn=0;
    // if using other line profile/frequency basis function than gaussian, please replace the for(...) with the following line:
    // for (Size other_line_idx=0; other_line_idx<parameters.nlines(); other_line_idx++)
    // for (Size other_line_idx: other_frequency_basis_interacting_with[pointidx][lineidx])
    for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
    {
        Real delta_freq=0;
        if (USING_COMOVING_FRAME)
        {
            delta_freq=non_doppler_shifted_frequencies[freqidx]-model.lines.line[lineidx];
        }
        else
        {
            //local opacity should use the local frequencies, not the global ones...
            delta_freq=non_doppler_shifted_frequencies[freqidx]-model.lines.line[lineidx];
            // delta_freq=frequencies[rayidx][pointidx][freqidx]-model.lines.line[lineidx];
        }
        const Real inv_width=model.lines.inverse_width(pointidx,lineidx);
        toreturn+=INVERSE_SQRT_PI*inv_width*std::exp(-std::pow(delta_freq*inv_width,2))
                  *model.lines.line[lineidx]*model.lines.emissivity(pointidx, lineidx);
    }
    return toreturn;
}

inline void Collocation :: setup_rhs_Eigen(Model& model)
{
    rhs = VectorXr::Zero (collocation_mat.cols());
    for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
    {
        for (Size freqidx=0; freqidx<parameters.nfreqs(); freqidx++)
        {
            for (Size pointidx=0; pointidx<point_locations.size(); pointidx++)
            {
                const Size curr_mat_idx=get_mat_index(rayidx, freqidx, pointidx);
                //FIXME: check if boundary condition needs to apply
                Real local_velocity_gradient=0;
                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(pointidx))
                {
                    double dist=0;//these two variables are actually not used, but required in the function definition get_next
                    double ddist=0;

                    Size antipod=model.geometry.rays.antipod[rayidx];

                    const Size point_behind=model.geometry.get_next(pointidx, antipod, pointidx, dist, ddist);
                    const Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                    if (point_behind==parameters.npoints())
                    {//if there lies no point behind this point, we need to evaluate the boundary condition here
                        boundary_condition_required=true;

                        dist=0;//these two variables are actually not used, but required in the function definition get_next
                        ddist=0;
                        //also compute velocity gradient for the balancing factor
                        Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                        std::cout<<"next point: "<<next_point<<std::endl;
                        if (next_point!=parameters.npoints())//if there actually exists a next point
                        {
                            std::cout<<"computing velocity grad"<<std::endl;
                            local_velocity_gradient=model.geometry.rays.direction[rayidx].dot(model.geometry.points.velocity[next_point]-model.geometry.points.velocity[pointidx])/dist;
                        }
                    }
                }
                // boundary_condition_required=false;
                if (boundary_condition_required)
                {
                    // std::cout<<"local velocity grad: "<<local_velocity_gradient<<std::endl;
                    //This term makes sure that the magnitude of the matrix terms corresponding to this boundary condition is similar to the other nearby terms
                    const Real balancing_factor=compute_balancing_factor_boundary(rayidx, freqidx, pointidx, local_velocity_gradient, model);
                    rhs[curr_mat_idx]=balancing_factor*boundary_intensity(model, rayidx, freqidx, pointidx);
                }
                else
                {
                    //TODO: possibly some doppler shift required here? //no, we do not interpolate these values
                    rhs[curr_mat_idx]=get_emissivity(model, rayidx, freqidx, pointidx)/get_opacity(model, rayidx, freqidx, pointidx);
                }
            }
        }
    }
    std::cout<<"rhs: "<<rhs<<std::endl;

    rescale_matrix_and_rhs_Eigen(model);

    std::cout<<"rescaled rhs: "<<rhs<<std::endl;
}

// Rescale the matrix rows such that there are ones on the diagonal
// Practically, this is the same as applying a the inverse of the diagonal as left preconditioner
// TODO: the explicit multiplication factor for the boundary conditions might be a bit useless now
//Assumes that both the eigen matrix and rhs have already been constructed
inline void Collocation :: rescale_matrix_and_rhs_Eigen(Model& model)
{
    //get inverse of diagonal as diagonal matrix
    //FIXME: check whether we have zero on the diagonal; if yes, abort
    const auto diagonal_inverse=collocation_mat.diagonal().asDiagonal().inverse();
    rhs=diagonal_inverse*rhs;
    collocation_mat=diagonal_inverse*collocation_mat;

    //also printing out the rescaled version
    if (model.COMPUTING_SVD_AND_CONDITION_NUMBER)
    {
        Eigen::JacobiSVD<MatrixXr> svd(collocation_mat);
        Real cond = svd.singularValues()(0)
        / svd.singularValues()(svd.singularValues().size()-1);
        std::cout<<"condition number: "<<cond<<std::endl;
        // collocation_mat.data().squeeze();
        std::cout << MatrixXr(collocation_mat) << std::endl;
    }

}

///  Getter for the boundary conditions
///    @param[in] model  : reference to model object
///    @param[in] rayidx: ray index at which to evaluate the boundary condition
///    @param[in] freqidx: line frequency index at which to evaluate boundary condition
///    @param[in] pointidx: point index of the boundary point
///    @returns incoming radiation intensity at the boundary
////////////////////////////////////////////////////////////////////////////
///  Copied from solver.tpp and modified a bit
inline Real Collocation :: boundary_intensity (Model& model, Size rayidx, Size freqidx, Size pointidx)
{
    Real freq=0;
    if (USING_COMOVING_FRAME)
    {
        freq=non_doppler_shifted_frequencies[freqidx];
    }
    else
    {
        freq=frequencies[rayidx][pointidx][freqidx];
    }

    const Size bdy_id = model.geometry.boundary.point2boundary[pointidx];

    switch (model.geometry.boundary.boundary_condition[bdy_id])
    {
        case Zero    : return 0.0;
        case Thermal : return planck (model.geometry.boundary.boundary_temperature[bdy_id], freq);
        default      : return planck (T_CMB, freq);
    }
}


///  Planck function
///    @param[in] temp : temperature of the corresponding black body
///    @param[in] freq : frequency at which to evaluate the function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
///  Copied from solver.tpp
accel inline Real Collocation :: planck (Real temp, Real freq)
{
    return TWO_HH_OVER_CC_SQUARED * (freq*freq*freq) / expm1 (HH_OVER_KB*freq/temp);
}

// threaded_for (p, parameters.npoints(),
// {
//     for (Size l = 0; l < parameters.nlspecs(); l++)
//     {
//         for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
//         {
//             const Size lid = line_index (l, k);
//             emissivity (p, lid) = lineProducingSpecies[l].get_emissivity (p, k);
//                opacity (p, lid) = lineProducingSpecies[l].get_opacity    (p, k);
//         }
//     }
// })

inline void Collocation :: solve_collocation(Model& model)
{
    // std::cout<<"Eigen n threads: "<<Eigen::nbThreads()<<std::endl;
    // Eigen::BiCGSTAB<SparseMatrix<Real> > solver;
    // solver.compute(collocation_mat);
    // basis_coefficients = solver.solve(rhs);
    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error()      << std::endl;
    SparseLU <SparseMatrix<Real>, COLAMDOrdering<int>> solver;

    cout << "Analyzing system..."      << endl;

    solver.analyzePattern (collocation_mat);

    cout << "Factorizing system..."    << endl;

    solver.factorize (collocation_mat);

    if (solver.info() != Eigen::Success)
    {
        cout << "Factorization failed with error message:" << endl;
        cout << solver.lastErrorMessage()                  << endl;
        // cout << endl << RT << endl;

        throw std::runtime_error ("Eigen solver ERROR.");
    }

    cout << "Solving equations..." << endl;

    basis_coefficients = solver.solve(rhs);

    if (solver.info() != Eigen::Success)
    {
        cout << "Solving failed with error:" << endl;
        cout << solver.lastErrorMessage()    << endl;
        assert (false);
    }

    cout << "error: " << (collocation_mat*basis_coefficients-rhs).norm()/rhs.norm() << endl;

    // if (solver.info() != Eigen::Success)
    // {
    //     cout << "Solving failed with error:" << endl;
    //     cout << solver.lastErrorMessage()    << endl;
    //     assert (false);
    // }

    Size vectorsize=collocation_mat.cols();
    std::cout<<"computed coeffs"<<std::endl;
    for (Size i=0; i<vectorsize; i++)
    {
      std::cout<<"i: "<<i<<"value: "<<basis_coefficients(i)<<std::endl;
    }

}

// //TODO: calculate for all quadrature points the intensity; or exactly integrate the gaussian with another gaussian
inline void Collocation :: compute_J(Model& model)
{
  model.radiation.initialize_J();

  // //For testing purposes, initialize I too
  // for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
  // {
  //     for (Size pointidx=0; pointidx<parameters.npoints(); pointidx++)
  //     {
  //         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
  //         {
  //             model.radiation.I(rayidx, pointidx, lineidx)=0.0;
  //         }
  //     }
  // }
  //
  // for (Size rayidx=0; rayidx<parameters.hnrays(); rayidx++)
  // {
  //     for (Size pointidx=0; pointidx<parameters.npoints(); pointidx++)
  //     {
  //         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
  //         {
  //             model.radiation.u(rayidx, pointidx, lineidx)=0.0;
  //         }
  //     }
  // }

  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
      for (Size p=0; p<parameters.npoints(); p++)
      // threaded_for (p, parameters.npoints(),
      {
          // const Size p=points_in_grid[idx];

          for (Size k = 0; k < lspec.linedata.nrad; k++)
          {
              // const Size1 freq_nrs = lspec.nr_line[p][k];
              const Size lid = model.lines.line_index(l, k);
              const Size1 freq_nrs = lspec.nr_line[p][k];
              const Size middle_linefreqidx=freq_nrs[parameters.nquads()/2];//for determining which frequencies interact with the line frequency; FIXME: remove this hack by explicitly calculating this

              // Initialize values
              lspec.Jlin[p][k] = 0.0;

              // // Just sum for all directions over the non-zero basis functions
              // std::set<std::vector<Size>> nonzero_triplets;
              // get_nonzero_basis_triplets(nonzero_triplets, rayidx, lid, pointidx);


              //or just use a matrix-vector product
              for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
              {
                  Real lp_freq=0;
                  Real lp_freq_inv_width=0;
                  if (USING_COMOVING_FRAME)
                  {
                      // lp_freq=non_doppler_shifted_frequencies[lid];
                      // lp_freq_inv_width=non_doppler_shifted_inverse_widths[lid][p];
                      lp_freq=model.lines.line[lid];
                      lp_freq_inv_width=model.lines.inverse_width(lid,p);
                  }
                  else
                  {
                      const Real doppler_shift_factor=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[p],model.geometry.rays.direction[rayidx]);
                      std::cout<<"doppler shift factor: "<<doppler_shift_factor<<std::endl;
                      //FIXME (or delete the non-comoving version): add doppler shifts to line frequencies (and inverse widths)
                      // lp_freq=frequencies[rayidx][p][lid];
                      // lp_freq_inv_width=frequencies_inverse_widths[rayidx][p][lid];
                      lp_freq=model.lines.line[lid]*doppler_shift_factor;
                      lp_freq_inv_width=model.lines.inverse_width(lid,p)/doppler_shift_factor;
                  }
                  // Size hrayidx=rayidx;//rayidx for u
                  // if (hrayidx>=parameters.hnrays())
                  // {
                  //     hrayidx=model.geometry.rays.antipod[hrayidx];
                  // }

                  std::set<std::vector<Size>> nonzero_triplets;

                  get_nonzero_basis_triplets(nonzero_triplets, rayidx, middle_linefreqidx, p);
                  Real I=0;
                  Real temp_jlin=0;
                  for (std::vector<Size> triplet: nonzero_triplets)
                  // for (Size z = 0; z < parameters.nquads(); z++)
                  {
                      const Size triplet_basis_index=get_mat_index(triplet[0], triplet[1], triplet[2]);
                      if (USING_COMOVING_FRAME)
                      {
                          // const Real doppler_shifted_frequency=lp_freq*calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[p]
                          //                                                                                       , model.geometry.rays.direction[rayidx]);
                          // Real doppler_shifted_frequency=lp_freq;//try it out non-doppler shifted

                          I+=basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          lspec.Jlin[p][k]+=basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                          // I+=basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency);
                          // lspec.Jlin[p][k]+=basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency, lp_freq_inv_width);
                      }
                      else
                      {
                          //TODO: figure out which values we actually want to keep...; note to self: I and u use the frequency index (nlines*nquads); so we are not able to compute it here
                          // model.radiation.I(rayidx, p, lid) +=basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          // model.radiation.u(hrayidx, p, lid)+=1.0/2*basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          //FIXME: doppler shift the line profile frequency such that it corresponds to the local frequency!!!!!!!!!

                          I+=basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          lspec.Jlin[p][k]+=basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                          // std::cout<<"please be nonzero: "<<basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[2], point_locations[p])*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width)<<std::endl;

                          // lspec.Jlin[p][k]+=200000*basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p])*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                      }

                  }
                  std::cout<<"p: "<<p<<" I: "<<I<<std::endl;

              }

              // // Integrate over the line
              // for (Size z = 0; z < parameters.nquads(); z++)
              // {
              //     lspec.Jlin[p][k] += lspec.quadrature.weights[z] * radiation.J(p, freq_nrs[z]);
              // }


              double diff = 0.0;

              // // Collect the approximated part
              // for (Size m = 0; m < lspec.lambda.get_size(p,k); m++)
              // {
              //     const Size I = lspec.index(lspec.lambda.get_nr(p,k,m), lspec.linedata.irad[k]);
              //
              //     diff += lspec.lambda.get_Ls(p,k,m) * lspec.population[I];
              // }

              lspec.Jeff[p][k] = lspec.Jlin[p][k] - HH_OVER_FOUR_PI * diff;
              lspec.Jdif[p][k] = HH_OVER_FOUR_PI * diff;

              // if (lspec.Jeff[p][k]<0)
              // {
              //     lspec.Jeff[p][k]=0;
              // }
              // std::cout<<"Jeff: "<<lspec.Jeff[p][k]<<std::endl;
          }
      }//)
  }
  //print out j for first transition
  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
      for (Size pointidx=0; pointidx<parameters.npoints(); pointidx++)
      {
          std::cout<<"jeff: "<<lspec.Jeff[pointidx][0]<<std::endl;
      }
  }
}
