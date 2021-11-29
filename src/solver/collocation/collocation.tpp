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
        Vector3D temp_vector=(model.geometry.points.position[neighbors_coarser_grid.back()]-model.geometry.points.position[pointid]);
        rbf_radii[pointid]=std::sqrt(temp_vector.squaredNorm());
        //currently, i am choosing the radius such that the last point lies on the boundary of the domain of the radial basis function
        //Therefore, we do not need to include it in the vector of interacting basis functions
        vector<Size> subvector={neighbors_coarser_grid.begin(),neighbors_coarser_grid.end()-1};
        std::cout<<"subvector size should be 3/15 and is: "<<subvector.size()<<std::endl;


        //calculate the max doppler shift possible
        Real max_doppler_shift=1.0;
        // vector<Real> min_doppler_shift(parameters.nrays(),1.0);
        Vector3D velocity_curr_point=model.geometry.points.velocity[pointid];
        // for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
        // {
        for (Size neighbor: subvector)
        {
            max_doppler_shift=std::max(max_doppler_shift,Real(
        1+std::sqrt((model.geometry.points.velocity[neighbor]-velocity_curr_point).squaredNorm())));
            // min_doppler_shift[rayidx]=std::min(max_doppler_shift,
        // 1+std::sqrt((model.geometry.points.velocity[neighbor]-velocity_curr_point).dot(model.geometry.rays.direction[rayidx])/CC);
        }
        // }

        //Now converting the subvector to the right (basis) indices
        for (Size index=0; index<subvector.size(); index++)
        {
            subvector[index]=reverse_index_conversion[subvector[index]];
        }
        other_radial_basis_interacting_with[pointid]=subvector;//this also includes itself

        // For every frequency, calculate which frequencies lie nearby
        for (Size lineid=0; lineid<parameters.nlines(); lineid++)
        {
            //The line frequencies are not ordered, so just search the whole list...
            for (Size other_lineid; other_lineid<parameters.nlines(); other_lineid++)
            {
                // Somewhat loose uniform (over direction and neighbors) estimate
                // Start by bounding the apparent frequency nu' for the point (pointidx, lineid): nu/max_dpplr_shift<nu'<nu*max_dpplr_shift
                if ((model.lines.line[lineid]/max_doppler_shift<model.lines.line[other_lineid]+1/model.lines.inverse_width(pointid,other_lineid)*TRUNCATION_SIGMA*max_doppler_shift)
                &&  (model.lines.line[lineid]*max_doppler_shift>model.lines.line[other_lineid]-1/model.lines.inverse_width(pointid,other_lineid)*TRUNCATION_SIGMA*max_doppler_shift))
                {
                    other_frequency_basis_interacting_with[pointid][lineid].insert(other_lineid);
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


inline void Collocation :: setup_basis(Model& model)
{

    //then we also set all the radii and locations for the radial basis functions
    //currently just uses all the points // might need to be changed later to all the points in a certain grid // note: currently also adding the boundary points
    point_locations.resize(parameters.npoints());
    rbf_radii.resize(parameters.npoints());
    other_radial_basis_interacting_with.resize(parameters.npoints());
    //We need all frequecies at each point (can be different at each point due to doppler shift)
    frequencies.resize(parameters.nrays());
    frequencies_inverse_widths.resize(parameters.nrays());
    other_frequency_basis_interacting_with.resize(parameters.npoints());  //?with a bit more computational effort, this could be reduced a bit in size
    // for (Size i:model.geometry.points.multiscale.get_current_points_in_grid())
    for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    {
        point_locations[pointid]=model.geometry.points.position[pointid];
        other_frequency_basis_interacting_with[pointid].resize(parameters.nlines());
    }

    // We need the inverse line widths, so we might as well make sure that they are calculated now
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
            frequencies[rayid][pointid].resize(parameters.nlines());
            frequencies_inverse_widths[rayid][pointid].resize(parameters.nlines());

            for (Size lineid=0; lineid<parameters.nlines(); lineid++)
            {
                frequencies[rayid][pointid][lineid]=doppler_shift*model.lines.line[lineid];
                frequencies_inverse_widths[rayid][pointid][lineid]=model.lines.inverse_width(pointid,lineid)/doppler_shift;

            }

        }
    }

    // Also do not forget to initialize the basis coefficients // TODO: figure out which exact format (std::vector, Eigen, ...) to use
    basis_coefficients=vector<Real>(0.0,parameters.nrays()*parameters.nlines()*parameters.npoints());

    set_interacting_bases(model);
    //Now determine which basis's interact with eachother




}


inline bool Collocation :: freq_cutoff_condition(Real abs_freq_diff, Real freq_basis_inverse_width)
{
    return (abs_freq_diff*freq_basis_inverse_width>TRUNCATION_SIGMA);
}

/// The basis function associated with each direction
inline Real Collocation :: basis_direction(Size rayindex)
{
    //as we do not need directional derivatives, piecewise constant basis functions are sufficient;
    return 1;
}

///  The basis function associated with each frequency, evaluated at frequency currnu
inline Real Collocation :: basis_freq(Size rayidx, Size lineidx, Size pointidx, Real currfreq)
{
    //TODO: figure out what to use exactly
    //also implement some logical data structure (possibly first read out all the sorted frequencies, then construct the basis functions)
    Real abs_freq_diff=std::abs(currfreq-frequencies[rayidx][pointidx][lineidx]);
    if (freq_cutoff_condition(abs_freq_diff,frequencies_inverse_widths[rayidx][pointidx][lineidx]))
    {
        return 0;
    }
    else
    {
        return INVERSE_SQRT_PI*frequencies_inverse_widths[rayidx][lineidx][pointidx]*std::exp(-std::pow(abs_freq_diff*frequencies_inverse_widths[rayidx][lineidx][pointidx],2));
    }
}

// ///  The derivative of the frequency basis function
// inline Real Collocation :: basis_freq_der(Size lineidx, Real currfreq)
// {
//     Real abs_freq_diff=std::abs(currfreq-frequencies[lineidx]);
//     if (freq_cutoff_condition(abs_freq_diff,frequencies_widths[lineidx]))
//     {
//         return 0;
//     }
//     else
//     {
//         return -2*(abs_freq_diff/frequencies_widths[lineidx])*std::exp(-std::pow(abs_freq_diff/frequencies_widths[lineidx],2));
//     }
// }

///  The radial basis function for the position
inline Real Collocation :: basis_point(Size centerpoint, Vector3D& location)
{
    Vector3D diff_vector=location-point_locations[centerpoint];
    Real distance=std::sqrt(diff_vector.dot(diff_vector));
    Real radius=rbf_radii[centerpoint];
    //if the location is too far from the center, the evalutation of the compact radial basis function is 0.
    if (distance>=radius)
    {
        return 0;
    }
    else
    {
        Real rel_dist=distance/radius;
        //This basis function has nice continuous derivatives at r=0 and r=1 (being 0 at both places);
        return -4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
    }
}

///  The directional derivative of the position basis function
inline Real Collocation :: basis_point_der(Size centerpoint, Vector3D& location, Size rayindex, Geometry& geometry)
{
    Vector3D diff_vector=location-point_locations[centerpoint];
    Real distance=std::sqrt(diff_vector.dot(diff_vector));
    Real radius=rbf_radii[centerpoint];
    //if the location is too far from the center, the evalutation of the compact radial basis function is 0.
    if (distance>=radius)
    {
        return 0;
    }
    else
    {
        Real rel_dist=distance/radius;

        Vector3D raydirection=geometry.rays.direction[rayindex];
        //The derivatives in the perpendicular direction (to the radial direction) are 0, so we can calculate the derivative
        // by taking the radial derivative and multiplying it by cos(theta) in which theta is the angle between the radial and the actual direction
        //Below, we see the simple rule for calculating the cosine
        //Note: The ray directions are normalised, so we can omit the normalization of the ray direction
        Real costheta=raydirection.dot(diff_vector)/std::sqrt(diff_vector.squaredNorm());

        Real radial_derivative =  12*std::pow(rel_dist,2)/std::pow(1+std::pow(rel_dist,3),2)-12*rel_dist/std::pow((1+std::pow(rel_dist,2)),2);

        return costheta*radial_derivative;
    }
}

/// Fills the triplets which correspond to the nonzero basis function products at a given triplet location
///   @param[in/out] basis_triplets_to_fill: An empty set which will contain the basis triplets after execution
inline void Collocation :: get_nonzero_basis_triplets(std::set<std::vector<Size>>& basis_triplets_to_fill, Size rayidx, Size lineidx, Size pointidx)
{
    // std::set<std::vector<Size>> basis_triplets;
    // I currently assume that the different ray direction bases do not interact with eachother
    for (Size interacting_points: other_radial_basis_interacting_with[pointidx])
    {
        for (Size interacting_freqs: other_frequency_basis_interacting_with[pointidx][lineidx])
        {
            vector<Size> temp_triplet{rayidx, interacting_freqs, interacting_points};
            basis_triplets_to_fill.insert(temp_triplet);
        }
    }
    // return basis_triplets;
}

///  Returns the matrix index corresponding to the general index (rayidx, lineidx, pointidx)
inline Size Collocation :: get_mat_index(Size rayidx, Size lineidx, Size pointidx)
{
    return rayidx*frequencies.size()*point_locations.size()+lineidx*point_locations.size()+pointidx;
}

// inline Real Collocation :: calculate_doppler_shift(Size pointidx, Size rayidx, Geometry& geometry)
// {
//     Vector3D raydirection=model.geometry.rays.direction[rayindex];
//     return 1+std::sqrt(raydirection.dot())
// }

inline void Collocation :: setup_basis_matrix_Eigen(Model& model)
{
    Size mat_size=parameters.nrays()*frequencies.size()*point_locations.size();
    vector<Triplet<Real, Size>> eigen_triplets;
    // eigen_triplets.reserve(mat_size*16*nrays());//matrix size times number of interacting points times the number of rays (scattering)
    eigen_triplets.reserve(mat_size*16);//matrix size times number of interacting points (does not include scattering)

    for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
    {
        for (Size lineidx=0; lineidx<frequencies.size(); lineidx++)
        {


            //TODO: either store all doppler shifted frequencies or calculate some relative doppler shifts
            //The doppler shift term n.grad(n.v(x)/c)*nu does not depend on the basis functions. Thus we calculate it here
            //TODO: actually check how it's calculated in Geometry
            // Real doppler_shift_term=0;//temporary value FIXME

            for (Size pointidx=0; pointidx<point_locations.size(); pointidx++)
            {
                //Now creating all the triplets, thus we need all nonzero basis triplets
                std::set<std::vector<Size>> nonzero_triplets;
                get_nonzero_basis_triplets(nonzero_triplets, rayidx, lineidx, pointidx);

                const Size curr_mat_idx=get_mat_index(rayidx, lineidx, pointidx);

                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(pointidx))
                {
                    double dist=0;//these two variables are actually not used, but required in the function definition get_next
                    double ddist=0;

                    const Size point_behind=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                    if (point_behind==parameters.npoints())
                    {//if there lies no point behind this point, we need to evaluate the boundary condition here
                        boundary_condition_required=true;
                    }
                }
                if (boundary_condition_required)
                {
                    //Merely evaluating the intensity at the boundary
                    for (std::vector<Size> triplet: nonzero_triplets)
                    {
                        const Size other_mat_idx=get_mat_index(triplet[0],triplet[1],triplet[2]);

                        Real curr_freq=frequencies[rayidx][pointidx][lineidx];// TODO: add doppler shift
                        Vector3D curr_location=point_locations[pointidx];

                        //TODO create EIGEN triple TODO
                        //start with the directional derivative
                        Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                        //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                        eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, other_mat_idx, basis_eval));

                    }
                }
                else
                {
                    for (std::vector<Size> triplet: nonzero_triplets)
                    {
                        const Size other_mat_idx=get_mat_index(triplet[0],triplet[1],triplet[2]);

                        Real curr_freq=frequencies[rayidx][pointidx][lineidx];
                        Vector3D curr_location=point_locations[pointidx];

                        //Start with the directional derivative
                        Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry);
                        //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                        Real opacity=get_opacity(model, rayidx, pointidx, lineidx)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                        eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, other_mat_idx, directional_der+opacity));

                        //For every scattering direction! add scattering
                        for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                        {
                            const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                            //TODO: add toggle to whether we are actually using scattering
                            Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                            Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                            //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                            eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                        }

                    }
                }
            }
        }
    }


    //after all this, construct the matrix
    collocation_mat.resize (mat_size, mat_size);

    collocation_mat.setFromTriplets (eigen_triplets.begin(), eigen_triplets.end());
    collocation_mat.data().squeeze();
}


// Returns the opacity (takes into account the gaussian line profile)
inline Real Collocation :: get_opacity(Model& model, Size rayidx, Size lineidx, Size pointidx)
{
    Real toreturn=0;
    // if using other line profile/frequency basis function than gaussian, please replace the for(...) with the following line:
    // for (Size other_line_idx=0; other_line_idx<parameters.nlines(); other_line_idx++)
    for (Size other_line_idx: other_frequency_basis_interacting_with[pointidx][lineidx])
    {
        const Real delta_freq=frequencies[rayidx][pointidx][lineidx]-frequencies[rayidx][pointidx][other_line_idx];
        const Real inv_width=frequencies_inverse_widths[rayidx][other_line_idx][pointidx];
        toreturn+=INVERSE_SQRT_PI*inv_width*std::exp(-std::pow(delta_freq*inv_width,2))
                  *frequencies[rayidx][pointidx][other_line_idx]*model.lines.opacity(pointidx, other_line_idx);
    }
    return toreturn;
}

// Returns the emissivity (takes into account the gaussian line profile)
inline Real Collocation :: get_emissivity(Model& model, Size rayidx, Size lineidx, Size pointidx)
{
    Real toreturn=0;
    // if using other line profile/frequency basis function than gaussian, please replace the for(...) with the following line:
    // for (Size other_line_idx=0; other_line_idx<parameters.nlines(); other_line_idx++)
    for (Size other_line_idx: other_frequency_basis_interacting_with[pointidx][lineidx])
    {
        const Real delta_freq=frequencies[rayidx][pointidx][other_line_idx]-frequencies[rayidx][pointidx][lineidx];
        const Real inv_width=frequencies_inverse_widths[rayidx][other_line_idx][pointidx];
        toreturn+=INVERSE_SQRT_PI*inv_width*std::exp(-std::pow(delta_freq*inv_width,2))
                  *frequencies[rayidx][pointidx][other_line_idx]*model.lines.emissivity(pointidx, other_line_idx);
    }
    return toreturn;
}

inline void Collocation :: setup_rhs_Eigen(Model& model)
{
    rhs = VectorXr::Zero (collocation_mat.cols());
    for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
    {
        for (Size lineidx=0; lineidx<frequencies.size(); lineidx++)
        {
            for (Size pointidx=0; pointidx<point_locations.size(); pointidx++)
            {
                const Size curr_mat_idx=get_mat_index(rayidx, lineidx, pointidx);
                //FIXME: check if boundary condition needs to apply
                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(pointidx))
                {
                    double dist=0;//these two variables are actually not used, but required in the function definition get_next
                    double ddist=0;

                    const Size point_behind=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                    if (point_behind==parameters.npoints())
                    {//if there lies no point behind this point, we need to evaluate the boundary condition here
                        boundary_condition_required=true;
                    }
                }
                if (boundary_condition_required)
                {
                    //filling in boundary intensity
                    rhs[curr_mat_idx]=0;
                    // rhs[curr_mat_idx]=boundary_intensity(model, pointidx, frequencies[rayidx][pointidx][lineidx]);
                }
                else
                {
                    rhs[curr_mat_idx]=get_emissivity(model, rayidx, lineidx, pointidx);
                }
            }
        }
    }
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

// //TODO: calculate for all quadrature points the intensity; or exactly integrate the gaussian with another gaussian
// inline void Collocation :: compute_J(Model& model)
// {
//   model.radiation.initialize_J();
//   for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
//   {
//       threaded_for (p, parameters.npoints(),
//       {
//           // const Size p=points_in_grid[idx];
//
//           for (Size k = 0; k < lspec.linedata.nrad; k++)
//           {
//               const Size1 freq_nrs = lspec.nr_line[p][k];
//
//               // Initialize values
//               lspec.Jlin[p][k] = 0.0;
//
//               // Just sum for all directions over the non-zero basis functions
//
//
//               //or just use a matrix-vector product
//               for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
//               {
//                   std::set<std::vector<Size>> nonzero_triplets;
//                   get_nonzero_basis_triplets(nonzero_triplets, rayidx, lineidx, p);
//                   lspec.Jlin[p][k]
//               }
//
//               // // Integrate over the line
//               // for (Size z = 0; z < parameters.nquads(); z++)
//               // {
//               //     lspec.Jlin[p][k] += lspec.quadrature.weights[z] * radiation.J(p, freq_nrs[z]);
//               // }
//
//
//               double diff = 0.0;
//
//               // Collect the approximated part
//               for (Size m = 0; m < lspec.lambda.get_size(p,k); m++)
//               {
//                   const Size I = lspec.index(lspec.lambda.get_nr(p,k,m), lspec.linedata.irad[k]);
//
//                   diff += lspec.lambda.get_Ls(p,k,m) * lspec.population[I];
//               }
//
//               lspec.Jeff[p][k] = lspec.Jlin[p][k] - HH_OVER_FOUR_PI * diff;
//               lspec.Jdif[p][k] = HH_OVER_FOUR_PI * diff;
//           }
//       })
//   }
// }
