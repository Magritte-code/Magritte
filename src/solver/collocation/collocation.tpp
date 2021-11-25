#include <cmath>

inline void Collocation :: set_interacting_bases(Model& model)
{
  //Set the kd_tree
  Eigen::Matrix<double, Eigen::Dynamic, 3> temp_mat;
  Size1 index_conversion;
  auto tuple=model.create_mat_for_kd_tree_of_lvl(0);
  temp_mat=std::get<0>(tuple);
  index_conversion=std::get<1>(tuple);
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
        other_radial_basis_interacting_with[pointid]=subvector;//this also includes itself

        for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        {
            other_frequency_basis_interacting_with[rayid][pointid].resize(subvector.size());
            for (Size temp=0; temp<subvector.size(); temp++)
            {
                other_frequency_basis_interacting_with[rayid][pointid][temp].resize(parameters.nfreqs());
            }
        }
        // also set what interactions there are between the frequency basis functions, which unfortunately also depend on the interacting points
        Size reduced_index=0;
        for (Size other_point: subvector)
        {
            for (Size freqid=0; freqid<parameters.nfreqs(); freqid++)
            {
                for (Size rayid=0; rayid<parameters.nrays(); rayid++)
                {
                    Real curr_freq=frequencies[rayid][pointid][freqid];
                    vector<Real> searchfreqs=frequencies[rayid][other_point];
                    for (Size search_freq_idx=0; search_freq_idx<parameters.nfreqs(); search_freq_idx++)
                    {   // The frequency basis function at the other point and search frequency is only nonzero for the frequency frequencies[pointid][freqid] if the main frequency of that gaussian does not lie to far from it (compared to the width)
                        Real abs_freq_diff=std::abs(curr_freq-searchfreqs[search_freq_idx]);
                        if (!freq_cutoff_condition(abs_freq_diff,frequencies_widths[rayid][other_point][search_freq_idx]))
                        {
                            other_frequency_basis_interacting_with[rayid][pointid][reduced_index][freqid].insert(search_freq_idx);
                        }
                    }
                }
            }
            reduced_index++;
        }
    }

}



inline void Collocation :: setup_basis(Model& model)
{
    //first, we need all frequencies with their respective gaussians
    //TODO: ask frederik how these are stored (p,f)->individually for each point?

    //then we also set all the radii and locations for the radial basis functions
    //currently just uses all the points // might need to be changed later to all the points in a certain grid // note: currently also adding the boundary points
    point_locations.resize(parameters.npoints());
    rbf_radii.resize(parameters.npoints());
    other_radial_basis_interacting_with.resize(parameters.npoints());
    //We need all frequecies at each point (can be different at each point due to doppler shift)
    frequencies.resize(parameters.nrays());
    frequencies_widths.resize(parameters.nrays());
    other_frequency_basis_interacting_with.resize(parameters.nrays());  //?with a bit more computational effort, this could be reduced a bit in size
    // for (Size i:model.geometry.points.multiscale.get_current_points_in_grid())
    for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    {
        point_locations[pointid]=model.geometry.points.position[pointid];
    }

    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        frequencies[rayid].resize(parameters.npoints());
        frequencies_widths[rayid].resize(parameters.npoints());
        other_frequency_basis_interacting_with[rayid].resize(parameters.npoints());

        for (Size pointid=0; pointid<parameters.npoints(); pointid++)
        {
            frequencies[rayid][pointid].resize(parameters.nfreqs());
            frequencies_widths[rayid][pointid].resize(parameters.nfreqs());
        }
    }



    set_interacting_bases(model);
    //Now determine which basis's interact with eachother; Note to self: frequency basis might depend on the exact position in the future




}


inline bool Collocation :: freq_cutoff_condition(Real abs_freq_diff, Real freq_basis_width)
{
    return (abs_freq_diff>freq_basis_width*TRUNCATION_SIGMA);
}

/// The basis function associated with each direction
inline Real Collocation :: basis_direction(Size rayindex)
{
    //as we do not need directional derivatives, piecewise constant basis functions are sufficient;
    return 1;
}

///  The basis function associated with each frequency, evaluated at frequency currnu
inline Real Collocation :: basis_freq(Size rayidx, Size freqidx, Size pointidx, Real currfreq)
{
    //TODO: figure out what to use exactly
    //also implement some logical data structure (possibly first read out all the sorted frequencies, then construct the basis functions)
    Real abs_freq_diff=std::abs(currfreq-frequencies[rayidx][pointidx][freqidx]);
    if (freq_cutoff_condition(abs_freq_diff,frequencies_widths[rayidx][pointidx][freqidx]))
    {
        return 0;
    }
    else
    {
        return std::exp(-std::pow(abs_freq_diff/frequencies_widths[rayidx][freqidx][pointidx],2));
    }
}

// ///  The derivative of the frequency basis function
// inline Real Collocation :: basis_freq_der(Size freqidx, Real currfreq)
// {
//     Real abs_freq_diff=std::abs(currfreq-frequencies[freqidx]);
//     if (freq_cutoff_condition(abs_freq_diff,frequencies_widths[freqidx]))
//     {
//         return 0;
//     }
//     else
//     {
//         return -2*(abs_freq_diff/frequencies_widths[freqidx])*std::exp(-std::pow(abs_freq_diff/frequencies_widths[freqidx],2));
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
inline void Collocation :: get_nonzero_basis_triplets(std::set<std::vector<Size>>& basis_triplets_to_fill, Size rayidx, Size freqidx, Size pointidx)
{
    // std::set<std::vector<Size>> basis_triplets;
    // I currently assume that the different ray direction bases do not interact with eachother
    for (Size interacting_points: other_radial_basis_interacting_with[pointidx])
    {
        for (Size interacting_freqs: other_frequency_basis_interacting_with[rayidx][pointidx][interacting_points][freqidx])
        {
            vector<Size> temp_triplet{rayidx, interacting_freqs, interacting_points};
            basis_triplets_to_fill.insert(temp_triplet);
        }
    }
    // return basis_triplets;
}

///  Returns the matrix index corresponding to the general index (rayidx, freqidx, pointidx)
inline Size Collocation :: get_mat_index(Size rayidx, Size freqidx, Size pointidx)
{
    return rayidx*frequencies.size()*point_locations.size()+freqidx*point_locations.size()+pointidx;
}

// inline Real Collocation :: calculate_doppler_shift(Size pointidx, Size rayidx, Geometry& geometry)
// {
//     Vector3D raydirection=model.geometry.rays.direction[rayindex];
//     return 1+std::sqrt(raydirection.dot())
// }

inline void Collocation :: setup_basis_matrix_Eigen(Model& model)
{
    Size mat_size=parameters.nrays()*frequencies.size()*point_locations.size();

    for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
    {
        for (Size pointidx=0; pointidx<point_locations.size(); pointidx++)
        {//note: frequency and point ordering switched because the frequency widths might depend on the points

            //TODO: either store all doppler shifted frequencies or calculate some relative doppler shifts
            //The doppler shift term n.grad(n.v(x)/c)*nu does not depend on the basis functions. Thus we calculate it here
            //TODO: actually check how it's calculated in Geometry
            // Real doppler_shift_term=0;//temporary value FIXME


            for (Size freqidx=0; freqidx<frequencies.size(); freqidx++)
            {
                //Now creating all the triplets, thus we need all nonzero basis triplets
                std::set<std::vector<Size>> nonzero_triplets;
                get_nonzero_basis_triplets(nonzero_triplets, rayidx, freqidx, pointidx);
                for (std::vector<Size> triplet: nonzero_triplets)
                {
                    const Size mat_idx=get_mat_index(triplet[0],triplet[1],triplet[2]);

                    Real curr_freq=frequencies[rayidx][pointidx][freqidx];// TODO: add doppler shift
                    Vector3D curr_location=point_locations[pointidx];

                    //TODO create EIGEN triple TODO
                    //start with the directional derivative
                    Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry);

                    // Vector3D velocity_gradient=TODO;
                    //also add doppler shift term //TODO: calculate doppler shift

                    //For every scattering direction! add scattering//currently just assume all directions
                    for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                    {
                        //TODO: add toggle to whether we are actually using scattering
                        Real integral_scatt_redistr=0;//FIXME: add me
                        Real scattering=-basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                    }
                }


            }
        }
    }

}
