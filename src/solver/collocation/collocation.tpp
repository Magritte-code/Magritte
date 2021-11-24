#include <cmath>


inline void Collocation :: setup_basis(Model& model)
{
    //first, we need all frequencies with their respective gaussians
    //TODO: ask frederik how these are stored (p,f)->individually for each point?
    frequencies.resize(parameters.nfreqs());
    frequencies_widths.resize(parameters.nfreqs());
    //then we also set all the radii and locations for the radial basis functions
    //currently just uses all the points // might need to be changed later to all the points in a certain grid // note: currently also adding the boundary points
    //Size n_rbf_basisfuns=TODO;
    point_locations.resize(parameters.npoints());
    rbf_radii.resize(parameters.npoints());
    other_basis_interacting_with.resize(parameters.npoints());
    // for (Size i:model.geometry.points.multiscale.get_current_points_in_grid())
    for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    {
        point_locations[pointid]=model.geometry.points.position[pointid];
        //rbf_radii[pointid]=TODO just use the nanoflann package for determining the radii (get N+1 closest points)
        //other_basis_interacting_with[pointid]=TODO also use it for the actual interacting indices
    }
}



/// The basis function associated with each direction
inline Real Collocation :: basis_direction(Size rayindex)
{
    //as we do not need directional derivatives, piecewise constant basis functions are sufficient;
    return 1;
}

///  The basis function associated with each frequency, evaluated at frequency currnu
inline Real Collocation :: basis_freq(Size freqidx, Real currfreq)
{
    //TODO: figure out what to use exactly
    //also implement some logical data structure (possibly first read out all the sorted frequencies, then construct the basis functions)
    Real abs_freq_diff=std::abs(currfreq-frequencies[freqidx]);
    if (abs_freq_diff>frequencies_widths[freqidx]*TRUNCATION_SIGMA)
    {
        return 0;
    }
    else
    {
        return std::exp(-std::pow(abs_freq_diff/frequencies_widths[freqidx],2));
    }
}
///  The derivative of the frequency basis function
inline Real Collocation :: basis_freq_der(Size freqidx, Real currfreq)
{
    Real abs_freq_diff=std::abs(currfreq-frequencies[freqidx]);
    if (abs_freq_diff>frequencies_widths[freqidx]*TRUNCATION_SIGMA)
    {
        return 0;
    }
    else
    {
        return -2*(abs_freq_diff/frequencies_widths[freqidx])*std::exp(-std::pow(abs_freq_diff/frequencies_widths[freqidx],2));
    }
}

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
