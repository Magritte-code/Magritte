#pragma once

#include "model/model.hpp"
#include "model/parameters/parameters.hpp"
#include <set>
#include <map>
#include <Eigen/Core>

// Helper class for setting up the collocation solver
class Collocation
{

private:
    // Comoving frame formulation not necessary and complicates equations; therefore to be deprecated
    // const bool USING_COMOVING_FRAME=true;//toggle for using the comoving frame
    const bool USING_COMOVING_FRAME=false;
    const bool USING_LEAST_SQUARES=false;
    const bool PRUNING_POSITIONAL_BASIS_FUNCTIONS=true;

    //internal parameters for the basis functions
    const Real TRUNCATION_SIGMA=5.0; //< if Delta(freq)>trunc_sigma*freq_width just set the result from the gaussian basis function to zero (as it is almost zero)
    const Real SLOPE_FACTOR=1.0;//1.0;//For stabilizing the equations a bit (making them somewhat more diagonally dominant);
    const Real SLOPE_STEEPNESS=4.0;//4.0;//slope steepness, also for stabilizing the equations; will practically be divided by 4 at x=0.
    const Real FREQ_OFFSET=0.0;//-0.2;//offset for the frequency basis functions (in terms of the inverse line width (mind the sqrt(2)))-in comoving frame
    //FIXME: freq offset makes the freq cutoff condition a bit nonsymmetrical

    //TODO: remove slope, as it can cause unphysical limit behaviour due to numerical inaccuracies
    const Real SLOPE_THRESHOLD=7.0*2.0;//The threshold for determining whether we need to add a slope to the radial basis function
    //should be related to the max derivative of the rbf, the number of neighbors (for diagonal dominance reasons). Note: this comes from a global bound
    //if opacity*typical radius>SLOPE_THRESHOLD, then do not add the slope.
    // const Real QUADRATURE_MIN=-0.5;
    // const Real QUADRATURE_MAX=0.5;
    // TODO: define quadratures such that we can probe the entire space... abs max doppler shift-1 should not exceed the quadrature bounds
    const Real MIN_OPACITY=1.0*std::pow(10,-28);//FIXME: use instead some adaptive measure based on rbf_der/(opacity*rbf). of diagonal element

    vector<vector<Size>> remaining_points_after_pruning;
    // Size n_freq_bases=1;//denotes the number of frequency bases used
    // // not necessary, as the size of frequency_quadrature will tell us the same
    // // = number of evaluation frequencies

    vector<vector<vector<bool>>> rbf_using_slope; //For every ray, for every point, for every frequency, denotes whether radial basis functions at a given point need to include the slope
    //currently using some global bound; MAYBE TODO: use the exact formulation for each basis function; this will need some code duplication (explicit rbf derivative without slope), and will obviously also depend on the exact direction (hatn dot grad)
    //frequency dependence is due to multiplying with the opacity and angular dependence might come into play when the rbf basis function could depend on the angle

    // DEPRECATED: makes that the basis functions are too peaked for my liking
    // const Real DISTANCE_EXPONENT=1.0/5.0; //Exponent for redistributing the distances
    // const Real DISTANCE_EPS=0.1;//TODO: should actually be chosen adaptively (depending on how close the closest point lies)
    // vector<Real> half_min_dist; //Half the distance for

    const Real MAX_FREQ_QUAD_RATIO=0.7;//gives the maximum ratio between evaluation in the neighboring freq basis versus the evaluation in local frequency basis
    const Real SQRT_MIN_LN_RATIO=std::sqrt(std::log(1.0/MAX_FREQ_QUAD_RATIO));//assuming a gaussian frequency basis function, this is a useful coonstant for determining the inverse width
    //FIXME: actually compute the necessary value; eval of 0.7 in nearest frequency
    // const Real FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR=1/MAX_FREQ_QUAD_RATIO;//In order to sample the line wings sufficiently in the case of doppler shifts, we increase the distance between the frequencies on the line wings.
    const Size FREQ_QUADRATURE_CENTER_N_ELEMENTS=7;//the number of quadrature elements which should be in the line center
    //TODO?: let this depend on parameters.nquads()?
    const Real FREQ_QUADRATURE_EDGE=1.5;//In order to make sure to sample the boundary condition well for all points, a safety param is included (enlarges the bounds for the frequency quadrature respectively subtracting and adding this to the min and max doppler shifts)
    // const Real FREQ_QUADRATURE_CENTER_SEPARATION_DISTANCE=(2.0*FREQ_QUADRATURE_EDGE)/(FREQ_QUADRATURE_CENTER_N_ELEMENTS-1);//The separation distance of the central frequency quadrature points (in terms of line width)

    Size n_freq_bases=0;//denotes the total number of frequencies we evaluate (=size frequency_quadrature=size frequency_inverse_width_quadrature)
    // vector<Real> min_raydir_doppler_shift;//for every direction, contains the minimum doppler shift in that direction
    // vector<Real> max_raydir_doppler_shift;//for every direction, contains the maximum doppler shift in that direction
    vector<Size> center_line_freq_indices;//for every line, stores the frequency index corresponding to the center of the line.
    vector<Real> frequency_quadrature;//for every line frequency, contains the frequency quadrature
    vector<Real> frequency_inverse_width_quadrature;//for every line frequency, contains the typical inverse width
    // vector<Real> frequency_quadrature_in_line_widths;//for every line frequency, contains the frequency quadrature in terms of line widths compare to the center line frequency
    //not a good idea, might be far too complex and might ruin the exact influence region of boundary conditions...

    //TODO: adaptively define this such that diagonal remains somewhat ?constant...?
    //ideas: should drastically reduce ill-conditionedness of a matrix, so should (implicitly) depend on the number of points/distance of the closest point/...

    //NOTE TO SELF: FREQUENCIES (ESPECIALLY WIDTH) MIGHT CHANGE TO ALSO BE ON THE POINT ID (DUE TO TEMPERATURE DIFFERENCES)
    vector<vector<vector<Real>>> frequencies;//For every direction, for each point, the DOPPLER SHIFTED frequencies of the line transitions
    // THESE ARE NOT ORDERED
    vector<vector<vector<Real>>> frequencies_inverse_widths; //the corresponding inverse widths of the line transition gaussians
    // Dev note: it is defined as 1/(sqrt(2)*sigma) for use in the traditional gaussian

    //not the most efficient way of accessing the default frequencies and widths, but convient.
    vector<Real> non_doppler_shifted_frequencies;//for every frequency, the non doppler shifted line frequency//currently equal to the frequency_quadrature
    vector<Real> non_doppler_shifted_inverse_widths;//for every frequency, the non doppler shifted inverse width//currently equal to frequency_inverse_width_quadrature

    // //TODO: this thing is of the size of the actual matrix; maybe just make an approximation with respect to the actual doppler shift compared to the neighbors...
    // // so simply calculate the max velocity difference, then the max doppler shift is bounded (multiplicatively) by ~1+||Delta v||/c
    // // the advantage will be that we no longer care about the exact direction, nor about the exact other location
    // vector<vector<vector<vector<std::set<Size>>>>> other_frequency_basis_interacting_with; //< contains for every direction, every point, for every interacting rbf with that point, for each frequency basis, the indices of the frequency basis functions interacting with it of the interacting point

    // We make an approsimation by simply calculating the max velocity difference, then the max doppler shift is bounded (multiplicatively) by ~1+||Delta v||/c
    // the advantage will be that we no longer care about the exact direction, nor about the exact other location
    vector<vector<std::set<Size>>> other_frequency_basis_interacting_with; //< contains for every point, for each frequency basis, the indices of the frequency basis functions interacting with it of the interacting point
    // NOTE TO DEV: might change to be also dependent on the exact point (temperature differences)
    // Also, this (and the same thing for the rbf) could be replaced by a vector<set<>>; as we do not particularly care for the order

    const Size N_POINTS_IN_RADIAL_BASIS=8;//the number of points in each radial basis function
    const Real MIN_OPTICAL_DEPTH_DIFF=0.0001;//the minimal optical depth difference between basis functions
    const Real boundary_condition_LS_importance_factor=std::pow(10,4);//

    vector<Size> points_in_grid;//the indices of the points we use in the grid
    Size n_points_in_grid;
    // std::map<Size,Size> indexmap;//maps the compacted indices [0,...,points_in_grid.size()-1] to the actual indice
    vector<Size> n_point_bases;// For all frequencies, stores the number of point bases
    vector<Size> cum_n_point_bases;// For all frequencies, stores the cumulative number of point bases
    // This vector has n_basis_freq+1 elements; starts with 0 and the last element is the total number of point bases.
    vector<vector<Size>> basis_point_idx;// For all frequencies, contains the indices of all positional basis functions
    vector<vector<Size>> local_basis_point_idx;// For all frequencies, contains the indices of all positional basis functions in the collocation grid
    vector<Vector3D> point_locations;// Locations of all points in the collocation grid
    // TODO: when properly implementing these things: check everywhere to add this
    vector<Size> index_conversion; //< for transforming back to the regular point indices used in model (as we might only need a reduced set of all the points)
    // std::map<Size,Size> reverse_index_conversion; //< for transforming from the regular point indices to the point indices used here
    // vector<Real> rbf_radii;//< for every basis point, the typical radius
    vector<vector<Real>> rbf_radii; //<for every frequency, for every basis point, the typical radius
    // vector<std::set<Size>> other_radial_basis_interacting_with; //< contains for each rbf, the indices of the radial basis functions interacting with it (also includes itself)
    // //useful for transforming from basis space to actual solution space
    vector<vector<std::set<Size>>> radial_basis_functions_having_evaluation_at;//< contains for each evaluation frequency, for each point, which basis functions can have a non-zero evaluation at that point
    Eigen::SparseMatrix<Real> collocation_mat; //contains the entries of the sparse collocation matrix
    // Eigen::SparseMatrix<Real> basis_eval_matrix; //contains the entries of the sparse basis matrix for evaluating the intensity
    VectorXr rhs;// The right hand side of the collocation solver


    //for the directions, the basis functions are 1, so no need to store anything

    VectorXr basis_coefficients;
    Eigen :: SparseMatrix<Real> basis_renorm_factor;
    // VectorXr basis_renorm_factor;

    inline void prune_positional_bases(Model& model);
    inline void set_interacting_bases(Model& model);

    //returns whether the frequencies do not lie to far from eachother
    inline bool freq_cutoff_condition(Real abs_freq_diff, Real freq_basis_width);
    // //returns whether the points do not lie to far from eachother
    // inline


public:

    Parameters parameters;

    inline void setup_basis(Model& model);

    //FIXME: change the function definitions to be consistent with the new naming scheme
    //ALSO: check every evaluation of these functions (as some might have their arguments switch places)

    inline Real get_slope(Real x, Size basis_ray_id, Size basis_point_id, Size basis_freq_id);
    inline Real get_slope_derivative(Real x, Size basis_ray_id, Size basis_point_id, Size basis_freq_id);
// The complete radial basis function at an index (dir, nu, pid) is a product of the ones defined below
// basis functions for dir, nu, pid(implicitly x,y,z)
    inline Real basis_direction(Size basis_ray_id);
    inline Real basis_freq(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Real currfreq);
    inline Real basis_point(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry);
    inline Real basis_point_perp(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry);
    inline Real basis_point_symm(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry);
    // basis function derivatives for nu, position
    inline Real basis_freq_der(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Real currfreq);
    inline Real basis_point_der(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry);
    inline Real basis_point_der2(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry);
    inline Real basis_point_symm_der(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry);

    // Integral of the directional basis function over the solid angle
    inline Real basis_direction_int(Geometry& geometry, Size basis_ray_id);
    //Solution of the integral of the basis function together with the line profile function
    inline Real basis_freq_lp_int(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Real lp_freq, Real lp_freq_inv_width);

    //DEPRECATED
    // inline Real distance_manip(Size pointid, Real original_distance);
    // inline Real distance_manip_der(Size pointid, Real original_distance);
// somewhere to store the sparse matrix and the decomposition TODO: figure out what we actually want


// some function to set the giant matrix, starting from the compact basis functions
// psuedocode: first setup the sparse variant:
//  -if using the Phi^T*Phi variant, then for every row/column depending on preferred orientation, do the following:
//    -get the full set of indices which may interact with the basis's which are required at the location (nhat, nu, position);
//      -for simplicity, first define a list of all the indices which interact with with a given index (nhat, nu, position) and augment it with the coefficient.
//      -other option, do the matrix multiplication yourself/in eigen https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
//      - then just read out the resulting sparse matrix into h2lib and try to recompress
//      - also read out sparsity structure (all zero blocks)
//    -kernel matrix approach?
//    -or simpler, just let the software calculate Phi^T*Phi//does not seem to work for H2Lib
//

    // fills the basis triplets with the nonzeros triplets associated with the given point (rayidx, lineidx, pointidx)
    inline void get_nonzero_basis_triplets(std::set<std::vector<Size>>& basis_triplets_to_fill, Size rayid, Size freqid, Size pointid);

    // converts (rayidx, lineidx, pointidx) to the 1D matrix index
    // inline Size get_mat_index(Size rayidx, Size freqidx, Size pointidx);
    inline Size get_mat_row_index(Size rayid, Size freqid, Size pointid);
    inline Size get_mat_col_index(Size basis_ray_id, Size basis_freq_id, Size basis_point_id);
    inline Size get_mat_row_index_2nd_feautrier(Size rayid, Size freqid, Size pointid, bool is_v_eq);
    inline Size get_mat_col_index_2nd_feautrier(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, bool using_v);


// and maybe some helper methods for determining what is actually zero

    //DEPRECATED
    // inline Real compute_balancing_factor_boundary(Size rayidx, Size freqidx, Size pointidx, Real local_velocity_gradient, Model& model);

    inline bool is_boundary_condition_needed_for_direction_on_boundary(Model& model, Size rayid, Size point_on_boundary);

    inline void setup_basis_matrix_Eigen(Model& model);
    inline void setup_rhs_Eigen(Model& model);//Note: assumes one has first setup the basis matrix

    inline void setup_basis_square_matrix_Eigen(Model& model);
    inline void setup_rhs_square_Eigen(Model& model);//Note: assumes one has first setup the basis matrix

    inline void setup_basis_matrix_Eigen_2nd_feautrier(Model& model);
    inline void setup_rhs_Eigen_2nd_feautrier(Model& model);

    inline void rescale_matrix_and_rhs_Eigen(Model& model);

    inline vector<Size> get_two_nearby_points_on_ray(Model& model, Size rayidx, Size pointidx);

    inline Real get_opacity(Model& model, Size rayid, Size freqid, Size pointid);
    inline Real get_opacity_bound(Model& model, Size pointid, Size freqid, Size eval_pointid);
    inline Real get_opacity_grad(Model& model, Size rayid, Size freqid, Size pointid);
    inline Real get_emissivity(Model& model, Size rayid, Size freqid, Size pointid);
    inline Real get_source_grad(Model& model, Size rayid, Size freqid, Size pointid);

    // Return the boundary intensity at the frequency corresponding to the indices
    inline Real boundary_intensity (Model& model, Size rayid, Size freqid, Size pointid);
    // copied from solver
    inline Real planck (Real temp, Real freq);

// Obiviously, we also need some solving methods
    inline void solve_collocation(Model& model);
    //and maybe some variants: Eigen or h2lib

    // After computing the basis coefficients, we can obviously compute the intensity
    // TODO: should be calculated in comoving frame (or you may also doppler shift the line transition frequencies; doppler shifts are the same either way for the line transitions and the other frequencies (if using the same point and ray))
    inline void compute_J(Model& model);
    inline void compute_J_2nd_feautrier(Model& model);



};

#include "collocation.tpp"
