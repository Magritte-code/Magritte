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
    //internal parameters for the basis functions
    const Real TRUNCATION_SIGMA=5.0; //< if Delta(freq)>trunc_sigma*freq_width just set the result from the gaussian basis function to zero (as it is almost zero)
    const Real SLOPE_FACTOR=1000.0;//0.0;

    //NOTE TO SELF: FREQUENCIES (ESPECIALLY WIDTH) MIGHT CHANGE TO ALSO BE ON THE POINT ID (DUE TO TEMPERATURE DIFFERENCES)
    vector<vector<vector<Real>>> frequencies;//For every direction, for each point, the DOPPLER SHIFTED frequencies of the line transitions
    // THESE ARE NOT ORDERED
    vector<vector<vector<Real>>> frequencies_inverse_widths; //the corresponding inverse widths of the line transition gaussians
    // Dev note: it is defined as 1/(sqrt(2)*sigma) for use in the traditional gaussian
    vector<Real> non_doppler_shifted_frequencies;

    // //TODO: this thing is of the size of the actual matrix; maybe just make an approximation with respect to the actual doppler shift compared to the neighbors...
    // // so simply calculate the max velocity difference, then the max doppler shift is bounded (multiplicatively) by ~1+||Delta v||/c
    // // the advantage will be that we no longer care about the exact direction, nor about the exact other location
    // vector<vector<vector<vector<std::set<Size>>>>> other_frequency_basis_interacting_with; //< contains for every direction, every point, for every interacting rbf with that point, for each frequency basis, the indices of the frequency basis functions interacting with it of the interacting point

    // We make an approsimation by simply calculating the max velocity difference, then the max doppler shift is bounded (multiplicatively) by ~1+||Delta v||/c
    // the advantage will be that we no longer care about the exact direction, nor about the exact other location
    vector<vector<std::set<Size>>> other_frequency_basis_interacting_with; //< contains for every point, for each frequency basis, the indices of the frequency basis functions interacting with it of the interacting point
    // NOTE TO DEV: might change to be also dependent on the exact point (temperature differences)
    // Also, this (and the same thing for the rbf) could be replaced by a vector<set<>>; as we do not particularly care for the order

    const Size N_POINTS_IN_RADIAL_BASIS=16;

    vector<Vector3D> point_locations;
    vector<Size> index_conversion; //< for transforming back to the regular point indices used in model (as we might only need a reduced set of all the points)
    std::map<Size,Size> reverse_index_conversion; //< for transforming from the regular point indices to the point indices used here
    vector<Real> rbf_radii;
    vector<std::set<Size>> other_radial_basis_interacting_with; //< contains for each rbf, the indices of the radial basis functions interacting with it (also includes itself)
    //useful for transforming from basis space to actual solution space

    Eigen::SparseMatrix<Real> collocation_mat; //contains the entries of the sparse collocation matrix
    // Eigen::SparseMatrix<Real> basis_eval_matrix; //contains the entries of the sparse basis matrix for evaluating the intensity
    VectorXr rhs;// The right hand side of the collocation solver


    //for the directions, the basis functions are 1, so no need to store anything

    VectorXr basis_coefficients;


    inline void set_interacting_bases(Model& model);

    //returns whether the frequencies do not lie to far from eachother
    inline bool freq_cutoff_condition(Real abs_freq_diff, Real freq_basis_width);
    // //returns whether the points do not lie to far from eachother
    // inline


public:

    Parameters parameters;

    inline void setup_basis(Model& model);

// The complete radial basis function at an index (dir, nu, pid) is a product of the ones defined below
// basis functions for dir, nu, pid(implicitly x,y,z)
    inline Real basis_direction(Size rayindex);
    inline Real basis_freq(Size rayidx, Size lineidx, Size pointidx, Real currfreq);
    inline Real basis_point(Size pointid, Vector3D& location, Size rayindex, Geometry& geometry);

    // basis function derivatives for nu, position
    inline Real basis_freq_der(Size rayidx, Size lineidx, Size pointidx, Real currfreq);
    inline Real basis_point_der(Size centerpoint, Vector3D& location, Size rayindex, Geometry& geometry);

    // Integral of the directional basis function over the solid angle
    inline Real basis_direction_int(Geometry& geometry, Size rayindex);
    //Solution of the integral of the basis function together with the line profile function
    inline Real basis_freq_lp_int(Size rayidx, Size lineidx, Size pointidx, Real lp_freq, Real lp_freq_inv_width);

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
    inline void get_nonzero_basis_triplets(std::set<std::vector<Size>>& basis_triplets_to_fill, Size rayidx, Size lineidx, Size pointidx);

    // converts (rayidx, lineidx, pointidx) to the 1D matrix index
    inline Size get_mat_index(Size rayidx, Size lineidx, Size pointidx);
// and maybe some helper methods for determining what is actually zero

    inline void setup_basis_matrix_Eigen(Model& model);
    inline void setup_rhs_Eigen(Model& model);//Note: assumes one has first setup the basis matrix

    inline Real get_opacity(Model& model, Size rayidx, Size lineidx, Size pointidx);
    inline Real get_emissivity(Model& model, Size rayidx, Size lineidx, Size pointidx);

    // Return the boundary intensity at the frequency corresponding to the indices
    inline Real boundary_intensity (Model& model, Size rayidx, Size freqidx, Size pointidx);
    // copied from solver
    inline Real planck (Real temp, Real freq);

// Obiviously, we also need some solving methods
    inline void solve_collocation(Model& model);
    //and maybe some variants: Eigen or h2lib

    // After computing the basis coefficients, we can obviously compute the intensity
    // TODO: should be calculated in comoving frame (or you may also doppler shift the line transition frequencies; doppler shifts are the same either way for the line transitions and the other frequencies (if using the same point and ray))
    inline void compute_J(Model& model);



};

#include "collocation.tpp"
