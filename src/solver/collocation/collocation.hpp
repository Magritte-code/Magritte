#pragma once

#include "model/model.hpp"
#include "model/parameters/parameters.hpp"

// Helper class for setting up the collocation solver
class Collocation
{

private:
    //internal parameters for the basis functions
    const Real TRUNCATION_SIGMA=5.0; //< if Delta(freq)>trunc_sigma*freq_width just set the result from the gaussian basis function to zero (as it is almost zero)

    //NOTE TO SELF: FREQUENCIES (ESPECIALLY WIDTH) MIGHT CHANGE TO ALSO BE ON THE POINT ID (DUE TO TEMPERATURE DIFFERENCES)
    vector<Real> frequencies;
    vector<Real> frequencies_widths;

    const Size N_POINTS_IN_RADIAL_BASIS=16;

    vector<Vector3D> point_locations;
    vector<Real> rbf_radii;
    vector<vector<Size>> other_basis_interacting_with; //< contains for each rbf, the indices of the radial basis functions interacting with it (also includes itself)
    //useful for transforming from basis space to actual solution space

    //for the directions, the basis functions are 1, so no need to store anything


public:

    Parameters parameters;

    inline void setup_basis(Model& model);

// The complete radial basis function at an index (dir, nu, pid) is a product of the ones defined below
// basis functions for dir, nu, pid(implicitly x,y,z)
    inline Real basis_direction(Size rayindex);
    inline Real basis_freq(Size freqidx, Real currnu);
    inline Real basis_point(Size pointid, Vector3D& location);

// basis function derivatives for nu, position
    inline Real basis_freq_der(Size freqidx, Real currfreq);
    inline Real basis_point_der(Size centerpoint, Vector3D& location, Size rayindex, Geometry& geometry);

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



// and maybe some helper methods for determining what is actually zero



// Obiviously, we also need some solving methods
    inline void solve_collocation(Model model);
    //and maybe some variants: Eigen or h2lib



};

#include "collocation.tpp"
