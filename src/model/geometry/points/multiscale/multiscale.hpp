#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include <set>
///First implementation of the multiscale class
/// This class is meant to easily coarsen the grid stored in it.
struct Multiscale
{
    Parameters parameters;
    // Mask to indicate which points are still active at
    // a certain level of coarsening. True if still
    // in mesh, i.e. for the finest mesh (l=0) this is all true.
    // First index over levels of coarsening (l),
    // second index over points (p).
    // Assume the finest mesh (l=0) is given.
    vector<vector<bool>> mask;

    // Nearest neighbors for each point in the geometry,
    // first index over levels of coarsening (l),
    // second index over points (p),
    // third index over neighbors (n).
    // Assume the finest mesh (l=0) is given.
    vector<vector<std::set<Size>>> neighbors;

    // Function to check whether two points are similar enough to be able to coarsen them
    // Should be a lambda function supplied by something else.
    // Maybe also add tolerance to it TODO
    // Should return true when we may coarsen the grid
    // Should also take into account wehther there are boundary points involved
    std::function<bool(Size, Size)> points_are_similar;


    // Coarsen the mesh,
    // i.e. add another layer of coarsening.
    inline void coarsen ();

    // Returns whether the mesh at a point (p) can be coarsened.
    inline bool can_be_coarsened (const Size p);

    // Coarsens the neighbors of p and updates the neighbors of p and neighbors of the neighbors of neighbors
    inline void coarsen_around_point (const Size p);

    //returns the current coarsening level (=size neighbors-1)
    inline Size get_curr_coars_lvl();

    //initializes the multiscale class using the given n_neighbors and neighbors vectors
    inline void set_all_neighbors(vector <Size>& new_n_neighbors, vector <Size>& new_neigbours);
    //Maybe TODO: also add support for other ways to initialize this class
    //Also maybe TODO: add support for copying this struct (if we would ever need it)

    //sets the comparison function to see whether points are similar enough
    inline void set_comparison_fun(std::function<bool(Size, Size)>);

    //Returns the neighbors of point p, at coarsening level coars_lvl
    inline std::set<Size> get_neighbors(Size p, Size coars_lvl);
};

#include "multiscale.tpp"
