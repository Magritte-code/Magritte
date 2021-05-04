#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include <set>
#include <map>
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

    // Maps a point to the point which deleted it
    // Needed for fallback plan for interpolation
    std::map<Size,Size> point_deleted_map;

    ///DEPRECATED
    // // Function to check whether two points are similar enough to be able to coarsen them
    // // Should be a lambda function supplied by something else.
    // // Maybe also add tolerance to it TODO
    // // Should return true when we may coarsen the grid
    // // Should also take into account wehther there are boundary points involved
    // std::function<bool(Size, Size)> points_are_similar;
    //
    // ///FIXME: also take care to NOT delete boundary points; maybe just take the same reference as geometry
    //
    // std::function<bool(Size)> not_on_boundary;

    //current coarsening level
    Size curr_coarsening_lvl=0;

    ///DEPRECATED
    // // Just some safety precations that should be true before calling coarsen
    // bool boundary_set=false;
    // bool comparison_set=false;

    ///DEPRECATED, moved to Model
    // // Coarsen the mesh,
    // // i.e. add another layer of coarsening.
    // inline void coarsen ();
    //
    // // Returns whether the mesh at a point (p) can be coarsened.
    // inline bool can_be_coarsened (const Size p, std::set<Size>& points_coarsened_around);
    //
    // // Coarsens the neighbors of p and updates the neighbors of p and neighbors of the neighbors of neighbors
    // inline void coarsen_around_point (const Size p);

    //returns the current max coarsening level (=size neighbors-1)
    inline Size get_max_coars_lvl();
    //returns the current coarsening level
    inline Size get_curr_coars_lvl() const;
    //sets the current coarsening level
    inline void set_curr_coars_lvl(Size lvl);

    //initializes the multiscale class using the given n_neighbors and neighbors vectors
    inline void set_all_neighbors(vector <Size>& new_n_neighbors, vector <Size>& new_neigbours);
    //Maybe TODO: also add support for other ways to initialize this class
    //Also maybe TODO: add support for copying this struct (if we would ever need it)

    // //sets the comparison function to see whether points are similar enough
    // inline void set_comparison_fun(std::function<bool(Size, Size)> func);
    // //sets function that checks whether a point does not lie on the boundary
    // inline void set_not_on_boundary_fun(std::function<bool(Size)> func);

    //Returns the neighbors of point p, at coarsening level coars_lvl
    inline std::set<Size> get_neighbors(const Size p, const Size coars_lvl) const;
    //Returns the neighbors of point p, at current coarsening level
    inline std::set<Size> get_neighbors(const Size p) const;
    //Returns all neighbors on the finest level as a single vector
    inline vector<Size> get_all_neighbors_as_vector() const;
    //Returns the number of neighbors on the finest level as vector
    inline vector<Size> get_all_nb_neighbors() const;
    //Returns the number of neighbors od point p at coarsening level coars_lvl
    inline Size get_nb_neighbors(const Size p, const Size coars_lvl) const;
    //Returns the number of neighbors od point p at the current coarsening level
    inline Size get_nb_neighbors(const Size p) const;
    //Returns the mask of points still in the grid at level curr_level
    inline vector<bool> get_mask(const Size curr_lvl) const;
    //Returns the total number of points remaining at level curr_level
    inline Size get_total_points(const Size curr_lvl);
    //Returns the points in the current grid
    inline vector<Size> get_current_points_in_grid();

};

#include "multiscale.tpp"
