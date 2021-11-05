#include "io/io.hpp"
#include "model/parameters/parameters.hpp"
#include <set>
#include <map>
#include <tuple>
#include <nanoflann.hpp>

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
    vector<Vector<unsigned char>> mask;
    // Note: Vector<bool> doesnt work due to not having .data()

    // The gpu compatible way of storing the neighbors
    vector<Vector<Size>>     intern_cum_n_neighbors;   ///< cumulative number of neighbors (coars lvl, point)
    vector<Vector<Size>>         intern_n_neighbors;   ///< number of neighbors (coars lvl, point)
    vector<Vector<Size>>           intern_neighbors;   ///< neighbors of each point (coars lvl, point)


    // Sets the neighbors compatible with gpu usage from the other neighbors
    inline void set_intern_neighbors();

    // Returns a reference to the first neighbors in the linearized list in which the neighbors are stores and the amount of neighbors the point has
    inline std::tuple<Size*,Size> get_intern_neighbors(const Size p) const;
    inline std::tuple<Size*,Size> get_intern_neighbors(const Size p, const Size coars_lvl) const;


    //Although the next datastructures are more easy to work with, there is currently no gpu support for vectors of vectors

    // Nearest neighbors for each point in the geometry,
    // first index over levels of coarsening (l),
    // second index over points (p),
    // third index over neighbors (n).
    // Assume the finest mesh (l=0) is given.
    vector<vector<std::set<Size>>> neighbors;

    // Maps a point to the point which deleted it
    // Needed for fallback plan for interpolation
    std::map<Size,Size> point_deleted_map;

    // The current coarsening level
    Size curr_coarsening_lvl=0;
    Size curr_coarsening_lvl_copy=0;//temporary copy for when temporarily returning 0 as level
    bool temporarily_returning_0_as_lvl=false;
    // // For allowing to temporarily return the base level as level to use
    // bool returning_0_as_level=false;


    //returns the current max coarsening level (=size neighbors-1)
    inline Size get_max_coars_lvl() const;
    //returns the current coarsening level
    inline Size get_curr_coars_lvl() const;
    //sets the current coarsening level
    inline void set_curr_coars_lvl(Size lvl);

    //initializes the multiscale class using the given n_neighbors and neighbors vectors
    inline void set_all_neighbors(vector <Size>& new_n_neighbors, vector <Size>& new_neigbours);
    //Maybe TODO: also add support for other ways to initialize this class
    //Also maybe TODO: add support for copying this struct (if we would ever need it)


    //Returns the neighbors of point p, at coarsening level coars_lvl
    inline std::set<Size> get_neighbors(const Size p, const Size coars_lvl) const;
    //Returns the neighbors of point p, at current coarsening level
    inline std::set<Size> get_neighbors(const Size p) const;
    //Returns all neighbors on the specified level as a single vector
    inline Size1 get_all_neighbors_as_vector(const Size coars_lvl) const;
    //Returns all neighbors on the finest level as a single vector
    inline Size1 get_all_neighbors_as_vector() const;
    //Returns the number of neighbors on the finest level as vector
    inline Size1 get_all_n_neighbors(const Size coars_lvl) const;
    //Returns the number of neighbors on the finest level as vector
    inline Size1 get_all_n_neighbors() const;
    //Returns the number of neighbors od point p at coarsening level coars_lvl
    inline Size get_n_neighbors(const Size p, const Size coars_lvl) const;
    //Returns the number of neighbors od point p at the current coarsening level
    inline Size get_n_neighbors(const Size p) const;
    // //Returns the pointer to the start of the mask of points still in the grid at the specified level
    // inline Vector<unsigned char> get_mask(const Size lvl) const;
    //Returns the total number of points remaining at level curr_level
    inline Size get_total_points(const Size lvl);
    //Returns the points in the current grid
    inline Size1 get_current_points_in_grid();
    //Returns the points in the grid at the specified level
    inline Size1 get_points_at_lvl(const Size lvl);


    ///Temporarily sets the level to the the original grid. Important: do not forget to call the next function when done
    inline void temporary_set_level_to_original_grid();
    ///Resets the level to what it should be
    inline void reset_temporary_level();

};

#include "multiscale.tpp"
