#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"
#include <set>
// Coarsen the mesh,
// i.e. add another layer of coarsening.
inline void Multiscale::coarsen()
{
    mask     .push_back(vector<bool>          (mask     .back()));//should be deep copy
    neighbors.push_back(vector<std::set<Size>>(neighbors.back()));//should be deep copy
    for (Size p = 0; p < parameters.npoints(); p++)
    {
        if (can_be_coarsened(p)) {coarsen_around_point(p);}
    }
}


// Returns whether the mesh at a point (p) can be coarsened.
inline bool Multiscale::can_be_coarsened (const Size p)
{
    if (!mask.back()[p]) return false;//if point no longer in grid, do not coarsen

    for (const Size n : neighbors.back()[p])
    {
        // Do not coarsen if a neighbor was already coarsend at this level,
        // this avoids creating large holes in the mesh.
        if (!mask.back()[n]) return false;

        // Do not coarsen if the required coarsening criterion does not hold.
        // TODO:maybe add some kind of tolerance
        if(!points_are_similar(p,n)) return false;
    }

    // Otherwise coarsen.
    return true;
}


inline void Multiscale::coarsen_around_point (const Size p)
{
    // New neighbors for p.
    std::set<Size> new_neighbors;

    // Identify all neighbors with the current point (p),
    // i.e. remove its neighbors from the mesh by masking,
    // TODO: figure out whether it is ok/necessary to empty their neighbors
    for (const Size n : neighbors.back()[p])
    {
        mask.back()[n] = false;
    }

    // ..., and replace the neighbors of p (which are removed)
    // in their neighbors by p, using the symmetry of neighbors.
    std::set<Size> neighbors_of_neighbors;
    for (const Size n : neighbors.back()[p])
    {
      for (const Size n_n : neighbors.back()[n])
      {
        if (mask.back()[n_n]&&n_n!=p)//if neighbor of neighbor is still in the grid, i.e. not just a neighbor; also, do never try to add a point as its own neighbor!!!!
        {
          neighbors_of_neighbors.insert(n_n);
          neighbors.back()[n_n].erase(n);//remove n from neighbors of n_n
        }
      }
    }
    for (const Size n_n:neighbors_of_neighbors)
    {
      // Replace the removed points by p
      neighbors.back()[n_n].insert(p);

      // And add the others as new neighbors.
      new_neighbors.insert(n_n);
    }

    // Set the new nearest neighbors.
    neighbors.back()[p] = new_neighbors;
}

/// Sets all neighbors and initializes the data structure
///   @param[in]  n_neighbors: Determines how much neighbors each point has
///   @param[in]  neigbours: A 1D array which contains all neighbors of all points (in order)
/// assumed length = sum of new_n_neighbors
/////////////////////////////////
//TODO: add some form of validation to sizes of new_neighbors and new_n_neighbors
inline void Multiscale::set_all_neighbors(vector<Size>& n_neighbors, vector<Size>& new_neighbors)
{
    neighbors.resize(1);
    mask.resize(1);
    std::vector<bool> alltrue(parameters.npoints(), true);
    mask[0]=alltrue;

    //neighbors[0].resize(parameters.npoints());
    vector<std::set<Size>> temp_neighbors;
    temp_neighbors.resize(parameters.npoints())
    Size curr_index=0;
    auto beginvect=std::begin(new_neighbors);
    for (Size i=0; i<parameters.npoints(); i++)
    {
      // std::set<Size> curr_neighbors;
      std::set<Size> curr_neighbors(beginvect+curr_index,beginvect+curr_index+n_neighbors[i])
      curr_index+=n_neighbors[i];
      // for (Size j=0; j<n_neighbors[i]; j++)
      // {
      //   curr_neighbors.insert(new_neighbors[curr_index])
      //   curr_index++;
      // }
      temp_neighbors[i]=curr_neighbors;
    }
    neighbors[0]=temp_neighbors;
}

inline void Multiscale::set_comparison_fun(std::function<bool(Size,Size)> func)
{
    points_are_similar=func;
}

inline std::set<Size> Multiscale::get_neighbors(Size p, Size coars_lvl)
{//TODO: add check for whether p is still in grid, or just set their neighbors to empty during coarsening
    return neighbors[coars_lvl][p];
}
