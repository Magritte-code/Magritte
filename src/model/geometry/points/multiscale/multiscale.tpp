#include "tools/types.hpp"
#include <set>
// Coarsen the mesh,
// i.e. add another layer of coarsening.
inline void Multiscale::coarsen()
{//FIXME: also add check whether everything is setup correctly (both functions should be defined)
    mask     .push_back(vector<bool>          (mask     .back()));//should be deep copy
    neighbors.push_back(vector<std::set<Size>>(neighbors.back()));//should be deep copy
    std::set<Size> points_coarsened_around;
    curr_coarsening_lvl=get_max_coars_lvl();

    for (Size p = 0; p < parameters.npoints(); p++)
    {
        if (can_be_coarsened(p, points_coarsened_around))
        {
          coarsen_around_point(p);
          std::cout << "Deleted around point: " << p <<std::endl;
          points_coarsened_around.insert(p);
        }
    }
}


// Returns whether the mesh at a point (p) can be coarsened.
inline bool Multiscale::can_be_coarsened (const Size p, std::set<Size>& points_coarsened_around)
{
    if (!mask.back()[p]||!not_on_boundary(p)) {return false;}//if point no longer in grid, do not coarsen
    //if the point lies on the boundary, do not waste time trying to coarsen around it

    for (const Size n : neighbors.back()[p])
    {
        // Do not coarsen if a neighbor was already coarsend at this level,
        // this avoids creating large holes in the mesh.
        // TODO: replace this....
        // if (!mask.back()[n]) {return false;}
        if (points_coarsened_around.find(n)!=points_coarsened_around.end()) {return false;}

        // Do not coarsen if the required coarsening criterion does not hold.
        // TODO:maybe add some kind of tolerance
        if(!points_are_similar(p,n)) {return false;}
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
      if (not_on_boundary(n))//boundary points will NEVER get removed
        {mask.back()[n] = false;}
    }

    // ..., and replace the neighbors of p (which are removed)
    // in their neighbors by p, using the symmetry of neighbors.
    std::set<Size> neighbors_of_neighbors;
    for (const Size n : neighbors.back()[p])
    {
      if (not_on_boundary(n))
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
    }
    std::cout << "Size neighbors_of_neighbors: " << neighbors_of_neighbors.size() << std::endl;
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
    vector<bool> alltrue(parameters.npoints(), true);
    mask[0]=alltrue;

    //neighbors[0].resize(parameters.npoints());
    vector<std::set<Size>> temp_neighbors;
    temp_neighbors.resize(parameters.npoints());
    Size curr_index=0;
    auto beginvect=std::begin(new_neighbors);//is iterator
    for (Size i=0; i<parameters.npoints(); i++)
    {
      // std::set<Size> curr_neighbors;
      std::set<Size> curr_neighbors(beginvect+curr_index,beginvect+curr_index+n_neighbors[i]);
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

//sets function that checks whether a point does not lie on the boundary
inline void Multiscale::set_not_on_boundary_fun(std::function<bool(Size)> func)
{
  not_on_boundary=func;
}

inline std::set<Size> Multiscale::get_neighbors(const Size p, const Size coars_lvl) const
{//TODO: add check for whether p is still in grid, or just set their neighbors to empty during coarsening
  return neighbors[coars_lvl][p];
}

inline std::set<Size> Multiscale::get_neighbors(const Size p) const
{
  return neighbors[curr_coarsening_lvl][p];
}

inline Size Multiscale::get_nb_neighbors(const Size p, const Size coars_lvl) const
{//TODO: add check for whether p is still in grid, or just set their neighbors to empty during coarsening
  return neighbors[coars_lvl][p].size();
}

inline Size Multiscale::get_nb_neighbors(const Size p) const
{
  return neighbors[curr_coarsening_lvl][p].size();
}

/// Returns the current coarsening level
inline Size Multiscale::get_max_coars_lvl()
{
  if (mask.size()>0)//if already initialized
  {
    return mask.size()-1;
  }
  return 0;//TODO throw error; because not initialized
}

/// Returns the current coarsening level
inline Size Multiscale::get_curr_coars_lvl()
{return curr_coarsening_lvl;}

/// Sets the current coarsening level to lvl (if lvl<=get_max_coars_lvl)
inline void Multiscale::set_curr_coars_lvl(Size lvl)
{
  if (lvl<=get_max_coars_lvl())
  {curr_coarsening_lvl=lvl;}
}

inline vector<bool> Multiscale::get_mask(const Size curr_lvl)
{
  return mask[curr_lvl];
}
