#include "tools/types.hpp"
#include <set>
#include <algorithm>


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

    vector<std::set<Size>> temp_neighbors;
    temp_neighbors.resize(parameters.npoints());
    Size curr_index=0;
    auto beginvect=std::begin(new_neighbors);//is iterator
    for (Size i=0; i<parameters.npoints(); i++)
    {
      std::set<Size> curr_neighbors(beginvect+curr_index,beginvect+curr_index+n_neighbors[i]);
      curr_index+=n_neighbors[i];
      temp_neighbors[i]=curr_neighbors;
    }
    neighbors[0]=temp_neighbors;
}

/// Returns the neighbors of a point at the given coarsening level
///   @param[in]  p: Index of the point
///   @param[in]  coars_lvl: The index of the coarsening level
//////////////////////////////////////////////////////////////////
inline std::set<Size> Multiscale::get_neighbors(const Size p, const Size coars_lvl) const
{
    return neighbors[coars_lvl][p];
}

///Returns the neighbors of a point at the current coarsening level
///   @param[in]  p: Index of the point
///////////////////////////////////////////////////////////////////
inline std::set<Size> Multiscale::get_neighbors(const Size p) const
{
    return neighbors[curr_coarsening_lvl][p];
}

///Returns all neighbors at the finest level as a single vector
///   @param[out] all_neighbors: The linearized vector of all neighbors
///////////////////////////////////////////////////////////////////////
inline Size1 Multiscale::get_all_neighbors_as_vector() const
{
    Size nb_points=parameters.npoints();
    Size tot_n_nbs=0;
    for (Size point=0; point<nb_points; point++)
    {
        tot_n_nbs+=get_nb_neighbors(point,0);
    }
        vector<Size> all_neighbors;
        all_neighbors.reserve(tot_n_nbs);
    for (Size point=0; point<nb_points; point++)
    {
        std::set<Size> temp_neighbors=get_neighbors(point,0);
        all_neighbors.insert(std::end(all_neighbors), std::begin(temp_neighbors), std::end(temp_neighbors));
    }
    return all_neighbors;
}

/// Returns the number of neighbors of each point as vector at the finest level
///   @param[out] nb_neighbors: The number of neighbors at each point at the finest level
/////////////////////////////////////////////////////////////////////////////////////////
inline Size1 Multiscale::get_all_nb_neighbors() const
{
    Size nb_points=parameters.npoints();
    vector<Size> nb_neighbors;
    nb_neighbors.resize(nb_points);
    for (Size point=0; point<nb_points; point++)
    {
        nb_neighbors[point]=get_nb_neighbors(point,0);
    }
    return nb_neighbors;
}

/// Returns the number of neighbors of a point at the given coarsening level
///   @param[in]  p: The index of the point
///   @param[in]  coars_lvl: The index of the coarsening level
////////////////////////////////////////////////////////////////////////////
inline Size Multiscale::get_nb_neighbors(const Size p, const Size coars_lvl) const
{
    return neighbors[coars_lvl][p].size();
}

/// Returns the number of neighbors of a point at the current coarsening level
///   @param[in]  p: The index of the point
//////////////////////////////////////////////////////////////////////////////
inline Size Multiscale::get_nb_neighbors(const Size p) const
{
    return neighbors[curr_coarsening_lvl][p].size();
}

/// Returns the current coarsening level
////////////////////////////////////////
inline Size Multiscale::get_max_coars_lvl()
{
    if (mask.size()>0)//if already initialized
    {
        return mask.size()-1;
    }
    return 0;//If not initialized, the max coarsening level is 0
}

/// Returns the current coarsening level
////////////////////////////////////////
inline Size Multiscale::get_curr_coars_lvl() const
{
    return curr_coarsening_lvl;
}

/// Sets the current coarsening level to the specified lvl (if below the maximum coarsening level)
///   @param[in]  lvl: The coarsening level to set to
//////////////////////////////////////////////////////////////////////////////////////////////////
inline void Multiscale::set_curr_coars_lvl(Size lvl)
{
    if (lvl<=get_max_coars_lvl())
    {
        curr_coarsening_lvl=lvl;
    }
}

/// Returns the mask at the given coarsening level
///   @param[in]  lvl: The index of the coarsening level
////////////////////////////////////////////////////////
inline Bool1 Multiscale::get_mask(const Size lvl) const
{
    return mask[lvl];
}
/// Returns the total number of points remaining at the given coarsening level
///   @param[in]  lvl: The index of the coarsening level
//////////////////////////////////////////////////////////////////////////////
inline Size Multiscale::get_total_points(const Size lvl)
{//Maybe TODO: if useful, just store this value
    return std::count(mask[lvl].begin(), mask[lvl].end(), true);
}

/// Returns the points in the current grid (ordered from low to high)
/////////////////////////////////////////////////////////////////////
inline Size1 Multiscale::get_current_points_in_grid()
{
    Size1 toreturn;
    for (Size idx=0; idx<parameters.npoints(); idx++)
    {
        if (mask[curr_coarsening_lvl][idx])
        {
            toreturn.push_back(idx);
        }
    }
    return toreturn;
}
