#include "tools/types.hpp"
#include <set>
#include <algorithm>
#include <tuple>


///  Sets all neighbors and initializes the data structure
///    @param[in]  n_neighbors: Determines how much neighbors each point has
///    @param[in]  neigbours: A 1D array which contains all neighbors of all points (in order)
///  assumed length = sum of new_n_neighbors
/////////////////////////////////
//TODO: add some form of validation to sizes of new_neighbors and new_n_neighbors
inline void Multiscale::set_all_neighbors(vector<Size>& n_neighbors, vector<Size>& new_neighbors)
{
    neighbors.resize(1);
    mask.resize(1);
    mask[0].resize(parameters.npoints());
    for (Size i=0; i<parameters.npoints(); i++)
    {
        mask[0][i]=true;
    }
    // vector<bool> alltrue(parameters.npoints(), true);
    // mask[0]=alltrue;

    vector<std::set<Size>> temp_neighbors;
    temp_neighbors.resize(parameters.npoints());
    Size curr_index=0;
    auto beginvect=std::begin(new_neighbors);//is iterator
    for (Size i=0; i<parameters.npoints(); i++)
    {
        std::cout<<"mask: "<<(mask[0][i]==true)<<std::endl;
        std::set<Size> curr_neighbors(beginvect+curr_index,beginvect+curr_index+n_neighbors[i]);
        curr_index+=n_neighbors[i];
        temp_neighbors[i]=curr_neighbors;
    }
    neighbors[0]=temp_neighbors;

    set_gpu_neighbors();
}

///  Returns the gpu_compatible reference to the neighbors of a point in the current grid and the amount of neighbors it has
///    @param[in] p: The index of the point to get its neighbors
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline std::tuple<Size*,Size> Multiscale::get_gpu_neighbors(const Size p) const
{
    return std::make_tuple(gpu_neighbors[get_curr_coars_lvl()].dat+gpu_cum_n_neighbors[get_curr_coars_lvl()][p],
                           gpu_n_neighbors[get_curr_coars_lvl()][p]);
}


///  Returns the gpu_compatible reference to the neighbors of a point in the specified grid and the amount of neighbors it has
///    @param[in] p: The index of the point to get its neighbors
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline std::tuple<Size*,Size> Multiscale::get_gpu_neighbors(const Size p, const Size coars_lvl) const
{
    return std::make_tuple(gpu_neighbors[coars_lvl].dat+gpu_cum_n_neighbors[coars_lvl][p],
                           gpu_n_neighbors[coars_lvl][p]);
}


///  Returns the neighbors of a point at the given coarsening level
///    @param[in]  p: Index of the point
///    @param[in]  coars_lvl: The index of the coarsening level
//////////////////////////////////////////////////////////////////
inline std::set<Size> Multiscale::get_neighbors(const Size p, const Size coars_lvl) const
{
    return neighbors[coars_lvl][p];
}

///  Returns the neighbors of a point at the current coarsening level
///    @param[in]  p: Index of the point
///////////////////////////////////////////////////////////////////
inline std::set<Size> Multiscale::get_neighbors(const Size p) const
{
    return neighbors[curr_coarsening_lvl][p];
}

///  Sets the gpu compatible neighbors from the other neighbors
///////////////////////////////////////////////////////////////
inline void Multiscale::set_gpu_neighbors()
{
    gpu_cum_n_neighbors.resize(get_max_coars_lvl()+1);
    gpu_n_neighbors.resize(get_max_coars_lvl()+1);
    gpu_neighbors.resize(get_max_coars_lvl()+1);
    for (Size lvl=0; lvl<=get_max_coars_lvl(); lvl++)
    {
        vector<Size> temp_lin_neighbors=get_all_neighbors_as_vector(lvl);
        gpu_neighbors[lvl].resize(temp_lin_neighbors.size());
        for (Size i=0; i<temp_lin_neighbors.size(); i++)
        {
          gpu_neighbors[lvl][i]=temp_lin_neighbors[i];
        }
        // gpu_neighbors[lvl]=Vector<Size>(get_all_neighbors_as_vector(lvl));
        vector<Size> temp_n_neighbors=get_all_nb_neighbors(lvl);
        gpu_n_neighbors[lvl].resize(parameters.npoints());
        for (Size i=0; i<parameters.npoints(); i++)
        {
          gpu_n_neighbors[lvl][i]=temp_n_neighbors[i];
        }
        // gpu_n_neighbors[lvl]=Vector<Size>(get_all_nb_neighbors(lvl));

        gpu_cum_n_neighbors[lvl].resize(parameters.npoints());
        gpu_cum_n_neighbors[lvl][0] = 0;
        //And finally calculating the cumulative number of neighbors
        for (Size point = 1; point < parameters.npoints(); point++)
        {
            gpu_cum_n_neighbors[lvl][point] = gpu_cum_n_neighbors[lvl][point-1] + gpu_n_neighbors[lvl][point-1];
        }

        gpu_cum_n_neighbors[lvl].copy_vec_to_ptr();
        gpu_n_neighbors[lvl].copy_vec_to_ptr();
        gpu_neighbors[lvl].copy_vec_to_ptr();
    }
}



///  Returns all neighbors at the specified level as a single vector
///    @param[in] coars_lvl: the specified coarsening level
///    @throw invalid_argument: when coars_lvl is higher than get_max_coars_lvl()
///    @return The linearized vector of all neighbors
///////////////////////////////////////////////////////////////////////
inline Size1 Multiscale::get_all_neighbors_as_vector(const Size coars_lvl) const
{
    if (coars_lvl>get_max_coars_lvl())
    {
      throw std::invalid_argument("coarsening level higher than the maximum coarsening level");
    }
    Size nb_points=parameters.npoints();
    Size tot_n_nbs=0;
    for (Size point=0; point<nb_points; point++)
    {
        tot_n_nbs+=get_nb_neighbors(point,coars_lvl);
    }

    Size1 all_neighbors;
    all_neighbors.reserve(tot_n_nbs);
    for (Size point=0; point<nb_points; point++)
    {
        std::set<Size> temp_neighbors=get_neighbors(point,coars_lvl);
        all_neighbors.insert(std::end(all_neighbors), std::begin(temp_neighbors), std::end(temp_neighbors));
    }
    return all_neighbors;
}

///  Returns all neighbors at the finest level as a single vector
///    @return The linearized vector of all neighbors
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



///  Returns the number of neighbors of each point as vector at the specified level
///    @param[in] coars_lvl: the specified coarsening level
///    @throw invalid_argument: when coars_lvl is higher than get_max_coars_lvl()
///    @return The number of neighbors at each point at the finest level
///////////////////////////////////////////////////////////////////////////////////
inline Size1 Multiscale::get_all_nb_neighbors(const Size coars_lvl) const
{
    if (coars_lvl>get_max_coars_lvl())
    {
        throw std::invalid_argument("coarsening level higher than the maximum coarsening level");
    }
    Size nb_points=parameters.npoints();
    Size1 nb_neighbors;
    nb_neighbors.resize(nb_points);
    for (Size point=0; point<nb_points; point++)
    {
        nb_neighbors[point]=get_nb_neighbors(point,coars_lvl);
    }
    return nb_neighbors;
}


///  Returns the number of neighbors of each point as vector at the finest level
///    @param[out] nb_neighbors: The number of neighbors at each point at the finest level
/////////////////////////////////////////////////////////////////////////////////////////
inline Size1 Multiscale::get_all_nb_neighbors() const
{
    Size nb_points=parameters.npoints();
    Size1 nb_neighbors;
    nb_neighbors.resize(nb_points);
    for (Size point=0; point<nb_points; point++)
    {
        nb_neighbors[point]=get_nb_neighbors(point,0);
    }
    return nb_neighbors;
}

///  Returns the number of neighbors of a point at the given coarsening level
///    @param[in]  p: The index of the point
///    @param[in]  coars_lvl: The index of the coarsening level
////////////////////////////////////////////////////////////////////////////
inline Size Multiscale::get_nb_neighbors(const Size p, const Size coars_lvl) const
{
    return neighbors[coars_lvl][p].size();
}

///  Returns the number of neighbors of a point at the current coarsening level
///    @param[in]  p: The index of the point
//////////////////////////////////////////////////////////////////////////////
inline Size Multiscale::get_nb_neighbors(const Size p) const
{
    return neighbors[curr_coarsening_lvl][p].size();
}

///  Returns the current coarsening level
////////////////////////////////////////
inline Size Multiscale::get_max_coars_lvl() const
{
    if (mask.size()>0)//if already initialized
    {
        return mask.size()-1;
    }
    return 0;//If not initialized, the max coarsening level is 0
}

///  Returns the current coarsening level
////////////////////////////////////////
inline Size Multiscale::get_curr_coars_lvl() const
{
    return curr_coarsening_lvl;
}

///  Sets the current coarsening level to the specified lvl (if below the maximum coarsening level)
///    @param[in]  lvl: The coarsening level to set to
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
inline Vector<unsigned char> Multiscale::get_mask(const Size lvl) const
{
    return mask[lvl];
}
/// Returns the total number of points remaining at the given coarsening level
///   @param[in]  lvl: The index of the coarsening level
//////////////////////////////////////////////////////////////////////////////
inline Size Multiscale::get_total_points(const Size lvl)
{//Maybe TODO: if useful, just store this value
    return std::count(mask[lvl].dat, mask[lvl].dat+parameters.npoints(), true);
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
