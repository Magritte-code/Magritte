#pragma once


#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"

// Structure for neighbors which allows simple insertion and
// deletion of neighbors
struct Neighbors
{
    Parameters parameters;
//TODO set private variables to private
    vector <Size> n_neighbors;  ///< number of neighbors each point has
    Size2          neighbors;   ///< 2d array which contains the neigbors of each point

    inline void delete_single_neighbor(int point, int neighbor);
    inline void delete_all_neighbors(int point);
    inline void add_single_neighbor(int point, int neighbor);

//TODO set to private (probably), not intended for direct use
    inline void set_all_neighbors(vector <Size>& new_n_neighbors,
       vector <Size>& new_neigbours);

    inline vector <Size> get_neighbors (int point) const;
    inline int get_n_neighbors (int point) const;
    inline vector <Size> get_flattened_neigbors_list() const;

//maybe TODO also add constructor that needs 2d vector as input for new_neigbours
Neighbors()=default;

//TODO check whether it is actually a deep copy
Neighbors (const Neighbors& other):
 neighbors(other.neighbors), n_neighbors(other.n_neighbors), parameters(other.parameters){}//deep copy of neighbors

}
;


#include "neighbors.tpp"
