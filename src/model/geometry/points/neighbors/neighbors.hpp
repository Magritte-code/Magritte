#pragma once


#include "model/parameters/parameters.hpp"
#include "tools/types.hpp"

// Structure for neighbors which allows simple insertion and
// deletion of neighbors
struct Neighbors
{
    Parameters parameters;
//TODO set private variables to private
    Vector <Size> n_neighbors;  ///< number of neighbors each point has
    Size2          neighbors;  ///

    inline void delete_single_neighbor(int point, int neighbor);
    inline void delete_all_neighbors(int point);
    inline void add_single_neighbor(int point, int neighbor);

//TODO set to private (probably), not intended for direct use
    private inline void set_all_neighbors(Vector <Size> new_n_neighbors,
       Vector <Size> new_neigbours);

    inline Vector <Size> get_neighbors (int point);
    inline int get_n_neighbors (int point);
    inline Vector <Size> get_flattened_neigbors_list();

//constuctor
    Neighbors(Vector <Size> new_n_neighbors, Vector <Size> new_neigbours);
//maybe TODO also add constructor that needs 2d vector as input for new_neigbours

//TODO: add constructor
    Neighbors (Neighbors& other);//copy neighbors

}



#include "neighbors.tpp"
