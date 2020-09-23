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
    Size_2          neighbors;  ///


    void delete_single_neighbor(int point, int neighbor);
    void delete_all_neighbors(int point);
    void add_single_neighbor(int point, int neighbor);

//TODO set to private (probably), not intended for direct use
    void set_all_neighbors(Vector <Size> new_n_neighbors, Vector <Size> new_neigbours);


    Vector <Size> get_neighbors (int point);

//TODO: add constructor
    Neighbors (Neighbors& other);//copy neighbors

}



#include "neighbors.tpp"