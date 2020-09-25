#include <iterator>
#include "tools/types.hpp"

/// Deletes a single neighbor
///   @param[in]  point: the point from which to delete the neighbor
///   @param[in]  neighbor: the neighbor to delete
////////////////////
inline void Neighbors :: delete_single_neighbor(int point, int neighbor)
{
  if ((neighbor<parameters.npoints())&&(point<parameters.npoints()))
  {
    neighbors[point].erase(std::remove(neighbors[point].begin(), neighbors[point].end(), neighbor), neighbors[point].end());
    n_neighbors[point]=neighbors[point].size();
    //TODO maybe make assumption that only one neighbor is deleted
  }
  else
  {
    //Fixme: throw exception
  }
}


/// Deletes a neighbors of a single point
///   @param[in]  point: the point from which to delete all neighbors
////////////////////
inline void Neighbors :: delete_all_neighbors(int point)
{
  if (point<parameters.npoints())
  {
    neighbors[point].clear();
    n_neighbors[point]=0;
  }
  else
  {
    //Fixme: throw exception
  }
}


/// Adds a single neighbor to a point
///   @param[in]  point: the point to which to add the neighbor
///   @param[in]  neighbor: the neighbor to add
////////////////////
inline void Neighbors :: add_single_neighbor(int point, int neighbor)
{
  if (neighbor<parameters.npoints()&&point<parameters.npoints())
  {
    neighbors[point].push_back(neighbor);
    n_neighbors[point]++;
  }
  else
  {
    //Fixme: throw exception
  }
}

/// Sets all neighbors
///   @param[in]  new_n_neighbors: Determines how much neighbors each point has
///   @param[in]  new_neigbours: A 1D array which contains all neighbors of all points (in order)
/// assumed length = sum of new_n_neighbors
/////////////////////////////////
inline void Neighbors :: set_all_neighbors(Vector <Size> new_n_neighbors, Vector <Size> new_neighbors)
{
// Ignoring checking whether the lengths match for now: The paracabs vector does not support begin nor end
//  auto length_of_list=accumulate(std::begin(new_n_neighbors), std::end(new_n_neighbors), 0, std::plus<Size>())
//  if (length_of_list==std::size(new_neigbours))
//  {
    cout << "flag5" << endl;
    n_neighbors.resize(parameters.npoints());
    cout << "flag5" << endl;
    cout << new_neighbors[0] << endl;
    cout << "flag5" << endl;
    for (int i=0; i<parameters.npoints(); i++)
    {
      n_neighbors[i]=new_n_neighbors[i];
    }// I cant figure out why the allocation fails....
    //n_neighbors=new_n_neighbors.dat;
    cout << "flag6" << endl;
    neighbors.resize(parameters.npoints());
    cout << "flag7" << endl;
    int curr_index=0;
    for (Size i=0; i<parameters.npoints(); i++)
    {
      cout << "flag8" << endl;
      cout << n_neighbors[0] << endl;
      cout << n_neighbors[i] << endl;
      neighbors[i].resize(n_neighbors[i]);
          cout << "flag8" << endl;
      for (Size j=0; j<n_neighbors[i]; j++)
      {
        neighbors[i][j]=new_neighbors[curr_index];
        curr_index++;
      }
      //neighbors[i]=std::sub(&new_neigbours[curr_index],&new_neigbours[curr_index+n_neighbors[i]]);
      //curr_index+=n_neighbors[i];
    }
//  }
//  else
//  {
      //Fixme: throw exception, because size of new_neigbours is not correct
//  }
}


/// Returns the neigbours of a point
///   @param[in]  point: the point for which to return its neighbors
///   @return: the vector of neighbors of the point
//////////////////
inline Vector <Size> Neighbors :: get_neighbors (int point) const
{
  return neighbors[point];
}

/// Returns the number of neigbours of a point
///   @param[in]  point: the point for which to return its #neighbors
///   @return: number of neighbors of the point
//////////////////
inline int Neighbors :: get_n_neighbors (int point) const
{
  return n_neighbors[point];
}

/// Returns the flattened neighbors list
///   @return: the vector of neighbors of all points
//////////////////
inline Vector <Size> Neighbors :: get_flattened_neigbors_list() const
{
  Size total_size = 0;
  for (const auto& part : neighbors)
      {total_size += part.size();}//precalc the size of the resulting array
  Vector<Size> result;
  result.resize(total_size);
  int curr_index=0;
  for (const auto& part : neighbors)
      {for (int j=0; j<part.size(); j++)
        {result[curr_index]=part[j];
        curr_index++;}
      }//Vector does not support .begin(), nor .end()
      //{result.insert(result.end(), sub.begin(), sub.end());}
  return result;
}


/// Deep copy of the Neighbors construct
/////////////////////
/*Neighbors::Neighbors(Neighbors& other)
{
  auto n_neighbors(other.n_neighbors);
  auto temp_neighbors(other.neighbors)
  n_neighbors=temp_n_neighbors;
  neighbors=temp_neighbors;
}*/
