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


/// Adds a neighbor to a point and vice versa
///   @param[in]  point: the point to which to add the neighbor
///   @param[in]  neighbor: the neighbor to add
////////////////////
inline void Neighbors :: add_neighbor(int point, int neighbor)
{
  if (neighbor<parameters.npoints()&&point<parameters.npoints())
  {
    neighbors[point].push_back(neighbor);
    n_neighbors[point]++;
    neighbors[neighbor].push_back(point);
    n_neighbors[neighbor]++;
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
//currently has a workaround for points having multiple times the same neighbor
inline void Neighbors :: set_all_neighbors(vector <Size> &new_n_neighbors, vector <Size> &new_neighbors)
{
// Ignoring checking whether the lengths match for now: The paracabs vector does not support begin nor end
//  auto length_of_list=accumulate(std::begin(new_n_neighbors), std::end(new_n_neighbors), 0, std::plus<Size>())
//  if (length_of_list==std::size(new_neigbours))
//  {
    n_neighbors=new_n_neighbors;
    neighbors.resize(parameters.npoints());
    int curr_index=0;
    for (Size i=0; i<parameters.npoints(); i++)
    {
      //neighbors[i].resize(n_neighbors[i]);
      for (Size j=0; j<n_neighbors[i]; j++)
      {//workaround
        if(std::find(neighbors[i].begin(), neighbors[i].end(), new_neighbors[curr_index]) == neighbors[i].end()) {
          neighbors[i].push_back(new_neighbors[curr_index]);//neighbors[i][j]=...
        }
        curr_index++;
      }
      n_neighbors[i]=neighbors[i].size();
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
inline vector <Size> Neighbors :: get_neighbors (int point) const
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
inline vector <Size> Neighbors :: get_flattened_neigbors_list() const
{
  Size total_size = 0;
  for (const auto& part : neighbors)
      {total_size += part.size();}//precalc the size of the resulting array
  vector<Size> result;
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


/// Default constructor
//////////////////////
// Neighbors::Neighbors()
// {
//   Parameters parameters=Parameters();
//   vector <Size> n_neighbors= vector<Size>();
//   Size2           neighbors= Size2();
//
// }
/// Deep copy of the Neighbors construct
/////////////////////
/*Neighbors::Neighbors(const Neighbors& other)
{
  Size2 temp_neighbors(other.neighbors);
  Size1 temp_n_neighbors(other.n_neighbors);
  n_neighbors=temp_n_neighbors;
  neighbors=temp_neighbors;
  Parameters temp_params(other.parameters);
  parameters=temp_params;
  // n_neighbors=std::copy(other.n_neighbors.begin(), other.n_neighbors.end(), back_inserter(n_neighbors));
  // neighbors=std::copy(other.neighbors.begin(), other.neighbors.end(), back_inserter(neighbors));;
  // parameters= new Parameters();
}*/
