#include <iterator>

/// Deletes a single neighbor
///   @param[in]  point: the point from which to delete the neighbor
///   @param[in]  neighbor: the neighbor to delete
////////////////////
inline void Neighbors :: delete_single_neighbor(int point, int neighbor):
{
  if (neighbor<parameters.npoints()&&point<parameters.npoints()):
  {
    self.neighbors[point].erase(std::remove(self.neighbors[point].begin(), self.neighbors[point].end(), neighbor), self.neighbors[point].end());
    self.n_neighbors[point]=vec.size());
    //TODO maybe make assumption that only one neighbor is deleted
  }
  else:
  {
    //Fixme: throw exception
  }

}


/// Deletes a neighbors of a single point
///   @param[in]  point: the point from which to delete all neighbors
////////////////////
inline void Neighbors :: delete_all_neighbors(int point)
{
  if (neighbor<parameters.npoints()&&point<parameters.npoints()):
  {
    self.neighbors[point].clear();
    self.n_neighbors[point]=0;
  }
  else:
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
  if (neighbor<parameters.npoints()&&point<parameters.npoints()):
  {
    self.neighbors[point].push_back(neighbor);
    self.n_neighbors[point]++;
  }
  else:
  {
    //Fixme: throw exception
  }
}

/// Sets all neighbors
///   @param[in]  new_n_neighbors: Determines how much neighbors each point has
///   @param[in]  new_neigbours: A 1D array which contains all neighbors of all points (in order)
/// assumed length = sum of new_n_neighbors
/////////////////////////////////
inline void Neighbors :: set_all_neighbors(Vector <Size> new_n_neighbors, Vector <Size> new_neigbours)
{
  auto length_of_list=accumulate(begin(new_n_neighbors), end(new_n_neighbors), 0, plus<int>())
  if (length_of_list==std::size(new_neigbours)):
  {
    self.n_neighbors=new_n_neighbors;
    curr_index=0;
    for (i=0, i<parameters.npoints(),i++):
    {
      self.neighbors[i]=sub(&new_neigbours[curr_index],&new_neigbours[curr_index+self.n_neighbors[i]]);
      curr_index+=self.n_neighbors[i]
    }
  }
  else:
  {
      //Fixme: throw exception, because size of new_neigbours is not correct
  }
}


/// Returns the neigbours of a point
///   @param[in]  point: the point for which to return its neighbors
///   @return: the vector of neighbors of the point
//////////////////
inline Vector<Size> Neighbors :: get_neighbors (int point)
{
  return self.neighbors[point];
}

/// Deep copy of the Neighbors construct
/////////////////////
Neighbors::Neighbors(Neighbors& other)
{
  auto n_neighbors(other.n_neighbors);
  auto temp_neighbors(other.neighbors)
  self.n_neighbors=temp_n_neighbors;
  self.temp_neighbors=temp_neighbors;
}


}
