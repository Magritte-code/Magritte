#include <cmath>

inline double Model :: calc_diff_abundance_with_neighbours(int point, int next_coars_lvl)
{
  double temp_max=0;//current implementation is probably slow
  for (auto& neigbor:neighbors_lists[next_coars_lvl].get_neighbors(point))
    {
      double abundance_self=chemistry.species.abundance[point][0];
      double abundance_other=chemistry.species.abundance[neigbor][0];
      double temp_diff=std::abs((abundance_self-abundance_other)/(abundance_self+abundance_other));
      temp_max=std::max(temp_max, temp_diff);
    }
    return temp_max;
}




inline void Model :: coursen_grid(const float perc_points_deleted)
{
  if (curr_coarsening_lvl==max_reached_coarsening_lvl)//if we truly need to refine the grid
  {
    if (max_reached_coarsening_lvl==0)//first time refining grid
      {
      current_nb_points=parameters.npoints();
      neighbors_lists.resize(2);
      neighbors_lists[0]=Neighbors(geometry.points.curr_neighbors);//TODO properly implement deep copy
      neighbors_lists[1]=Neighbors(geometry.points.curr_neighbors);//TODO properly implement deep copy
      density_diff_at_point.resize(parameters.npoints());
      for (int i=0; i<parameters.npoints(); i++)//calc abs(f-f_other)/(f+f_other) for all points
        {
          density_diff_at_point[i]=calc_diff_abundance_with_neighbours(i,1);
        }
      }
      else
      {
        neighbors_lists.resize(max_reached_coarsening_lvl+1);
        neighbors_lists[max_reached_coarsening_lvl+1]=Neighbors(neighbors_lists[max_reached_coarsening_lvl]);//should be deep copy
      }
  max_reached_coarsening_lvl++;
  curr_coarsening_lvl++;
  //TODO
  Size nb_points_to_remove=Size(perc_points_deleted*current_nb_points);
  for (Size i=0; i<nb_points_to_remove; i++)
  {
      //TODO remove points
  }
  //will be harder to parallelize
  }
  else
  {//TODO just change current Neighbors
    curr_coarsening_lvl++;
    geometry.points.curr_neighbors=neighbors_lists[1];//should be shallow copy
  }
}
