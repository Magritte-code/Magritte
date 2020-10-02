#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cfloat>
#include <set>
#include "tools/types.hpp"

//NOTE: i have mistakenly called tetrahedra triangles throughout this entire piece of code

// Calculates the relative difference of the density with respect to the neighbours
//    @Returns: abs((d-d_other)/(d+d_other)) such that the smallest values will be assigned
// to the points we want to delete
inline double Model :: calc_diff_abundance_with_neighbours(Size point, Size next_coars_lvl)
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


///Helper function for deleting stuff from the density_diff_maps
///   @Param[in/out] std_map: the map with points as keys
///   @Param[in/out] reverse_map: the map with points as values
///   @Param[in] point_to_del: the point to delete from both maps
inline void delete_int_from_both_maps(std::multimap<Size,double> &std_map, std::multimap<double,Size> &reverse_map, Size point_to_del)
{
  double curr_value=(*std_map.find(point_to_del)).second;
  std_map.erase(point_to_del);//key is unique, so no problem
  //there are multiple points that currently have the same diff value, so we must iterate
  typedef std::multimap<double,Size>::iterator iterator;
  std::pair<iterator, iterator> iterpair = reverse_map.equal_range(curr_value);
  iterator it = iterpair.first;
  for (; it != iterpair.second; ++it) {
    if (it->second == point_to_del) {
      reverse_map.erase(it);
      break;
    }
  }
}

///Helper function for deleting stuff from the ears_maps
///   @Param[in/out] std_map: the map with vectors as keys
///   @Param[in/out] reverse_map: the map with vectors as values
///   @Param[in] vect_to_del: the vector to delete from both maps
///   @Returns: Returns an iterator for the standard map (needed for being able to delete while iterating)
inline std::multimap<vector<Size>,double>::iterator delete_vector_from_both_maps(
  std::multimap<vector<Size>,double> &std_map, std::multimap<double,vector<Size>> &reverse_map, vector<Size> vect_to_del)
{
  double curr_value=(*std_map.find(vect_to_del)).second;
  std::multimap<vector<Size>,double>::iterator it2=std_map.erase(std_map.find(vect_to_del));//key is unique, so no problem
  //there are multiple points that currently have the same diff value, so we must iterate
  typedef std::multimap<double,vector<Size>>::iterator iterator;
  std::pair<iterator, iterator> iterpair = reverse_map.equal_range(curr_value);
  iterator it = iterpair.first;
  for (; it != iterpair.second; ++it) {
    if (it->second == vect_to_del) {
      reverse_map.erase(it);
      break;
    }
  }
  return it2;
}

///Helper function: checks whether a vector contains an element
///   @Param[in] vect: the vector to check for
///   @Param[in] element: the element to check for
///   @Returns bool: true if element in vect, false otherwise
inline bool vector_contains_element(const vector<Size> vect, Size element)
{
  return (std::find(vect.begin(), vect.end(), element) != vect.end());
}



/// Given a tetrahedron, this calculates the power of the point with respect to the circumsphere
/// @Parameter [in] triangle: the tetrahedron from which we use the circumsphere
/// @Para
/// returns a positive value if the point is inside the circumsphere and a negative value if the point lies outside the circumsphere (zero if on the circumsphere)
inline double Model :: calc_power(const vector<Size> &triangle, Size point){
  Vector3D pos1=geometry.points.position[triangle[0]];
  Vector3D pos2=geometry.points.position[triangle[1]];
  Vector3D pos3=geometry.points.position[triangle[2]];
  Vector3D pos4=geometry.points.position[triangle[3]];
  Vector3D posp=geometry.points.position[point];//position of point

  //dividing insphere test with orientation test
  Eigen::Matrix<double,5,5> insphere;
  insphere << pos1.x(),pos2.x(),pos3.x(),pos4.x(),posp.x(),
              pos1.y(),pos2.y(),pos3.y(),pos4.y(),posp.y(),
              pos1.z(),pos2.z(),pos3.z(),pos4.z(),posp.z(),
              pos1.squaredNorm(),pos2.squaredNorm(),pos3.squaredNorm(),pos4.squaredNorm(),posp.squaredNorm(),
              1,1,1,1,1;
  Eigen::Matrix<double,4,4> orient;
  orient << pos1.x(),pos2.x(),pos3.x(),pos4.x(),
            pos1.y(),pos2.y(),pos3.y(),pos4.y(),
            pos1.z(),pos2.z(),pos3.z(),pos4.z(),
            1,1,1,1;

  return insphere.determinant()/orient.determinant();
}


//generates the new ears and inserts them into the ears maps
//TODO add description
inline void Model::generate_new_ears(const vector<Size> &neighbors_of_point,const vector<Size> &plane,std::map<Size, std::set<Size>> &neighbor_map,
  std::multimap<vector<Size>,double> &ears_map, std::multimap<double,vector<Size>> &rev_ears_map, Size &curr_point)
{
  Size count_neighbours;
  for (Size temp_point: neighbors_of_point)
  { //we want to create new triangles, therefore a point in the plane is not useful
    if (!vector_contains_element(plane,temp_point))
    {
      //we are currently adding far too many new ears??
      count_neighbours=0;
      for (Size point_of_plane: plane)
      {
        if (neighbor_map[temp_point].count(point_of_plane)!=0)//if temp_point is neighbor of point_of_plane
        {count_neighbours++;}
      }
      //if the candidate point is good for creating an ear with plane1
      if (count_neighbours==2)
      {
        Size point_not_neighbor_of_plane;//the point of the plane which is not a neighbor
        vector<Size> points_neighbor_of_plane;
        for (Size point_of_plane: plane)
        {
          if (neighbor_map[temp_point].count(point_of_plane)==0)
          {point_not_neighbor_of_plane=point_of_plane;}
          else
          {points_neighbor_of_plane.push_back(point_of_plane);}
        }
        //insert newly generated ear in maps
        vector<Size> new_possible_ear{temp_point,point_not_neighbor_of_plane,
            points_neighbor_of_plane[0],points_neighbor_of_plane[1]};
        double power=calc_power(new_possible_ear,curr_point);
        std::cout << "Inserting new ear" << std::endl;
        ears_map.insert(std::make_pair(new_possible_ear,power));
        rev_ears_map.insert(std::make_pair(power,new_possible_ear));
      }
    }
  }
}

/// Returns the total number of different points of two vectors
inline Size calculate_total_points(vector<Size> &v1, vector<Size> &v2)
{
  std::set<Size> union_of_points;
  for (Size temp_point: v1)
  {union_of_points.insert(temp_point);}
  for (Size temp_point: v2)
  {union_of_points.insert(temp_point);}
  return union_of_points.size();
}

///checks whether two triangles share a plane
inline bool triangle_shares_plane_with(vector<Size> &triangle1, vector<Size> &triangle2)
{
//if it shares a plane (and therefor size of union is (8-3)), return true
return (calculate_total_points(triangle1,triangle2)==5);
}

///checks whether 2 planes are equal
inline bool planes_are_equal(vector<Size> &plane1, vector<Size> &plane2)
{
return (calculate_total_points(plane1,plane2)==3);
}

/// Returns all relevant planes, given the tetrahedra we are currently adding
inline vector<vector<Size>> return_all_relevant_planes(vector<vector<Size>> &triangles_to_work_with)
{
  vector<vector<Size>> toreturn;
  for (vector<Size> triangle: triangles_to_work_with)
  {
    vector<Size> newplane1{triangle[0],triangle[1],triangle[2]};
    bool newpl1_already_in=false;//for keeping track whether newplane1 has already been encountered (and therefor NOT needs to be returned)
    vector<Size> newplane2{triangle[0],triangle[1],triangle[3]};
    bool newpl2_already_in=false;
    for (Size i=0; i<toreturn.size();)//check if planes are already in the vector
    {
      vector<Size> curr_plane=toreturn[i];
      if (planes_are_equal(newplane1,curr_plane))
      {
        newpl1_already_in=true;
        toreturn.erase(toreturn.begin()+i);
      }else if (planes_are_equal(newplane2,curr_plane))
      {
        newpl2_already_in=true;
        toreturn.erase(toreturn.begin()+i);
      }else{++i;}
    }
    if (!newpl1_already_in)//if we truly have a new plane, also put it in
    {toreturn.push_back(newplane1);}
    if (!newpl2_already_in)
    {toreturn.push_back(newplane2);}
  }
  return toreturn;
}



/// Coarsens the grid
///   @Param[in] perc_points_deleted: if the grid has not yet been coarsened to the next level,
/// then it determines the percentage of points deleted, otherwise does nothing
inline void Model :: coarsen_grid(const float perc_points_deleted)
{
  if (curr_coarsening_lvl==max_reached_coarsening_lvl)//if we truly need to refine the grid
  {
    if (max_reached_coarsening_lvl==0)//first time refining grid
      {
      current_nb_points=parameters.npoints();
      neighbors_lists.resize(2);
      neighbors_lists[0]=Neighbors(geometry.points.curr_neighbors);//TODO properly implement deep copy
      neighbors_lists[1]=Neighbors(geometry.points.curr_neighbors);//TODO properly implement deep copy
      for (int i=0; i<parameters.npoints(); i++)//calc 1-abs(d-d_other)/(d+d_other) for all points
        {
          double curr_diff=calc_diff_abundance_with_neighbours(i,1);
          density_diff_map.insert(std::pair<Size,double>(i,curr_diff));
          rev_density_diff_map.insert(std::pair<double,Size>(curr_diff,i));
        }
      }
      else
      {
        neighbors_lists.resize(max_reached_coarsening_lvl+1);
        neighbors_lists[max_reached_coarsening_lvl+1]=Neighbors(neighbors_lists[max_reached_coarsening_lvl]);//should be deep copy
      }
      max_reached_coarsening_lvl++;
      curr_coarsening_lvl++;
      //repeat n times:
      Size nb_points_to_remove=Size(perc_points_deleted*current_nb_points);
      //TODO: remove infinite loop in inner for loop
      for (Size i=0; i<nb_points_to_remove; i++)
      {
        if (i%10==0)
        {std::cout << "Ten points deleted" << std::endl;}
        Size curr_point=(*rev_density_diff_map.begin()).second;//the current point to remove
        vector<Size> neighbors_of_point=neighbors_lists[curr_coarsening_lvl].get_neighbors(curr_point);
        for (Size neighbor :neighbors_of_point)
        {//deleting the point from its neighbors, the point itself still has its neighbors (for now)
          neighbors_lists[curr_coarsening_lvl].delete_single_neighbor(neighbor,curr_point);
        }

        //first calculate the relevant neighbors of the neighbors of the 'deleted' point
        std::map<Size, std::set<Size>> neighbor_map;//stores the relevant neighbors; WARNING does not get updated while recalculating mesh around the current point
        // use only for all
        vector<Size> cpy_of_neighbors=vector<Size>(neighbors_of_point);
        std::sort(cpy_of_neighbors.begin(),cpy_of_neighbors.end());

        for (Size neighbor: neighbors_of_point)
        {//intersecting the neighbors of the neighbor with the neighbors
          vector<Size> cpy_of_curr_neighbors=vector<Size>(
            neighbors_lists[curr_coarsening_lvl].get_neighbors(neighbor));
          std::sort(cpy_of_curr_neighbors.begin(),cpy_of_curr_neighbors.end());
          std::set<Size> rel_neigbors_of_neighbor;//relevant neighbors
          std::set_intersection(cpy_of_neighbors.begin(),cpy_of_neighbors.end(),
                      cpy_of_curr_neighbors.begin(),cpy_of_curr_neighbors.end(),
                      std::inserter(rel_neigbors_of_neighbor,rel_neigbors_of_neighbor.begin()));
          neighbor_map.insert(std::make_pair(neighbor,rel_neigbors_of_neighbor));
        }
        //iterating over all lines neighbors_of_point[i],neighbors_of_point[j]
        std::multimap<vector<Size>,double> ears_map;
        std::multimap<double,vector<Size>> rev_ears_map;

        for (Size i=0; i<neighbors_of_point.size(); i++)
        {
          for (Size j=0; j<i; j++)
          {
            if (neighbor_map[neighbors_of_point[i]].count(j)==0)//if those two points are not yet neighbors
            {
              std::set<Size> temp_intersection;
              std::set_intersection(neighbor_map[neighbors_of_point[i]].begin(),neighbor_map[neighbors_of_point[i]].end(),
                            neighbor_map[neighbors_of_point[j]].begin(),neighbor_map[neighbors_of_point[j]].end(),
                            std::inserter(temp_intersection,temp_intersection.begin()));
              //then for every pair in the intersection of the neighbors of both points, check if they are neighbors
              //a better implementation is probably possible
              for (Size point1: temp_intersection)
              {
                for (Size point2: temp_intersection)
                {
                  if (point1<point2 && //just such that we do not have duplicates
                  neighbor_map[point1].count(point2)!=0)///if point1 and 2 are neighbors
                  {
                    vector<Size> new_triangle=vector<Size>{i,j,point1,point2};
                    double power=calc_power(new_triangle,curr_point);
                    ears_map.insert(std::make_pair(new_triangle,power));
                    rev_ears_map.insert(std::make_pair(power,new_triangle));
                    //invariant: the first two element of the vector should correspond to the neighbors we want to add
                  }
                }
              }
            }
          }
        }
        std::cout << "number of neigbors: " << neighbors_of_point.size() << std::endl;
        std::cout << "ears_map size: " << ears_map.size() << std::endl;
        //now that we have all initial 'triangles', we can finally start to add them
        //infinite loop!!!! probably from deleting and readding some times
        //while(ears_map.size()>7)//in the last iteration, you should have an octohedron left: -> 8 possible ears
        while(!ears_map.empty())
        {
          std::cout << "Looping: current size of ears map: " << ears_map.size() << std::endl;
          vector<Size> triangle=(*rev_ears_map.begin()).second;//triangle which we are currently adding to the mesh
          vector<vector<Size>> triangles_to_work_with;//vector of triangles from which we will generate new triangles
          triangles_to_work_with.push_back(triangle);
          //if (!vector_contains_element(neighbors_lists[curr_coarsening_lvl].get_neighbors(triangle[0]),triangle[1]))
          //{//necessary because of 'hack' with setting absolute priority to ear which also tries to add the same neighbors
          neighbors_lists[curr_coarsening_lvl].add_single_neighbor(triangle[0], triangle[1]);
          neighbors_lists[curr_coarsening_lvl].add_single_neighbor(triangle[1], triangle[0]);
          neighbor_map[triangle[0]].insert(triangle[1]);
          neighbor_map[triangle[1]].insert(triangle[0]);
          //invariant: the first two element of the vector should correspond to the neighbors we want to add
          //}else{std::cout << "hack active" << std::endl;}
          delete_vector_from_both_maps(ears_map,rev_ears_map,triangle);
          //check which possible ears need to be deleted and which need to be added to the list from which to generate new ears
          for (auto it = ears_map.cbegin(); it != ears_map.cend();)
          {
            vector<Size> curr_triangle=(*it).first;//candidate triangle to check whether it needs to be changed/deleted
            //if this candidate triangle would try to add the same neighbors, we just delete it (although i'm not sure whether it can happen before the mesh is almost complete...)
            if ((curr_triangle[0]==triangle[0]&&curr_triangle[1]==triangle[1])||
                (curr_triangle[1]==triangle[0]&&curr_triangle[0]==triangle[1]))
            {
              triangles_to_work_with.push_back(curr_triangle);
              it=delete_vector_from_both_maps(ears_map,rev_ears_map,curr_triangle);
                // double curr_power=(*it).second;
                // if (curr_power!=-DBL_MAX)//infinite loop protection
                // {
                // std::cout << "Infinite loop prevention activated" << std::endl;
                // it = delete_vector_from_both_maps(ears_map, rev_ears_map, curr_triangle);
                // ears_map.insert(std::make_pair(curr_triangle,-DBL_MAX));
                // rev_ears_map.insert(std::make_pair(-DBL_MAX,curr_triangle));//giving maximum priority to this triangle
//              }//a better solution would be to also recalulate all other triangles which share a plane, but this requires some refactoring
//              else
//              {
//                ++it;
//              }
            }
            else
            {
              if (triangle_shares_plane_with(triangle, curr_triangle))
                // std::set<Size> union_of_points;
                // //insert all points in a set; easier to work with
                // for (Size temp_point: triangle)
                // {
                //   union_of_points.insert(temp_point);
                // }
                // for (Size temp_point: curr_triangle)
                // {
                //   union_of_points.insert(temp_point);
                // }
                // //if it shares a plane (and therefor size of union is (8-3)), delete them
                // if (union_of_points.size()==5)
              {
                it = delete_vector_from_both_maps(ears_map, rev_ears_map, curr_triangle);
              }
              else//we do nothing with the current entry
              {
                ++it;
              }
            }
          }

            // vector<Size> plane1{triangle[0],triangle[1],triangle[2]};
            // vector<Size> plane2{triangle[0],triangle[1],triangle[3]};
            // Size count_neighbours_pl1;//counts the number of neighbors with plane 1
            // Size count_neighbours_pl2;//counts the number of neighbors with plane 2
            // Size point_not_neighbors_pl1;//just to not need to recalc which point of pl1 is not a neighbor
            // Size point_not_neighbors_pl2;

          vector<vector<Size>> relevant_planes=return_all_relevant_planes(triangles_to_work_with);
          //TODO FIXME: check if the same ear can be generated twice and check if it would actually matter...
          for (vector<Size> plane: relevant_planes)
          {
            generate_new_ears(neighbors_of_point, plane, neighbor_map,
                ears_map, rev_ears_map, curr_point);
          }
            //   //maybe put this into a seperate function
            //   //now that we have deleted all redundant ears, it is time to readd some new ears based on the new planes created
            // for (Size temp_point: neighbors_of_point)
            // { //we want to create new triangles, not just the current triangle again...
            //   if (!vector_contains_element(triangle,temp_point))
            //   {
            //     //we are currently adding far too many new ears (sometimes)??
            //     count_neighbours_pl1=0;
            //     count_neighbours_pl2=0;
            //     //vector<Size> neighbor_list_of_temp_point=neighbors_lists[curr_coarsening_lvl].get_neighbors(temp_point);
            //     for (Size point_of_pl1: plane1)
            //     {
            //       if (neighbor_map[temp_point].count(point_of_pl1)!=0)//vector_contains_element(neighbor_list_of_temp_point,point_of_pl1))
            //       {count_neighbours_pl1++;}
            //     }
            //     for (Size point_of_pl2: plane2)
            //     {
            //       if (neighbor_map[temp_point].count(point_of_pl2)!=0)//vector_contains_element(neighbor_list_of_temp_point,point_of_pl2))
            //       {count_neighbours_pl2++;}
            //     }
            //     //if the candidate point is good for creating an ear with plane1
            //     if (count_neighbours_pl1==2)
            //     {
            //       Size point_not_neighbor_of_plane;//the point which is not a neighbor
            //       vector<Size> points_neighbor_of_plane;
            //       points_neighbor_of_plane.reserve(2);
            //       for (Size point_of_pl1: plane1)
            //       {
            //         if (neighbor_map[temp_point].count(point_of_pl1)==0)//!vector_contains_element(neighbor_list_of_temp_point,point_of_pl1))
            //         {point_not_neighbor_of_plane=point_of_pl1;}
            //         else
            //         {points_neighbor_of_plane.push_back(point_of_pl1);}
            //       }
            //       //insert newly generated ear in maps
            //       vector<Size> new_possible_ear{temp_point,point_not_neighbor_of_plane,
            //           points_neighbor_of_plane[0],points_neighbor_of_plane[1]};
            //       double power=calc_power(new_possible_ear,curr_point);//TODO calculate this!!!!!
            //       std::cout << "Inserting new ear 1" << std::endl;
            //       ears_map.insert(std::make_pair(new_possible_ear,power));
            //       rev_ears_map.insert(std::make_pair(power,new_possible_ear));
            //       }
            //
            //     //if the candidate point is good for creating an ear with plane1
            //     if (count_neighbours_pl2==2)
            //     {
            //       Size point_not_neighbor_of_plane;//the point which is not a neighbor
            //       vector<Size> points_neighbor_of_plane;
            //       points_neighbor_of_plane.reserve(2);
            //       for (Size point_of_pl2: plane2)
            //       {
            //         if (neighbor_map[temp_point].count(point_of_pl2)==0)//!vector_contains_element(neighbor_list_of_temp_point,point_of_pl2))
            //         {point_not_neighbor_of_plane=point_of_pl2;}
            //         else
            //         {points_neighbor_of_plane.push_back(point_of_pl2);}
            //       }
            //       //insert newly generated ear in maps
            //       vector<Size> new_possible_ear{temp_point,point_not_neighbor_of_plane,
            //           points_neighbor_of_plane[0],points_neighbor_of_plane[1]};
            //       double power=calc_power(new_possible_ear,curr_point);
            //       std::cout << "Inserting new ear 2" << std::endl;
            //       ears_map.insert(std::make_pair(new_possible_ear,power));
            //       rev_ears_map.insert(std::make_pair(power,new_possible_ear));
            //       }
            //
            //     }
            //   }

        }

        //recalc density diff of neighbors (due to new mesh (and ofc deletion of point))
        for (Size neighbor :neighbors_of_point)
        {
          delete_int_from_both_maps(density_diff_map,rev_density_diff_map,neighbor);

          double calc_diff=calc_diff_abundance_with_neighbours(neighbor,curr_coarsening_lvl);
          density_diff_map.insert(std::pair<Size,double>(neighbor,calc_diff));
          rev_density_diff_map.insert(std::pair<double,Size>(calc_diff,neighbor));
        }

      //finally, delete the neighbors of this point
      neighbors_lists[curr_coarsening_lvl].delete_all_neighbors(curr_point);
      //and remove this entry from density maps
      delete_int_from_both_maps(density_diff_map, rev_density_diff_map, curr_point);
    }
    current_nb_points-=nb_points_to_remove;//decrement the number of points remaining
  }
  else
  {//TODO check whether deep copy!!!!!!
    curr_coarsening_lvl++;
    geometry.points.curr_neighbors=neighbors_lists[curr_coarsening_lvl];//should be shallow copy
  }
}





/// Resets the grid
inline void Model :: reset_grid()
{
  curr_coarsening_lvl=0;
  max_reached_coarsening_lvl=0;
  density_diff_map=std::multimap<Size,double>();
  rev_density_diff_map=std::multimap<double,Size>();
  geometry.points.curr_neighbors=Neighbors(neighbors_lists[0]);
  neighbors_lists=vector<Neighbors>();
  current_nb_points=0;
}
