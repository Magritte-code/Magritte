#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cfloat>
#include <set>
#include "tools/types.hpp"
#include <limits>
#include <algorithm>

//NOTE: i have mistakenly called tetrahedra triangles throughout this entire piece of code

// Calculates the relative difference of the density with respect to the neighbours
//    @Returns: abs((d-d_other)/(d+d_other)) such that the smallest values will be assigned
// to the points we want to delete
///TODO : change to use only two points
// inline double Model :: calc_diff_abundance_with_neighbours(Size point, Size next_coars_lvl)
// {
//   double temp_max=0;//current implementation is probably slow
//   for (Size neigbor: neighbors_lists[next_coars_lvl].get_neighbors(point))
//     {
//       double abundance_self=chemistry.species.abundance[point][1];
//       double abundance_other=chemistry.species.abundance[neigbor][1];
//       //std::cout<< "abundance self: " << abundance_self << std::endl;
//       //std::cout<< "abundance other: " << abundance_other << std::endl;
//       double temp_diff=std::abs((abundance_self-abundance_other)/(abundance_self+abundance_other));
//       temp_max=std::max(temp_max, temp_diff);
//     }
//     return temp_max;
// }

// Calculates the relative difference of the density of point1 with respect to point2
//    @Returns: abs((d-d_other)/(d+d_other)) such that the smallest values will be assigned
inline double Model :: calc_diff_abundance_with_point(Size point1, Size point2)
{
    double abundance_self=chemistry.species.abundance[point1][1];
    double abundance_other=chemistry.species.abundance[point2][1];
    //std::cout<< "abundance self: " << abundance_self << std::endl;
    //std::cout<< "abundance other: " << abundance_other << std::endl;
    double temp_diff=std::abs((abundance_self-abundance_other)/(abundance_self+abundance_other));
    return temp_diff;
}

/// Returns function that determines whether two points are similar enough, according to the tolerance
inline std::function<bool(Size,Size)> Model::points_are_similar(double tolerance)
{
  std::function<bool(Size,Size)> fun_to_return = [&](Size p1, Size p2)
  {
    return (calc_diff_abundance_with_point(p1,p2)<tolerance);
  };
  return fun_to_return;
}


/// Deprecated
// ///Helper function for deleting stuff from the density_diff_maps
// ///   @Param[in/out] std_map: the map with points as keys
// ///   @Param[in/out] reverse_map: the map with points as values
// ///   @Param[in] point_to_del: the point to delete from both maps
// inline void delete_int_from_both_maps(std::multimap<Size,double> &std_map, std::multimap<double,Size> &reverse_map, Size point_to_del)
// {
//   double curr_value=(*std_map.find(point_to_del)).second;
//   std_map.erase(point_to_del);//key is unique, so no problem
//   //there are multiple points that currently have the same diff value, so we must iterate
//   typedef std::multimap<double,Size>::iterator iterator;
//   std::pair<iterator, iterator> iterpair = reverse_map.equal_range(curr_value);
//   iterator it = iterpair.first;
//   for (; it != iterpair.second; ++it) {
//     if (it->second == point_to_del) {
//       reverse_map.erase(it);
//       break;
//     }
//   }
// }
//
// ///Helper function for deleting stuff from the ears_maps
// ///   @Param[in/out] std_map: the map with vectors as keys
// ///   @Param[in/out] reverse_map: the map with vectors as values
// ///   @Param[in] vect_to_del: the vector to delete from both maps
// ///   @Returns: Returns an iterator for the standard map (needed for being able to delete while iterating)
// inline std::multimap<vector<Size>,double>::iterator delete_vector_from_both_maps(
//   std::multimap<vector<Size>,double> &std_map, std::multimap<double,vector<Size>> &reverse_map, vector<Size> vect_to_del)
// {
//   double curr_value=(*std_map.find(vect_to_del)).second;
//   std::multimap<vector<Size>,double>::iterator it2=std_map.erase(std_map.find(vect_to_del));//key is unique, so no problem
//   //there are multiple points that currently have the same diff value, so we must iterate
//   typedef std::multimap<double,vector<Size>>::iterator iterator;
//   std::pair<iterator, iterator> iterpair = reverse_map.equal_range(curr_value);
//   iterator it = iterpair.first;
//   for (; it != iterpair.second; ++it) {
//     if (it->second == vect_to_del) {
//       reverse_map.erase(it);
//       break;
//     }
//   }
//   return it2;
// }
//
// ///Helper function: checks whether a vector contains an element
// ///   @Param[in] vect: the vector to check for
// ///   @Param[in] element: the element to check for
// ///   @Returns bool: true if element in vect, false otherwise
// inline bool vector_contains_element(const vector<Size> &vect, Size element)
// {
//   return (std::find(vect.begin(), vect.end(), element) != vect.end());
// }
//
//
// /// Given a tetrahedron, this calculates the power of the point with respect to the circumsphere
// /// @Parameter [in] triangle: the tetrahedron from which we use the circumsphere
// /// @Parameter [in] point: the point the calculate the power with
// /// returns a positive value if the point is inside the circumsphere and a negative value if the point lies outside the circumsphere (zero if on the circumsphere)
// /// returns nan when the points of the 'tetrahedron' are (almost) coplanar
// inline double Model :: calc_power(const vector<Size> &triangle, Size point){
//   Vector3D pos1=geometry.points.position[triangle[0]];
//   Vector3D pos2=geometry.points.position[triangle[1]];
//   Vector3D pos3=geometry.points.position[triangle[2]];
//   Vector3D pos4=geometry.points.position[triangle[3]];
//   Vector3D posp=geometry.points.position[point];//position of point
//   //dividing insphere test with orientation test
//
//   double average_dist_sq=((pos1-posp).squaredNorm()+(pos2-posp).squaredNorm()+(pos3-posp).squaredNorm()+(pos4-posp).squaredNorm())/4;
//   //std::cout << "Average distance squared: " << average_dist_sq << std::endl;
//
//   Eigen::Matrix<double,4,4> orient;
//   orient << pos1.x(),pos2.x(),pos3.x(),pos4.x(),
//             pos1.y(),pos2.y(),pos3.y(),pos4.y(),
//             pos1.z(),pos2.z(),pos3.z(),pos4.z(),
//             1,1,1,1;
//   //result is nan or stupidly large when the tetrahedra is (almost) coplanar
//   Eigen::Matrix<double,5,5> insphere;
//   insphere << pos1.x(),pos2.x(),pos3.x(),pos4.x(),posp.x(),
//               pos1.y(),pos2.y(),pos3.y(),pos4.y(),posp.y(),
//               pos1.z(),pos2.z(),pos3.z(),pos4.z(),posp.z(),
//               pos1.squaredNorm(),pos2.squaredNorm(),pos3.squaredNorm(),pos4.squaredNorm(),posp.squaredNorm(),
//               1,1,1,1,1;
//   //dividing insphere test with orientation test
//   double result=insphere.determinant()/orient.determinant();
//   //TODO: change this test to something reasonable instead (eg divide the result by the average distance squared to point and change to pow(10,5))
//   if (std::isnan(result)||std::abs(result)/average_dist_sq>pow(10,5))//10^5 is chosen arbitrarily
//   {
//     //std::cout << "tetradron is coplanar; result = " << result/average_dist_sq << std::endl;
//     //std::cout << "orient determinant = " << orient.determinant() << std::endl;
//     return std::numeric_limits<double>::quiet_NaN();
//   }
//   else
//   {
//     //std::cout << "reasonable tetrahedron; result = " << result/average_dist_sq << std::endl;
//     //std::cout << "orient determinant = " << orient.determinant() << std::endl;
//     return result;
//   }
// }
//
// //TODO: maybe just return the sign (because possible overflow issues if the number get much larger)
// inline double Model::orientation(const vector<Size> &plane, Size point){
//   Vector3D pos1=geometry.points.position[plane[0]];
//   Vector3D pos2=geometry.points.position[plane[1]];
//   Vector3D pos3=geometry.points.position[plane[2]];
//   Vector3D posp=geometry.points.position[point];//position of point
//
//   Eigen::Matrix<double,4,4> orient;
//   orient << pos1.x(),pos2.x(),pos3.x(),posp.x(),
//             pos1.y(),pos2.y(),pos3.y(),posp.y(),
//             pos1.z(),pos2.z(),pos3.z(),posp.z(),
//             1,1,1,1;
//   return orient.determinant();
// }
//
//
// // /// Returns the total number of different points of two vectors
// // inline Size calculate_total_points(vector<Size> &v1, vector<Size> &v2)
// // {
// //   std::set<Size> union_of_points;
// //   for (Size temp_point: v1)
// //   {union_of_points.insert(temp_point);}
// //   for (Size temp_point: v2)
// //   {union_of_points.insert(temp_point);}
// //   return union_of_points.size();
// // }
//
//
// /// Returns the total number of different points of two vectors //attempt 2: now no construction of sets needed (should be faster for small vectors)
// inline Size calculate_total_points(vector<Size> &v1, vector<Size> &v2)
// {
//   Size nb_elements_in_v2_not_in_v1=0;
//   for (Size element:v2)
//   {
//     if (!vector_contains_element(v1,element))
//     {
//       nb_elements_in_v2_not_in_v1++;
//     }
//   }
//   return v1.size()+nb_elements_in_v2_not_in_v1;
// }
//
// /// Checks whether a tetrahedron can be formed given 2 lines, by drawing exactly one extra line
// ///   TODO!!!
// ///   @Returns: tetrahedron if a tetrahedron can be formed, else returns empty vector.
// inline vector<Size> lines_form_tetrahedron(vector<Size> &v1, vector<Size> &v2, std::map<Size, std::set<Size>> &neighbor_map)
// {//if we have 4 distinct points
//   if (calculate_total_points(v1,v2)==4)
//   {//check if the total number of neighbors between the two lines is equal to 3
//     Size nb_count=0;
//     Size idx_v1_not_neighbor;//keeps track of which index of v1 has no neighbor with a index in v2
//     Size idx_v2_not_neighbor;//same, but vice versa
//     Size idx_count=0;
//     for (Size point2: v2)
//     {
//       if (neighbor_map[v1[0]].count(point2)!=0)
//       {nb_count++;idx_count++;}
//       else{
//       idx_v1_not_neighbor=0;
//       idx_v2_not_neighbor=idx_count;
//       }
//     }
//     idx_count=0;
//     for (Size point2: v2)
//     {
//       if (neighbor_map[v1[1]].count(point2)!=0)
//       {nb_count++;idx_count++;}
//       else{
//         idx_v1_not_neighbor=1;
//         idx_v2_not_neighbor=idx_count;
//       }
//     }
//     if (nb_count==3)
//     {
//       return vector<Size>{v1[idx_v1_not_neighbor],v2[idx_v2_not_neighbor],
//         v1[1-idx_v1_not_neighbor],v2[1-idx_v2_not_neighbor]};
//     }
//   }
//   return vector<Size>();
// }
//
//
// //Checks whether a proposed ear contains a forbidden plane
// inline bool has_forbidden_plane(vector<Size> &ear,std::set<std::vector<Size>> &forbidden_planes)
// {
//   for (vector<Size> forbidden_plane:forbidden_planes)
//   {//if triangle has forbidden plane
//     if (calculate_total_points(forbidden_plane,ear)==4)
//     {
//       // std::cout << "Forbidden plane detected" << std::endl;
//       return true;
//     }
//   }
//   return false;
// }
//
// //checks if a plane is an edge plane
// inline bool is_edge_plane(vector<Size> &plane,std::set<std::vector<Size>> &edge_planes)
// {
//   for (vector<Size> edge_plane:edge_planes)
//   {//if plane is edge_plane
//     if (calculate_total_points(edge_plane,plane)==3)
//     {
//       return true;
//     }
//   }
//   return false;
// }
// //checks if plane is a forbidden plane; does the same thing as is_edge_plane
// inline bool is_forbidden_plane(vector<Size> &plane,std::set<std::vector<Size>> &forbidden_planes)
// {return is_edge_plane(plane,forbidden_planes);}
//
// //note: it should be possible to cache some results
// // Checks whether a point is surrounded by tetrahedra //attempt 2: changed order of if statements to reduce computational overhead
// inline bool point_surrounded_by_tetras(Size &point, std::map<vector<Size>,Size> &plane_counter, std::set<std::vector<Size>> &edge_planes)
// {
//   for (auto it=plane_counter.begin(); it!=plane_counter.end(); ++it)
//   {
//     vector<Size> curr_plane=(*it).first;
//     Size nb_tetra_of_plane=(*it).second;
//     if (vector_contains_element(curr_plane, point))//check if the plane we are currently looking at
//     {
//       if (nb_tetra_of_plane<2)
//       {
//       if (!is_edge_plane(curr_plane,edge_planes))
//       {
//         // std::cout<<"Inner plane not fully surrounded: "<<curr_plane[0]<<", "<<curr_plane[1]<<", "<<curr_plane[2]<<std::endl;
//         return false;
//       }else{
//         if (nb_tetra_of_plane<1)
//         { //std::cout<<"Outer plane not fully surrounded: "<<curr_plane[0]<<", "<<curr_plane[1]<<", "<<curr_plane[2]<<std::endl;
//           return false;}
//       }
//       }
//       // if (!is_edge_plane(curr_plane,edge_planes))
//       // {
//       //   if (nb_tetra_of_plane<2)
//       //   { std::cout<<"Inner plane not fully surrounded: "<<curr_plane[0]<<", "<<curr_plane[1]<<", "<<curr_plane[2]<<std::endl;
//       //     return false;}
//       // }else{
//       //   if (nb_tetra_of_plane<1)
//       //   { std::cout<<"Outer plane not fully surrounded: "<<curr_plane[0]<<", "<<curr_plane[1]<<", "<<curr_plane[2]<<std::endl;
//       //     return false;}
//       // }
//     }
//   }
//   return true;
// }
//
// // Adds a single plane to plane_counter
// inline void add_plane_to_counter(vector<Size> &plane, std::map<vector<Size>,Size> &plane_counter)
// {
//   std::sort(plane.begin(),plane.end());//order of elements in plane can change!
//   if (plane_counter.find(plane)!=plane_counter.end()) {
//     plane_counter[plane]++;
//   } else {
//     plane_counter.insert(std::pair<vector<Size>,Size>(plane,1));
//   }
// }
//
// // Adds the planes of a tetrahedron to the planes counter
// inline void add_planes_of_tetra_to_counter(vector<Size> &tetra, std::map<vector<Size>,Size> &plane_counter)
// {
//   vector<Size> plane1{tetra[0],tetra[1],tetra[2]};
//   vector<Size> plane2{tetra[0],tetra[1],tetra[3]};
//   vector<Size> plane3{tetra[0],tetra[2],tetra[3]};
//   vector<Size> plane4{tetra[1],tetra[2],tetra[3]};
//   add_plane_to_counter(plane1, plane_counter);
//   add_plane_to_counter(plane2, plane_counter);
//   add_plane_to_counter(plane3, plane_counter);
//   add_plane_to_counter(plane4, plane_counter);
// }
//
// //generates the new ears and inserts them into the ears maps
// //TODO: fix this, wrong logic!!!!!!!!!
// // inline void Model::generate_new_ears(const std::set<vector<Size>> &neighbor_lines, vector<Size> &new_line, std::map<Size, std::set<Size>> &neighbor_map,
// //  std::multimap<vector<Size>,double> &ears_map, std::multimap<double,vector<Size>> &rev_ears_map, Size &curr_point)
// inline void Model::generate_new_ears(const vector<Size> &neighbors_of_point, const vector<Size> &plane, std::map<Size, std::set<Size>> &neighbor_map,
//  std::multimap<vector<Size>,double> &ears_map, std::multimap<double,vector<Size>> &rev_ears_map, Size &curr_point, Size &orient_point, std::map<vector<Size>,Size> &plane_counter, std::set<std::vector<Size>> &edge_planes, std::set<vector<Size>> forbidden_planes)
// {
//   // std::cout << "plane is: " << plane[0]<< ", " << plane[1] << ", " << plane[2] << std::endl;
//   // std::cout << "orient point is: " << orient_point << std::endl;
// //   std::cout << "line is: " << new_line[0]<< ", " << new_line[1] << std::endl;
//   Size count_neighbours;
//
//   for (Size temp_point: neighbors_of_point)
//   { //we want to create new triangles, therefore a point in the plane is not useful
//     if (!vector_contains_element(plane,temp_point))
//     {
//       //we are currently adding sometimes far too many new ears??
//       count_neighbours=0;
//       for (Size point_of_plane: plane)
//       {
//         if (neighbor_map[temp_point].count(point_of_plane)!=0)//if temp_point is neighbor of point_of_plane
//         {count_neighbours++;}
//       }
//
//       //if the candidate point is good for creating an ear with plane
//       if (count_neighbours==2)
//       {
//         double prod_of_orient=orientation(plane,temp_point)*orientation(plane,orient_point);//should be negative if temp_point lies on the opposite side of the plane as the orientation point
//         if (prod_of_orient<0)
//         {
//           // std::cout << "neighbors of temp_point are: ";
//           // for (std::set<Size>::iterator it=neighbor_map[temp_point].begin(); it!=neighbor_map[temp_point].end(); ++it)
//           //   std::cout << ' ' << *it;
//           // std::cout << std::endl;
//           // std::cout << "temp_point is: " << temp_point << std::endl;
//
//           Size point_not_neighbor_of_plane;//the point of the plane which is not a neighbor
//           vector<Size> points_neighbor_of_plane;
//           for (Size point_of_plane: plane)
//           {
//             if (neighbor_map[temp_point].count(point_of_plane)==0)
//             {point_not_neighbor_of_plane=point_of_plane;}
//             else
//             {points_neighbor_of_plane.push_back(point_of_plane);}
//           }
//           //insert newly generated ear in maps; as usual, the first two points denote the proposed line
//           vector<Size> new_possible_ear{temp_point,point_not_neighbor_of_plane,
//               points_neighbor_of_plane[0],points_neighbor_of_plane[1]};
//           //checks whether we either have a fobidden plane, or that we do not need to add the line at all
//           if (!has_forbidden_plane(new_possible_ear,forbidden_planes)&&((!point_surrounded_by_tetras(temp_point,plane_counter,edge_planes))
//           &&(!point_surrounded_by_tetras(point_not_neighbor_of_plane,plane_counter,edge_planes))))
//           {
//             double power=calc_power(new_possible_ear,curr_point);
//             if (!std::isnan(power))
//             {//otherwise we would be proposing a coplanar tetrahedron, which would be ridiculous
//               // std::cout << "Inserting new ear; power is: " << power << std::endl;
//               ears_map.insert(std::make_pair(new_possible_ear,power));
//               rev_ears_map.insert(std::make_pair(power,new_possible_ear));
//             }
//           }
//         }
//       }
//     }
//   }
// }
//
// ///checks whether two triangles share a plane
// inline bool triangle_shares_plane_with(vector<Size> &triangle1, vector<Size> &triangle2)
// {
// //if it shares a plane (and therefor size of union is (8-3)), return true
// return (calculate_total_points(triangle1,triangle2)==5);
// }
//
// ///checks if two tetrahedra are equal
// inline bool triangles_are_equal(vector<Size> &triangle1, vector<Size> &triangle2)
// {
// return (calculate_total_points(triangle1,triangle2)==4);
// }
//
// ///checks whether 2 planes are equal
// inline bool planes_are_equal(vector<Size> &plane1, vector<Size> &plane2)
// {
// return (calculate_total_points(plane1,plane2)==3);
// }
//
// /// Returns all relevant planes (and their orientation points), given the tetrahedra we are currently adding
// inline std::pair<vector<Size>,vector<vector<Size>>> return_all_relevant_planes(vector<vector<Size>> &triangles_to_work_with, std::set<std::vector<Size>> &forbidden_planes)
// {
//   vector<Size> orient_points;//orientation points to make sure that it expands inwards
//   vector<vector<Size>> toreturn;
//   for (vector<Size> triangle: triangles_to_work_with)
//   {
//     vector<Size> newplane1{triangle[0],triangle[1],triangle[2]};
//     bool newpl1_already_in=false;//for keeping track whether newplane1 has already been encountered (and therefor NOT needs to be returned)
//     vector<Size> newplane2{triangle[0],triangle[1],triangle[3]};
//     bool newpl2_already_in=false;
//     for (Size i=0; i<toreturn.size();)//check if planes are already in the vector
//     {
//       vector<Size> curr_plane=toreturn[i];
//       if (planes_are_equal(newplane1,curr_plane))
//       {
//         newpl1_already_in=true;
//         // forbidden_planes.insert(newplane1);
//         toreturn.erase(toreturn.begin()+i);
//         orient_points.erase(orient_points.begin()+i);
//       }else if (planes_are_equal(newplane2,curr_plane))
//       {
//         newpl2_already_in=true;
//         // forbidden_planes.insert(newplane2);
//         toreturn.erase(toreturn.begin()+i);
//         orient_points.erase(orient_points.begin()+i);
//       }else{++i;}
//     }
//     if (!newpl1_already_in)//if we truly have a new plane, add it and the orientation point
//     {toreturn.push_back(newplane1);orient_points.push_back(triangle[3]);}
//     if (!newpl2_already_in)
//     {toreturn.push_back(newplane2);orient_points.push_back(triangle[2]);}
//   }
//
//
//
//   return std::pair<vector<Size>,vector<vector<Size>>>(orient_points,toreturn);
// }
//
// /// Coarsens the grid
// ///   @Param[in] perc_points_deleted: if the grid has not yet been coarsened to the next level,
// /// then it determines the percentage of points deleted, otherwise does nothing
// inline void Model :: coarsen_grid(float perc_points_deleted)
// {
//   if (curr_coarsening_lvl==max_reached_coarsening_lvl)//if we truly need to refine the grid
//   {
//     if (max_reached_coarsening_lvl==0)//first time refining grid; thus setting heaps of values
//       {
//       current_nb_points=parameters.npoints();
//       neighbors_lists.resize(2);
//       neighbors_lists[0]=Neighbors(geometry.points.curr_neighbors);
//       neighbors_lists[1]=Neighbors(geometry.points.curr_neighbors);
//       mask_list.resize(2);
//       mask_list[0]=std::vector<bool>(parameters.npoints(),true);
//       mask_list[1]=std::vector<bool>(parameters.npoints(),true);
//       nb_points_at_lvl.push_back(parameters.npoints());
//       for (Size i=0; i<parameters.npoints(); i++)//calc 1-abs(d-d_other)/(d+d_other) for all points
//         {//we will NEVER delete boundary points, so we do not add them to the queue
//           if(geometry.not_on_boundary(i))
//           {
//           double curr_diff=calc_diff_abundance_with_neighbours(i,1);
//           density_diff_map.insert(std::pair<Size,double>(i,curr_diff));
//           rev_density_diff_map.insert(std::pair<double,Size>(curr_diff,i));
//           //std::cout << "point: " << i << " density diff: " << curr_diff<< std::endl;
//           }
//         }
//       }
//       else
//       {
//         neighbors_lists.resize(max_reached_coarsening_lvl+1);
//         neighbors_lists[max_reached_coarsening_lvl+1]=Neighbors(neighbors_lists[max_reached_coarsening_lvl]);//should be deep copy
//         mask_list.resize(max_reached_coarsening_lvl+1);
//         mask_list[max_reached_coarsening_lvl+1]=std::vector<bool>(mask_list[max_reached_coarsening_lvl]);//should be a deep copy
//       }
//       max_reached_coarsening_lvl++;
//       curr_coarsening_lvl++;
//       //repeat n times:
//       Size nb_points_to_remove=Size(perc_points_deleted*current_nb_points);//TODO: change to curr nb - nb boundary
//       //debug vars
//       reduced_neighbors_before.resize(nb_points_to_remove);
//       reduced_neighbors_after.resize(nb_points_to_remove);
//       deleted_points.resize(nb_points_to_remove);
//       added_lines.resize(nb_points_to_remove);
//       added_tetras.resize(nb_points_to_remove);
//
//       for (Size i=0; i<nb_points_to_remove; i++)
//       {
//         if (rev_density_diff_map.empty())
//         {
//           // std::cout << "You are trying to delete more points than possible; stuff will break" << std::endl;
//           break;
//         }
//         Size curr_point=(*rev_density_diff_map.begin()).second;//the current point to remove
//         mask_list[curr_coarsening_lvl][curr_point]=false;
//         //std::cout << "current density diff: " << (*rev_density_diff_map.begin()).first << std::endl;
//         // std::cout << "current point: " << curr_point << std::endl;
//
//         vector<Size> neighbors_of_point=neighbors_lists[curr_coarsening_lvl].get_neighbors(curr_point);
//         // std::cout << "neighbors of current point: ";
//         // for (auto it = neighbors_of_point.begin(); it!=neighbors_of_point.end(); it++)
//         // {
//         //   std::cout << *it << " ";
//         // }
//         // std::cout << std::endl;
//
//         for (Size neighbor :neighbors_of_point)
//         {//deleting the point from its neighbors, the point itself still has its neighbors (for now)
//           neighbors_lists[curr_coarsening_lvl].delete_single_neighbor(neighbor,curr_point);
//         }
//
//         //first calculate the relevant neighbors of the neighbors of the 'deleted' point
//         std::map<Size, std::set<Size>> neighbor_map;//stores the relevant neighbors;
//         vector<Size> cpy_of_neighbors=vector<Size>(neighbors_of_point);
//         std::sort(cpy_of_neighbors.begin(),cpy_of_neighbors.end());
//
//
//         for (Size neighbor: neighbors_of_point)
//         {//intersecting the neighbors of the neighbor with the neighbors
//           vector<Size> cpy_of_curr_neighbors=vector<Size>(
//             neighbors_lists[curr_coarsening_lvl].get_neighbors(neighbor));
//           std::sort(cpy_of_curr_neighbors.begin(),cpy_of_curr_neighbors.end());
//           std::set<Size> rel_neigbors_of_neighbor;//relevant neighbors
//           std::set_intersection(cpy_of_neighbors.begin(),cpy_of_neighbors.end(),
//                       cpy_of_curr_neighbors.begin(),cpy_of_curr_neighbors.end(),
//                       std::inserter(rel_neigbors_of_neighbor,rel_neigbors_of_neighbor.begin()));
//           neighbor_map.insert(std::make_pair(neighbor,rel_neigbors_of_neighbor));
//         }
//
//         //debug stuff
//         reduced_neighbors_before[i]=std::map<Size, std::set<Size>>(neighbor_map);
//         deleted_points[i]=curr_point;
//         //end debug stuff
//
//
//         //iterating over all lines neighbors_of_point[i],neighbors_of_point[j]
//         std::multimap<vector<Size>,double> ears_map;
//         std::multimap<double,vector<Size>> rev_ears_map;
//
//
//         std::set<vector<Size>> all_added_tetras;
//         std::set<vector<Size>> edge_planes;//the planes which lie on the edge of the relevant domain
//         std::set<vector<Size>> forbidden_planes;//contains the planes that should not be used during initial generation
//         std::map<vector<Size>,Size> plane_counter;//counts the number of tetrahedra each plane has
//         // finding all forbidden planes
//         for (Size pointi:neighbors_of_point)
//         {
//           for (Size pointj:neighbors_of_point)
//           {
//             if (pointi<pointj && neighbor_map[pointi].count(pointj)!=0)//if those two points are neighbors
//             {
//               std::set<Size> temp_intersection;
//               std::set_intersection(neighbor_map[pointi].begin(),neighbor_map[pointi].end(),
//                             neighbor_map[pointj].begin(),neighbor_map[pointj].end(),
//                             std::inserter(temp_intersection,temp_intersection.begin()));
//               //then for every pair in the intersection of the neighbors of both points, check if they are neighbors
//               //a better implementation is probably possible
//
//               for (Size point1: temp_intersection)
//               {
//                 if (pointj<point1)//enforcing strict order
//                 {
//                   vector<Size> curr_plane{pointi,pointj,point1};
//                   edge_planes.insert(curr_plane);//add plane to edge planes
//                   // std::cout << "Edge plane found: " << curr_plane[0] << "," << curr_plane[1] << "," << curr_plane[2] << std::endl;
//                   for (Size point2: temp_intersection)
//                   {//adding strict order in order to remove duplicates
//                     if (point1<point2 && //just such that we do not have duplicates
//                     neighbor_map[point1].count(point2)!=0)///if point1 and 2 are neighbors
//                     {
//                       // std::cout << "Found initial tetrahedron" << std::endl;
//                       vector<vector<Size>> tetra_planes;
//                       vector<Size> orient_points{pointi,pointj,point1,point2};
//                       tetra_planes.resize(4);
//                       tetra_planes[0]=vector<Size>{pointj,point1,point2};
//                       tetra_planes[1]=vector<Size>{pointi,point1,point2};
//                       tetra_planes[2]=vector<Size>{pointi,pointj,point2};
//                       tetra_planes[3]=vector<Size>{pointi,pointj,point1};
//                       for (Size k=0; k<4; k++)
//                       {
//                         double ref_orient=orientation(tetra_planes[k],curr_point);
//                         double comp_orient=orientation(tetra_planes[k],orient_points[k]);
//                         if (ref_orient*comp_orient>=0)//then the last point of the tetrahedron lies on the same side as the deleted point, so this plane is unusable
//                         {
//                           forbidden_planes.insert(tetra_planes[k]);
//                           // std::cout << "Forbidden plane found: " << tetra_planes[k][0] << "," << tetra_planes[k][1] << "," << tetra_planes[k][2] << std::endl;
//                         }
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//         //finally adding edge planes to plane_counter
//         for (vector<Size> edge_plane:edge_planes)
//         {
//           if (!is_forbidden_plane(edge_plane,forbidden_planes))
//           {
//             plane_counter.insert(std::pair<vector<Size>,Size>(edge_plane,0));//Initially, the edge planes are not adjacent to any tetraherda (ignoring ofcourse those forbidden planes)
//           }
//         }
//
//         //std::set<vector<Size>> neighbor_lines;//this set contains all valid possible lines; will be used to construct new tetrahedra
//         // invariant: lines..[0]<lines..[1]
//         //TODO: replace this entire thing by first calculating all planes and then try adding a single point to it (and check orientation)
//
//         // vector<vector<Size>> planes;
//         // //TODO calc planes
//         //
//         // for (vector<Size> plane: planes){
//         //   generate_initial_ears(neighbors_of_point, plane, neighbor_map,
//         //    ears_map, rev_ears_map, curr_point);
//         // }
//
//         for (Size i=0; i<neighbors_of_point.size(); i++)
//         {
//           for (Size j=0; j<i; j++)
//           {
//             if (neighbor_map[neighbors_of_point[i]].count(neighbors_of_point[j])==0)//if those two points are not yet neighbors
//             {
//               std::set<Size> temp_intersection;
//               std::set_intersection(neighbor_map[neighbors_of_point[i]].begin(),neighbor_map[neighbors_of_point[i]].end(),
//                             neighbor_map[neighbors_of_point[j]].begin(),neighbor_map[neighbors_of_point[j]].end(),
//                             std::inserter(temp_intersection,temp_intersection.begin()));
//               //then for every pair in the intersection of the neighbors of both points, check if they are neighbors
//               //a better implementation is probably possible
//
//               for (Size point1: temp_intersection)
//               {
//                 for (Size point2: temp_intersection)
//                 {
//                   if (point1<point2 && //just such that we do not have duplicates
//                   neighbor_map[point1].count(point2)!=0)///if point1 and 2 are neighbors
//                   {//TODO: this is not guaranteed to be a convex hull!!!
//                     //THEREFORE check if this proposed ear actually lies INSIDE the hull!!!!!
//                     vector<Size> new_triangle=vector<Size>{neighbors_of_point[i],neighbors_of_point[j],point1,point2};
//                     // bool has_forbidden_plane=false;
//                     // for (vector<Size> forbidden_plane:forbidden_planes)
//                     // {//if triangle has forbidden plane
//                     //   if (calculate_total_points(forbidden_plane,new_triangle)==4)
//                     //   {
//                     //     has_forbidden_plane=true;
//                     //     std::cout << "Forbidden plane detected" << std::endl;
//                     //     break;
//                     //   }
//                     // }
//                     if (!has_forbidden_plane(new_triangle,forbidden_planes))//if we are still allowed to add this triangle to the ears list
//                     {
//                       double orient_pl1=orientation(vector<Size>{new_triangle[1],new_triangle[2],new_triangle[3]},new_triangle[0]);
//                       double orient_pl1_ref=orientation(vector<Size>{new_triangle[1],new_triangle[2],new_triangle[3]},curr_point);
//
//                       double orient_pl2=orientation(vector<Size>{new_triangle[0],new_triangle[2],new_triangle[3]},new_triangle[1]);
//                       double orient_pl2_ref=orientation(vector<Size>{new_triangle[0],new_triangle[2],new_triangle[3]},curr_point);
//                       // std::cout << "Current triangle: "<< new_triangle[0] << ", " << new_triangle[1] << ", " << new_triangle[2] << ", " << new_triangle[3] << std::endl;
//                       // std::cout << "orientation of triangle 1: " << orient_pl1 << std::endl;
//                       // std::cout << "reference orientation of triangle 1: " << orient_pl1_ref << std::endl;
//                       // std::cout << "orientation of triangle 2: " << orient_pl2 << std::endl;
//                       // std::cout << "reference orientation of triangle 2: " << orient_pl2_ref << std::endl;
//                       if (((orient_pl1*orient_pl1_ref)>=0)&&((orient_pl2*orient_pl2_ref)>=0))
//                       {
//                         double power=calc_power(new_triangle,curr_point);
//                         if (!std::isnan(power))
//                         {//otherwise we would be proposing a coplanar tetrahedron, which would be ridiculous
//                           ears_map.insert(std::make_pair(new_triangle,power));
//                           rev_ears_map.insert(std::make_pair(power,new_triangle));
//                           // std::cout << "Creating ear: "<< new_triangle[0] << ", " << new_triangle[1] << ", " << new_triangle[2] << ", " << new_triangle[3] << std::endl;
//                           //invariant: the first two element of the vector should correspond to the neighbors we want to add
//                         }
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//         // for (Size point1: neighbors_of_point)
//         // {
//         //   for (Size point2: neighbor_map[point1])
//         //   {
//         //     if (point1<point2)
//         //     {
//         //       neighbor_lines.insert(vector<Size>{point1,point2});
//         //     }
//         //   }
//         // }
//         //or just print all lines out...
//         //std::cout << "Initial lines size: " << neighbor_lines.size();
//         //std::cout << "number of neigbors: " << neighbors_of_point.size() << std::endl;
//         //std::cout << "ears_map size: " << ears_map.size() << std::endl;
//         //now that we have all initial 'triangles', we can finally start to add them
//         //while(ears_map.size()>7)//in the last iteration, you should have an octohedron left: -> 8 possible ears
//         while(!ears_map.empty())
//         {
//           // std::cout << "Looping: current size of ears map: " << ears_map.size() << std::endl;
//           // std::cout << "Looping: current size of lines: " << neighbor_lines.size() << std::endl;
//           vector<Size> triangle=(*rev_ears_map.begin()).second;//triangle which we are currently adding to the mesh
//           //i we somehow end up with a line we no longer need to add, ignore it and repeat the loop
//           if ((point_surrounded_by_tetras(triangle[0],plane_counter,edge_planes))||(point_surrounded_by_tetras(triangle[1],plane_counter,edge_planes)))
//           {
//             // std::cout<<"Proposed ear no longer useful"<<std::endl;
//             // std::cout << "Deleted ear: " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << ", " << triangle[3] << std::endl;
//             delete_vector_from_both_maps(ears_map,rev_ears_map,triangle);
//             continue;
//           }
//           // std::cout << "Currently adding ear: (" << triangle[0] << "," << triangle[1] << "," << triangle[2] << "," << triangle[3] << ")" << std::endl;
//           vector<vector<Size>> triangles_to_work_with;//vector of triangles from which we will generate  new triangles
//           // triangles_to_work_with.push_back(triangle);
//           neighbors_lists[curr_coarsening_lvl].add_neighbor(triangle[0], triangle[1]);
//
//           //neighbors_lists[curr_coarsening_lvl].add_single_neighbor(triangle[1], triangle[0]); Changed the function, now also adds the point as neighbor of the neighbor
//           neighbor_map[triangle[0]].insert(triangle[1]);
//           neighbor_map[triangle[1]].insert(triangle[0]);
//           // forbidden_planes.insert(vector<Size>{triangle[0],triangle[2],triangle[3]});
//           // forbidden_planes.insert(vector<Size>{triangle[1],triangle[2],triangle[3]});
//           //debug stuff
//           added_lines[i].push_back(vector<Size>{triangle[0],triangle[1]});
//           vector<vector<Size>> debug_tetra_to_add;
//           // debug_tetra_to_add.push_back(vector<Size>{triangle});
//
//           //invariant: the first two element of the vector should correspond to the neighbors we want to add
//           // std::cout << "Deleted ear: " << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << ", " << triangle[3] << std::endl;
//           delete_vector_from_both_maps(ears_map,rev_ears_map,triangle);
//           //Delete corresponding line, which may no longer be used
//           // vector<Size> line_to_delete{triangle[2],triangle[3]};
//           // std::sort(line_to_delete.begin(),line_to_delete.end());
//           //deleting the old line from the lines and the neighbors_map
//           //neighbor_lines.erase(line_to_delete);
//           // neighbor_map[line_to_delete[0]].erase(line_to_delete[1]);
//           // neighbor_map[line_to_delete[1]].erase(line_to_delete[0]);
//           //check which possible ears need to be deleted and which need to be added to the list from which to generate new ears
//
//           std::set<vector<Size>> temp_collection_tetras;//temporarily keeps the tetras to finally add them to all_added_tetras
//           //calculates exactly which tetrahedra are formed with the added line and adds them to triangles_to_work_with
//           std::set<Size> temp_intersection;
//           std::set_intersection(neighbor_map[triangle[0]].begin(),neighbor_map[triangle[0]].end(),
//                         neighbor_map[triangle[1]].begin(),neighbor_map[triangle[1]].end(),
//                         std::inserter(temp_intersection,temp_intersection.begin()));
//           for (Size temp_point:temp_intersection)
//           {
//             for (Size temp_point2:temp_intersection)
//             {
//               if ((temp_point<temp_point2)&&(neighbor_map[temp_point].count(temp_point2)!=0))
//               {
//                 vector<Size> curr_tetra{triangle[0],triangle[1],temp_point,temp_point2};
//                 // if (!has_forbidden_plane(curr_tetra,forbidden_planes))//as usual, the forbidden planes are still not allowed
//                   // {
//                   bool may_add=true;//now we have an edge case
//                   for (vector<Size> old_tetra: all_added_tetras)
//                   {
//                     //if shares plane with tetra already added,
//                     if (triangle_shares_plane_with(old_tetra, curr_tetra))
//                     {
//                       //find point of old tetra not in curr tetra and vice versa; also get shared triangle
//                       Size old_point=0;
//                       vector<Size> shared_triangle;
//                       for (Size temp_point: old_tetra)
//                       {if(!vector_contains_element(curr_tetra, temp_point)){old_point=temp_point;}else{shared_triangle.push_back(temp_point);}}
//                       Size new_point=0;
//                       for (Size temp_point: curr_tetra)
//                       {if(!vector_contains_element(old_tetra, temp_point)){new_point=temp_point;}}
//                       // std::cout << "shared triangle: " << shared_triangle[0] << ", " << shared_triangle[1] << ", " << shared_triangle[2] << std::endl;
//                       // std::cout << "old point: " << old_point << std::endl;
//                       // std::cout << "new point: " << new_point << std::endl;
//                       //if they both lie on the same side, ignore the new triangle
//                       if (orientation(shared_triangle,old_point)*orientation(shared_triangle,new_point)>=0)
//                       {
//                         // std::cout << "tetra not added" << std::endl;
//                         may_add=false;
//                         break;
//                       }
//                     }
//                   }
//                   if (may_add)
//                   {
//                     triangles_to_work_with.push_back(curr_tetra);
//                     debug_tetra_to_add.push_back(curr_tetra);
//                     add_planes_of_tetra_to_counter(curr_tetra,plane_counter);
//                     temp_collection_tetras.insert(curr_tetra);
//                     // all_added_tetras.insert(curr_tetra);
//                   }
//                 // }
//               }
//             }
//           }
//           for (vector<Size> curr_tetra:temp_collection_tetras)//and finally adding the new batch of tetras
//           {all_added_tetras.insert(curr_tetra);}
//
//           //TODO remove this section and just calculate which tetrahedra are formed!
//           // for (auto it = ears_map.cbegin(); it != ears_map.cend();)
//           // {
//           //   vector<Size> curr_triangle=(*it).first;//candidate triangle to check whether it needs to be changed/deleted
//           //   //if this candidate triangle would try to add the same neighbors, we just delete it (although i'm not sure whether it can happen before the mesh is almost complete...)
//           //   if ((curr_triangle[0]==triangle[0]&&curr_triangle[1]==triangle[1])||
//           //       (curr_triangle[1]==triangle[0]&&curr_triangle[0]==triangle[1]))
//           //   {
//           //     // bool already_in=false;
//           //     // //checks whether triangles_to_work_with already has a permutation of  curr_triangle
//           //     // for (vector<Size> temp_triangle:triangles_to_work_with)
//           //     // {
//           //     //   if (triangles_are_equal(temp_triangle,curr_triangle))
//           //     //   {
//           //     //     already_in=true;
//           //     //     break;
//           //     //   }
//           //     // }
//           //     // //if we did not already put some permutation of the triangle in triangles_to_work_with, we put it in and update the forbidden planes
//           //     // if (!already_in)
//           //     // {
//           //     //   triangles_to_work_with.push_back(curr_triangle);
//           //     //   // forbidden_planes.insert(vector<Size>{curr_triangle[0],curr_triangle[2],curr_triangle[3]});
//           //     //   // forbidden_planes.insert(vector<Size>{curr_triangle[1],curr_triangle[2],curr_triangle[3]});
//           //     //   //debug stuff
//           //     //   debug_tetra_to_add.push_back(curr_triangle);
//           //     // }
//           //
//           //     std::cout << "Deleted ear: " << curr_triangle[0] << ", " << curr_triangle[1] << ", " << curr_triangle[2] << ", " << curr_triangle[3] << std::endl;
//           //     it=delete_vector_from_both_maps(ears_map,rev_ears_map,curr_triangle);
//           //     //Delete corresponding line, which may no longer be used
//           //     // vector<Size> line_to_delete{curr_triangle[2],curr_triangle[3]};
//           //     // std::sort(line_to_delete.begin(),line_to_delete.end());
//           //     // neighbor_lines.erase(line_to_delete);
//           //     // neighbor_map[line_to_delete[0]].erase(line_to_delete[1]);
//           //     // neighbor_map[line_to_delete[1]].erase(line_to_delete[0]);
//           //   }
//           //   else//we do nothing with the current entry for now
//           //   {
//           //     ++it;
//           //   }
//           // }
//           //debug stuff
//           added_tetras[i].push_back(debug_tetra_to_add);
//
//           //Now we are deleting all ears that share a plane with a triangle in triangles_to_work_with
//           //TODO: maybe check if this here can be made faster
//           for (auto it = ears_map.cbegin(); it != ears_map.cend();)
//           {
//             vector<Size> curr_triangle=(*it).first;
//             bool deleted=false;
//             for (vector<Size> triangle: triangles_to_work_with)
//             {
//               if (triangle_shares_plane_with(triangle, curr_triangle))
//               {
//                 // std::cout << "Deleted ear: " << curr_triangle[0] << ", " << curr_triangle[1] << ", " << curr_triangle[2] << ", " << curr_triangle[3] << std::endl;
//                 it = delete_vector_from_both_maps(ears_map, rev_ears_map, curr_triangle);
//                 deleted=true;
//                 break;
//               }
//             }
//             if(!deleted)
//             {
//               ++it;
//             }
//           }
//
//           //std::cout << "nb triangles to work with: " << triangles_to_work_with.size() << std::endl;
//           std::pair<vector<Size>,vector<vector<Size>>> relevant_planes_pair=return_all_relevant_planes(triangles_to_work_with,forbidden_planes);
//           vector<Size> orient_points=relevant_planes_pair.first;
//           vector<vector<Size>> relevant_planes=relevant_planes_pair.second;
//           //std::cout << "number of planes: " << relevant_planes.size() << std::endl;
//           //TODO FIXME: check if the same ear can be generated twice and check if it would actually matter...
//           // std::cout << "relevant planes size: "<<relevant_planes.size() << std::endl;
//
//           // //if we have two planes that already form a tetrahedron, we can stop adding new ears
//           // if (relevant_planes.size()==2)
//           // {
//           //   // first two elements of the planes correspond to the line added and are therefore shared
//           //   vector<Size> first_plane=relevant_planes[0];
//           //   vector<Size> second_plane=relevant_planes[1];
//           //   Size point_to_check=second_plane[2];
//           //   Size count_neighbors=0;
//           //   for (Size temp_point: first_plane)
//           //   {
//           //     if (neighbor_map[temp_point].count(point_to_check)!=0)//if temp_point is neighbor of point_of_plane
//           //     {count_neighbors++;}
//           //   }
//           //   double prod_of_orient=orientation(first_plane,point_to_check)*orientation(first_plane,orient_points[0]);//should be negative if temp_point lies on the opposite side of the plane as the orientation point
//           //   if ((count_neighbors==3)&&(prod_of_orient<0))
//           //   {
//           //     std::cout << "Two proposed planes already form a new tetrahedron. Therefore we will not use them." << std::endl;
//           //     relevant_planes.clear();
//           //     orient_points.clear();
//           //
//           //   }
//           // }
//           // else if (relevant_planes.size()==0)//i think that our iterations must stop then... i consider this a hack for now
//           // {
//           //   ears_map.clear();
//           //   rev_ears_map.clear();
//           // }
//
//
//           //if for some reason (namely delaunay grids being weird), we still havent added a plane squished between two completed tetrahedra
//           // we add it now TODO check if this works...
//           // std::set<Size> temp_intersection;
//           // std::set_intersection(neighbor_map[triangle[0]].begin(),neighbor_map[triangle[0]].end(),
//           //               neighbor_map[triangle[1]].begin(),neighbor_map[triangle[1]].end(),
//           //               std::inserter(temp_intersection,temp_intersection.begin()));
//           // for (Size temp_point:temp_intersection)
//           // {
//           //   std::set<Size> temp_intersection2;//intersection of neighbors of temp_point with
//           //   std::set_intersection(temp_intersection.begin(),temp_intersection.end(),
//           //                 neighbor_map[temp_point].begin(),neighbor_map[temp_point].end(),
//           //                 std::inserter(temp_intersection2,temp_intersection2.begin()));
//           //   if (temp_intersection2.size()>=2){//note that i do not know how it could be larger than 2, but i'm taking no chances
//           //     forbidden_planes.insert(vector<Size>{triangle[0],triangle[1],temp_point});
//           //   }
//           // }
//           Size count=0;
//           for (vector<Size> plane: relevant_planes)
//           {//TODO: figure out whether new ears must be generated at all, because no new ears ever seem to be inserted...
//           generate_new_ears(neighbors_of_point, plane, neighbor_map,
//                 ears_map, rev_ears_map, curr_point, orient_points[count], plane_counter, edge_planes, forbidden_planes);
//           count++;
//           }
//
//
//
//           //and finally delete those lines from neighbor_map
//           //delete_all_useless_lines(triangles_to_work_with,neighbor_map);
//           // generate_new_ears(neighbor_lines, new_line, neighbor_map,
//           //       ears_map, rev_ears_map, curr_point);
//
//           //and finally add the new line to neighbor_lines
//           //neighbor_lines.insert(new_line);
//         }
//
//         //DEBUG STUFF
//         //Check if there is a full triangulisation
//         if (debug_mode){
//           std::set<vector<Size>> wrong_planes;
//           for (Size point: neighbors_of_point)
//           {
//             if (!point_surrounded_by_tetras(point, plane_counter, edge_planes))
//             {
//               // std::cout<<"Iteration: "<<i<<" point: "<<point<<std::endl;
//               for (auto pair:plane_counter)
//               {
//                 vector<Size> curr_plane=pair.first;
//                 Size count=pair.second;
//                 bool is_edgepl=is_edge_plane(curr_plane,edge_planes);
//                 if ((count<2&&!is_edgepl)||(count<1&&is_edgepl))
//                 {
//                 wrong_planes.insert(curr_plane);
//                 // std::cout<<curr_plane[0]<<", "<<curr_plane[1]<<", "<<curr_plane[2]<<" counted nb times: "<<count<<std::endl;
//                 }
//               }
//             }
//           }
//           if (wrong_planes.size()>0)
//           {
//             bool trulywrong=false;//our heuristic can be wrong at times, so we sometimes need to recount things
//             for (vector<Size> curr_plane:wrong_planes)
//             {
//               Size recount=0;
//               // std::cout<<curr_plane[0]<<", "<<curr_plane[1]<<", "<<curr_plane[2]<<std::endl;
//               for (Size temp_point:neighbors_of_point)
//               { if (neighbor_map[temp_point].count(curr_plane[0])!=0&&neighbor_map[temp_point].count(curr_plane[1])!=0&&neighbor_map[temp_point].count(curr_plane[2])!=0)
//                 {recount++;}
//               }
//               if ((is_edge_plane(curr_plane,edge_planes)&&(recount<1))||(!is_edge_plane(curr_plane,edge_planes)&&(recount<2)))
//               {trulywrong=true;break;}
//             }
//             if (trulywrong){
//               throw "Not a full triangulation!!!";//will probably translate to unknown exception in python
//             }
//           }
//         }
//
//
//         reduced_neighbors_after[i]=std::map<Size, std::set<Size>>(neighbor_map);
//
//
//         //recalc density diff of neighbors (due to new mesh (and ofc deletion of point))
//         for (Size neighbor :neighbors_of_point)
//         {//again, we will never try to delete points on the boundary
//           if (geometry.not_on_boundary(neighbor))
//           {
//           delete_int_from_both_maps(density_diff_map,rev_density_diff_map,neighbor);
//
//           double calc_diff=calc_diff_abundance_with_neighbours(neighbor,curr_coarsening_lvl);
//           density_diff_map.insert(std::pair<Size,double>(neighbor,calc_diff));
//           rev_density_diff_map.insert(std::pair<double,Size>(calc_diff,neighbor));
//           }
//         }
//
//       //finally, delete the neighbors of this point
//       neighbors_lists[curr_coarsening_lvl].delete_all_neighbors(curr_point);
//       //and remove this entry from density maps
//       delete_int_from_both_maps(density_diff_map, rev_density_diff_map, curr_point);
//     }
//     current_nb_points-=nb_points_to_remove;//decrement the number of points remaining
//     nb_points_at_lvl.push_back(current_nb_points);
//   }
//   else
//   {//TODO check whether deep copy!!!!!!
//     curr_coarsening_lvl++;
//     geometry.points.curr_neighbors=neighbors_lists[curr_coarsening_lvl];//should be shallow copy
//   }
// }
//
//
// /// Resets the grid
// inline void Model :: reset_grid()
// {
//   curr_coarsening_lvl=0;
//   max_reached_coarsening_lvl=0;
//   density_diff_map=std::multimap<Size,double>();
//   rev_density_diff_map=std::multimap<double,Size>();
//   geometry.points.curr_neighbors=Neighbors(neighbors_lists[0]);
//   neighbors_lists=vector<Neighbors>();
//   current_nb_points=0;
// }
//
// /// Calculates all tetrahedra that contain a given point and given the coarsening level
// ///     @Returns: all tetrahedra (each tetrahedron has their points sorted)
// inline std::set<vector<Size>> Model :: calc_all_tetra_with_point(Size point, Size coars_lvl)
// {
//   std::set<vector<Size>> to_return;
//   vector<Size> curr_neighbors=vector<Size>(neighbors_lists[coars_lvl].get_neighbors(point));
//   std::sort(curr_neighbors.begin(),curr_neighbors.end());//copying and sorting our vectors, because we want to intersect them
//   vector<Size>::iterator it = curr_neighbors.begin();
//   while(it!=curr_neighbors.end())
//   {
//     Size neighbor=*it;
//     vector<Size> temp_intersection;
//     vector<Size> neighbors_of_neighbor=vector<Size>(neighbors_lists[coars_lvl].get_neighbors(neighbor));
//     std::sort(neighbors_of_neighbor.begin(),neighbors_of_neighbor.end());
//     std::set_intersection(curr_neighbors.begin(),curr_neighbors.end(),
//                           neighbors_of_neighbor.begin(),neighbors_of_neighbor.end(),
//                           back_inserter(temp_intersection));
//     if (temp_intersection.size()>=2)//then we check if we can truly form a tetrahedon (if we would not delete elements from curr_neighbors, this should always be true)
//     {
//       Size nb_elements=temp_intersection.size();
//       for (Size i=0; i<nb_elements; i++)
//       {
//         for (Size j=0; j<i; j++)
//         {//also check if these two points i and j are actually neighbors of eachother
//           if (vector_contains_element(neighbors_lists[coars_lvl].get_neighbors(i),j))
//           {
//           vector<Size> to_add({point, neighbor, temp_intersection[i], temp_intersection[j]});
//           std::sort(to_add.begin(),to_add.end());
//           to_return.insert(to_add);
//           }
//         }
//       }
//     }
//     // We are deleting the first element each time in order to not propose duplicates
//     //For the explanation, assume neighbors a<b<c: without this strange optimisation,
//     // we would try to add (before sorting them): (point,a,b,c), (point,b,a,c) and (point,c,a,b)
//     //Now these two last options are no longer there.
//     it=curr_neighbors.erase(it);
//   }
//   return to_return;
// }
//
//
// /// Calculates the barycentric coordinates of a point given the tetrahedron
// ///   @Parameter [in] triangle: the tetrahedron for which we calculate these coordinates
// ///   @Parameter [in] point: the point for which we calculate these coordinates
// ///   Returns: the barycentric coordinates of the point with respect to the tetrahedron
// inline Eigen::Vector<double,4> Model :: calc_barycentric_coords(const vector<Size> &triangle, Size point){
//   Vector3D pos1=geometry.points.position[triangle[0]];
//   Vector3D pos2=geometry.points.position[triangle[1]];
//   Vector3D pos3=geometry.points.position[triangle[2]];
//   Vector3D pos4=geometry.points.position[triangle[3]];
//   Vector3D posp=geometry.points.position[point];//position of point
//   Vector3D temp_diff=posp-pos4;
//
//   Eigen::Vector3d diff;
//   diff << temp_diff.x(), temp_diff.y(), temp_diff.z();
//
//   //see wikipedia for explanation (for 2D case, but easily generalizable to 3D)
//   Eigen::Matrix<double,3,3> mat;
//   mat << pos1.x()-pos4.x(),pos2.x()-pos4.x(),pos3.x()-pos4.x(),
//          pos1.y()-pos4.y(),pos2.y()-pos4.y(),pos3.y()-pos4.y(),
//          pos1.z()-pos4.z(),pos2.z()-pos4.z(),pos3.z()-pos4.z();
//
//   Eigen::Vector<double,4> result;
//   result << mat.colPivHouseholderQr().solve(diff);
//   result[3]=1-result[0]-result[1]-result[2];
//   return result;
// }

// // The radial basis function we apply
//   inline double rbf(double radius) // the functor we want to apply, note
//   {//TODO: use sensible function instead; also let it depend on some predertemined radius; or just divide all elements by that number
//     // return std::exp(-x);
//     return std::max(1-radius,0.0);
//   }
//
// /// Interpolates the vector with length nb_points_in_coarser level
// ///   @Parameter [in] coarser_lvl: the coarsening level of the coarser grid
// ///   @Parameter [in] finer_lvl: the coarsening level of the finer grid; should be smaller than coarser_lvl
// ///   @Parameter [in/out] to_interpolate: the vector with values at the points of the coarser grid (and some irrelevant values inbetween), has length equal to the total number of points
// ///   @Returns: A vector containing the interpolated values (has length equal to parameters.npoints())
// inline void Model::interpolate_vector(Size coarser_lvl, Size finer_lvl, vector<double> &to_interpolate)
// {//TODO: only use large vector and some masks
//   //FIXME: in the case that we somehow have nothing in diff_points, just quit!!!
//   Size nb_points=parameters.npoints();
//   std::vector<Size> all_points(nb_points);
//   std::iota(all_points.begin(), all_points.end(), 0); // all_points will become: [0..nb_points]
//   vector<bool> coarse_mask=geometry.points.multiscale.get_mask(coarser_lvl);
//   vector<bool> finer_mask=geometry.points.multiscale.get_mask(finer_lvl);
//
//   //TODO calculate the difference
//   vector<Size> diff_points;
//   vector<double> to_interpolate_values;//values of the coarse grid we need to interpolate (after applying the mask)
//   vector<Size> coarse_points;
//   for (Size point: all_points)
//   {
//     if(finer_mask[point])
//     {
//       if (coarse_mask[point])
//       {
//         coarse_points.push_back(point);
//         to_interpolate_values.push_back(to_interpolate[point]);
//       }else{diff_points.push_back(point);}
//     }
//   }
//
//   Size nb_coarse_points=coarse_points.size();
//   Size nb_diff_points=diff_points.size();
//
//   //if we have truly nothing to do, just do nothing
//   if (nb_diff_points==0)
//   {return;}
//
//
//   std::vector<Eigen::Triplet<double>> tripletList;
//   tripletList.reserve(50*nb_coarse_points);//TODO: get better estimate
//
//   //the matrix A for A*weights=to_interpolate
//   //Eigen::MatrixXd coarse_rbf_mat(nb_coarse_points, nb_coarse_points);
//   //possibly better implementation possible, this is far too slow!!
//   double maxdist=0;
//   for (Size row=0; row<nb_coarse_points; row++)
//   {
//     std::cout<<row<<std::endl;
//     std::set<Size> neighbors_of_row=geometry.points.multiscale.get_neighbors(row);
//     for (Size col=0; col<nb_coarse_points; col++)
//     {
//       double dist=(geometry.points.position[coarse_points[row]]-geometry.points.position[coarse_points[col]]).squaredNorm();
//       dist=sqrt(dist);
//         // also trying to calculate the maximum relevant distance between points (needed for rescale for radial basis function)
//       if (neighbors_of_row.find(col)!=neighbors_of_row.end()&&dist>maxdist)
//       {
//         maxdist=dist;
//       }
//     }
//   }
//
//   for (Size row=0; row<nb_coarse_points; row++)
//   {
//     std::cout<<row<<std::endl;
//     std::set<Size> neighbors_of_row=geometry.points.multiscale.get_neighbors(row);
//     for (Size col=0; col<nb_coarse_points; col++)
//     {
//       double dist=(geometry.points.position[coarse_points[row]]-geometry.points.position[coarse_points[col]]).squaredNorm();
//       //coarse_rbf_mat(row,col)=dist;
//       if (dist<maxdist)
//       {
//         dist=sqrt(dist);
//         tripletList.push_back(Eigen::Triplet<double>(row,col,dist));
//         // also trying to calculate the maximum relevant distance between points (needed for rescale for radial basis function)
//         if (neighbors_of_row.find(col)!=neighbors_of_row.end()&&dist>maxdist)
//         {
//           maxdist=dist;
//         }
//       }
//     }
//   }
//   Eigen::SparseMatrix<double> coarse_rbf_mat(nb_coarse_points,nb_coarse_points);
//   coarse_rbf_mat.setFromTriplets(tripletList.begin(), tripletList.end());
//
//   coarse_rbf_mat=coarse_rbf_mat/maxdist;
//   coarse_rbf_mat=coarse_rbf_mat.unaryExpr(std::ptr_fun(rbf));
//
//   //now solving A*weights=to_interpolate_values
//   //Note: A is symmetrical, positive definite
//   //NOTE: if this were to be too slow, switch to LLT instead (but that would be less accurate)
//   Eigen::VectorXd right_hand_side(nb_coarse_points);
//   for (Size coarse_idx=0; coarse_idx<nb_coarse_points; coarse_idx++)
//   {
//     right_hand_side(coarse_idx)=to_interpolate_values[coarse_points[coarse_idx]];
//   }
//
//   std::cout<<"probably expensive solve operation"<<std::endl;
//   Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> decomp(coarse_rbf_mat);  // performs a Cholesky factorization of A
//   Eigen::VectorXd weights = decomp.solve(right_hand_side);
//   //Eigen::VectorXd weights=coarse_rbf_mat.ldlt().solve(right_hand_side);
//
//   std::cout<<"end expensive solve operation"<<std::endl;
//
//   //now construct the matrix to interpolate the other values
//   std::vector<Eigen::Triplet<double>> tripletList_diff;
//   tripletList_diff.reserve(200*nb_diff_points);//TODO: get better estimate
//   //the matrix A for x=A*weights
//   // Eigen::MatrixXd diff_rbf_mat(nb_diff_points, nb_coarse_points);
//   //possibly better implementation possible, but this suffices for now
//   //can probably be improved...
//   for (Size row=0; row<nb_diff_points; row++)
//   {
//     for (Size col=0; col<nb_coarse_points; col++)
//     {
//       double dist=(geometry.points.position[diff_points[row]]-geometry.points.position[coarse_points[col]]).squaredNorm();
//       //diff_rbf_mat(row,col)=dist;
//       if (dist<maxdist)
//       {tripletList_diff.push_back(Eigen::Triplet<double>(row,col,sqrt(dist)));}
//     }
//   }
//   Eigen::SparseMatrix<double> diff_rbf_mat(nb_diff_points,nb_coarse_points);
//   diff_rbf_mat.setFromTriplets(tripletList_diff.begin(), tripletList_diff.end());
//
//   //normalize the distance again and applying the radial basis function
//   diff_rbf_mat=diff_rbf_mat/maxdist;
//   diff_rbf_mat=diff_rbf_mat.unaryExpr(std::ptr_fun(rbf));
//
//   VectorXd interpolated_values=diff_rbf_mat*weights;
//
//   for (Size diff_idx=0; diff_idx<nb_diff_points; diff_idx++)
//   {
//     to_interpolate[diff_points[diff_idx]]=interpolated_values[diff_idx];
//   }
// }
//

/// Initializes the multigrid
///   @Parameter[in]: min_nb_points: stop coarsening when number of points of finest layers is equal or less than this
///   @Parameter[in]: max_coars_lvl: The maximum coarsening level we allow
///   @Parameter[in]: tol: the tolerance level for which we still consider points similar enough
inline int Model::setup_multigrid(Size min_nb_points, Size max_coars_lvl, double tol)
{
  auto fun_to_del=points_are_similar(tol);//function that says if two points are similar enough
  geometry.points.multiscale.set_not_on_boundary_fun([&](Size p){return geometry.not_on_boundary(p);});//function that says whether a point lies on the boundary
  geometry.points.multiscale.set_comparison_fun(fun_to_del);
  //first, we coarsen the grid until we either have too few points left or have too many coarsening levels
  while((geometry.points.multiscale.get_max_coars_lvl()<max_coars_lvl)&&
  (geometry.points.multiscale.get_total_points(geometry.points.multiscale.get_max_coars_lvl())>min_nb_points))
  {
    geometry.points.multiscale.coarsen();
  }
}

/// Computes second order feautrier using multigrid
inline int Model::compute_feautrier_order_2_multigrid()
{
  //set curr coarsening level to max
  geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_max_coars_lvl());
  //solve and interpolate J for all coarser grids
  while(geometry.points.multiscale.get_curr_coars_lvl()>0)
  {
    compute_radiation_field_feautrier_order_2();
    //for all frequencies, interpolate J
    interpolate_matrix_local(geometry.points.multiscale.get_curr_coars_lvl(),radiation.J);
    geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_curr_coars_lvl()-1);
  }
  //finally, solve for the final grid
  compute_radiation_field_feautrier_order_2();
  return (0);
}


template <typename T>
inline T rbf_local(T radius)
{
  return std::exp(-std::pow(radius,2));
}

/// Interpolates the vector with length parameters.npoints() from the coarser level to the first finer level
///   @Parameter [in] coarser_lvl: the coarsening level of the coarser grid
///   @Parameter [in/out] to_interpolate: the vector with values at the points of the coarser grid (and some irrelevant values inbetween), has length equal to the total number of points
template <typename T>//T should be double, long or long double
inline void Model::interpolate_vector_local(Size coarser_lvl, vector<T> &to_interpolate)
{//TODO: only use large vector and some masks
  //we cannot interpolate from lvl 0 to some level -..., so just return
  if (coarser_lvl==0)
  {return;}
  Size nb_points=parameters.npoints();
  std::vector<Size> all_points(nb_points);
  std::iota(all_points.begin(), all_points.end(), 0); // all_points will become: [0..nb_points-1]
  vector<bool> coarse_mask=geometry.points.multiscale.get_mask(coarser_lvl);
  vector<bool> finer_mask=geometry.points.multiscale.get_mask(coarser_lvl-1);

  //TODO calculate the difference
  vector<Size> diff_points;
  // vector<double> to_interpolate_values;//values of the coarse grid we need to interpolate (after applying the mask)
  vector<Size> coarse_points;
  for (Size point: all_points)
  {
    if(finer_mask[point])
    {
      if (coarse_mask[point])
      {
        coarse_points.push_back(point);
        // to_interpolate_values.push_back(to_interpolate[point]);
      }else{diff_points.push_back(point);}
    }
  }

  Size nb_coarse_points=coarse_points.size();
  Size nb_diff_points=diff_points.size();

  //if we have truly nothing to do, just do nothing
  if (nb_diff_points==0)
  {return;}


  //for every points in diff_points, try to find interpolating value using rbf function
  for (Size diff_point: diff_points)
  {
    std::set<Size> curr_neighbors=geometry.points.multiscale.get_neighbors(diff_point,coarser_lvl-1);
    //get neighbors in coarser grid
    vector<Size> neighbors_coarser_grid;
    neighbors_coarser_grid.reserve(curr_neighbors.size());//i am overallocating a bit, but this variable is temporary anyway...
    for (Size neighbor: curr_neighbors)
    {
      if (coarse_mask[neighbor])
      {
        neighbors_coarser_grid.push_back(neighbor);
      }
    }
    std::cout << "number of neighbors: " << neighbors_coarser_grid.size() << std::endl;
    // //now we could also use the extremely fast method of interpolating: taking the average of the neighbors // i am not sure what actually gives the best results: rbf or just averaging
    // double interpolated_value=accumulate(neighbors_coarser_grid.begin(), neighbors_coarser_grid.end(), 0.0)/neighbors_coarser_grid.size();
    // to_interpolate[diff_point]=interpolated_value;


    // Note: using the current multigrid creation method, the number of neighbors in the coarse grid is always at least 1
    if (neighbors_coarser_grid.size()==1)//in the rare case that we only have a single neighbor in the coarse grid, just set the value to its value
    {
      to_interpolate[diff_point]=to_interpolate[neighbors_coarser_grid[0]];
    }
    else
    {//use rbf estimate
      Size nb_neighbors_coarser_grid=neighbors_coarser_grid.size();
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rbf_mat(nb_neighbors_coarser_grid, nb_neighbors_coarser_grid);
      Eigen::Vector<T,Eigen::Dynamic> right_hand_side(nb_neighbors_coarser_grid);
      Eigen::Matrix<T,1,Eigen::Dynamic> distance_with_neighbors(1,nb_neighbors_coarser_grid);
      for (Size idx=0; idx<nb_neighbors_coarser_grid; idx++)
      {
        right_hand_side(idx)=to_interpolate[idx];
        distance_with_neighbors(idx)=std::sqrt((geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[diff_point]).squaredNorm());
        rbf_mat(idx,idx)=0;//distance between point and itself is zero
        for (Size idx2=0; idx2<idx; idx2++)
        {
          //calculate radius
          T radius=std::sqrt((geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[neighbors_coarser_grid[idx2]]).squaredNorm());
          rbf_mat(idx,idx2)=radius;
          rbf_mat(idx2,idx)=radius;
        }
      }
      T maxdist=rbf_mat.maxCoeff()*5.0;//arbitrary number 5 to make the max_dist larger
      rbf_mat=rbf_mat/maxdist;
      rbf_mat=rbf_mat.unaryExpr(std::ptr_fun(rbf_local<T>));
      distance_with_neighbors=distance_with_neighbors/maxdist;
      distance_with_neighbors=distance_with_neighbors.unaryExpr(std::ptr_fun(rbf_local<T>));
      Eigen::Vector<T,Eigen::Dynamic> weights=rbf_mat.llt().solve(right_hand_side);
      T interpolated_value=(distance_with_neighbors*weights)(0,0);

      to_interpolate[diff_point]=interpolated_value;
    }

  }


}


/// Interpolates the paracabs matrix from the coarser level to the first finer level
///   @Parameter [in] coarser_lvl: the coarsening level of the coarser grid
///   @Parameter [in/out] to_interpolate: the paracabs matrix with values at the points of the coarser grid (and some irrelevant values inbetween), has length equal to the total number of points
template <typename T>//T should be double, long or long double
inline void Model::interpolate_matrix_local(Size coarser_lvl, Matrix<T> &to_interpolate)
{//TODO: only use large vector and some masks
  //we cannot interpolate from lvl 0 to some level -..., so just return
  if (coarser_lvl==0)
  {return;}
  Size nb_points=parameters.npoints();
  std::vector<Size> all_points(nb_points);
  std::iota(all_points.begin(), all_points.end(), 0); // all_points will become: [0..nb_points-1]
  vector<bool> coarse_mask=geometry.points.multiscale.get_mask(coarser_lvl);
  vector<bool> finer_mask=geometry.points.multiscale.get_mask(coarser_lvl-1);

  vector<Size> diff_points;
  // vector<double> to_interpolate_values;//values of the coarse grid we need to interpolate (after applying the mask)
  vector<Size> coarse_points;
  for (Size point: all_points)
  {
    if(finer_mask[point])
    {
      if (coarse_mask[point])
      {
        coarse_points.push_back(point);
        // to_interpolate_values.push_back(to_interpolate[point]);
      }else{diff_points.push_back(point);}
    }
  }

  Size nb_coarse_points=coarse_points.size();
  Size nb_diff_points=diff_points.size();

  //if we have truly nothing to do, just do nothing
  if (nb_diff_points==0)
  {return;}


  //for every points in diff_points, try to find interpolating value using rbf function
  for (Size diff_point: diff_points)
  {
    std::set<Size> curr_neighbors=geometry.points.multiscale.get_neighbors(diff_point,coarser_lvl-1);
    //get neighbors in coarser grid
    vector<Size> neighbors_coarser_grid;
    neighbors_coarser_grid.reserve(curr_neighbors.size());//i am overallocating a bit, but this variable is temporary anyway...
    for (Size neighbor: curr_neighbors)
    {
      if (coarse_mask[neighbor])
      {
        neighbors_coarser_grid.push_back(neighbor);
      }
    }
    // std::cout << "number of neighbors: " << neighbors_coarser_grid.size() << std::endl;
    // //now we could also use the extremely fast method of interpolating: taking the average of the neighbors // i am not sure what actually gives the best results: rbf or just averaging
    // double interpolated_value=accumulate(neighbors_coarser_grid.begin(), neighbors_coarser_grid.end(), 0.0)/neighbors_coarser_grid.size();
    // to_interpolate[diff_point]=interpolated_value;
    // TODO:maybe reorder these things for speed (we do not need to recalculate the radius every time)
    // Note: using the current multigrid creation method, the number of neighbors in the coarse grid is always at least 1
    if (neighbors_coarser_grid.size()==1)//in the rare case that we only have a single neighbor in the coarse grid, just set the value to its value
    {
      for (Size freqidx=0; parameters.nfreqs(); freqidx++)
      {
      to_interpolate[diff_point,freqidx]=to_interpolate[neighbors_coarser_grid[0],freqidx];
      }
    }
    else
    {//use rbf estimate
      Size nb_neighbors_coarser_grid=neighbors_coarser_grid.size();
      Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rbf_mat(nb_neighbors_coarser_grid, nb_neighbors_coarser_grid);
      Eigen::Matrix<T,1,Eigen::Dynamic> distance_with_neighbors(1,nb_neighbors_coarser_grid);
      for (Size idx=0; idx<nb_neighbors_coarser_grid; idx++)
      {
        distance_with_neighbors(idx)=std::sqrt((geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[diff_point]).squaredNorm());
        rbf_mat(idx,idx)=0;//distance between point and itself is zero
        for (Size idx2=0; idx2<idx; idx2++)
        {
          //calculate radius
          T radius=std::sqrt((geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[neighbors_coarser_grid[idx2]]).squaredNorm());
          rbf_mat(idx,idx2)=radius;
          rbf_mat(idx2,idx)=radius;
        }
      }
      T maxdist=rbf_mat.maxCoeff()*5.0;//arbitrary number 5 to make the max_dist larger
      rbf_mat=rbf_mat/maxdist;
      rbf_mat=rbf_mat.unaryExpr(std::ptr_fun(rbf_local<T>));
      distance_with_neighbors=distance_with_neighbors/maxdist;
      distance_with_neighbors=distance_with_neighbors.unaryExpr(std::ptr_fun(rbf_local<T>));
      Eigen::LLT<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>> lltdec(rbf_mat);

      //now that we have our llt decomposition, calculate the interpolated value
      for (Size freqidx=0; parameters.nfreqs(); freqidx++)
      {
        Eigen::Vector<T,Eigen::Dynamic> right_hand_side(nb_neighbors_coarser_grid);
        for (Size idx=0; idx<nb_neighbors_coarser_grid; idx++)
        {
          right_hand_side(idx)=to_interpolate[idx,freqidx];
        }
        Eigen::Vector<T,Eigen::Dynamic> weights=lltdec.solve(right_hand_side);
        T interpolated_value=(distance_with_neighbors*weights)(0,0);
        to_interpolate[diff_point,freqidx]=interpolated_value;
      }
    }

  }
}





  //for every element, calculate distance with other elements

  //Size count=0;
  // for (Size point=0; point<nb_points; point++)//counts number of points in finer grid
  // {if (mask_list[finer_lvl][point]){count++;}}//DONE: also save this number during coarsening, such that we do not have to recalc this every time
  // vector<double> toreturn;
  ///TODO TODO: implement radial basis fun based approach
  // toreturn.resize(nb_points_at_lvl[finer_lvl]);
  //
  // Size count=0;
  // std::map<Size,double> value_map;//maps points of coarser grid to their values
  // for (Size point=0; point<nb_points; point++)
  // {
  //   if (mask_list[coarser_lvl][point])//if point belongs to coarser grid
  //   {
  //   count++;
  //   value_map.insert(std::pair<Size,double>(point,to_interpolate[count]));
  //   }
  // }
  //
  // Size curr_count=0;//holds the index of the point for which we are currently interpolating (or just filling in)
  // for (Size point=0; point<nb_points; point++)
  // {
  // //if point in finer grid but not in coarser grid, we try to interpolate
  //   if (mask_list[finer_lvl][point])
  //   {
  //     if (!mask_list[coarser_lvl][point])
  //     {
  //       vector<Size> neighbors_in_coarse_grid;
  //       // std::set<Size> neighbors_in_coarse_grid;
  //       vector<Size> neighbors_of_point=neighbors_lists[finer_lvl].get_neighbors(point);
  //       vector<Size> neighbors_in_fine_grid;
  //       // std::set<Size> neighbors_in_fine_grid;//contains neighbors in fine grid, but not in coarse
  //       for (Size neighbor: neighbors_of_point)
  //       {
  //         if (mask_list[coarser_lvl][neighbor])
  //         {
  //           neighbors_in_coarse_grid.push_back(neighbor);
  //         }
  //         else
  //         {
  //           neighbors_in_fine_grid.push_back(neighbor);
  //         }
  //       }
  //
  //       while (neighbors_in_coarse_grid.empty())//please note that we would be very unlucky if we would get in this loop...
  //       {
  //         // std::set<Size> temp;//
  //         vector<Size> temp;
  //         for (Size fine_neighbor: neighbors_in_fine_grid)
  //         {
  //           vector<Size> neighbors_of_fine_neighbor=neighbors_lists[finer_lvl].get_neighbors(fine_neighbor);
  //           for (Size neighbor_of_neighbor: neighbors_of_fine_neighbor)
  //           {//if the neighbor of neighbor belongs to the coarse grid, just insert it
  //             if (value_map.count(neighbor_of_neighbor)>0)
  //             {neighbors_in_coarse_grid.push_back(neighbor_of_neighbor);}
  //             else          //TODO: add some protection such that we do not add the same neighbors again and again
  //             {temp.push_back(neighbor_of_neighbor);}
  //           }
  //         }
  //         neighbors_in_fine_grid=temp;
  //       }
  //       // std::set<vector<Size>> tetrahedra;
  //       vector<vector<Size>> tetrahedra;
  //       for (Size coarse_neighbor: neighbors_in_coarse_grid)
  //       {
  //         //For now, we just calculate all tetrahedra that have a point in neighbors_in_fine_grid
  //         //This is NOT EFFICIENT, but simpler to think about
  //         //Maybe TODO: first look at the points which have the maximal nb of neighbors in neighbors_in_fine_grid,
  //         // then construct tetrahedra using those points
  //         //Other possibility: add new datastructure which maps the points onto the tetrahedra
  //         // should be calculated after refining the grid
  //         std::set<vector<Size>> curr_tetrahedra=calc_all_tetra_with_point(point, coarser_lvl);//if we are using qhull or anything that gives us this, use that instead
  //         for (vector<Size> tetra_to_check: curr_tetrahedra)
  //         {//try to add all non-duplicates to the tetrahedra set
  //           tetrahedra.push_back(tetra_to_check);
  //         }
  //       }
  //       //TODO: using the tetrahedra we have found, we could also try to interpolate the neighboring points that we have found
  //       //  should be quite a bit faster. (if we do this, we also need a map that maps the points of the fine grid to the index for the output vector)
  //       //initialized with value much larger than anything we should encounter
  //       double curr_min_sum_abs=10000;//criterion for determining which tetrahedron is the best; is the sum of abs values of the barycentric coords (is >=1;should be minimized)
  //       double curr_best_interp=0;
  //       for (vector<Size> tetrahedron: tetrahedra)
  //       {
  //         Eigen::Vector<double,4> curr_coords=calc_barycentric_coords(tetrahedron, point);
  //         double curr_sum_abs=0;
  //         for (double coord: curr_coords)
  //         {curr_sum_abs+=std::abs(coord);}
  //         if (curr_sum_abs<curr_min_sum_abs)
  //         {
  //           Eigen::Vector<double,4> interp_values;
  //           interp_values << value_map[tetrahedron[0]],value_map[tetrahedron[1]],
  //                            value_map[tetrahedron[2]],value_map[tetrahedron[3]];
  //           curr_best_interp=interp_values.dot(curr_coords);
  //           //if all barycentric coordinates are greater than (or equal to) zero, this means that our point lies within the tetrahedron,
  //           // and therefore we have found the best interpolation. So we stop searching
  //           if (curr_coords[0]>=0&&curr_coords[1]>=0&&curr_coords[2]>=0&&curr_coords[3]>=0)
  //           {
  //             break;
  //           }
  //         }
  //       }
  //       toreturn[curr_count]=curr_best_interp;
  //     }
  //     else
  //     {
  //       toreturn[curr_count]=value_map[point];
  //     }//else, if in finer grid, still add to new vector
  //   curr_count++;
  //   }
  // }
  // return toreturn;
// }
