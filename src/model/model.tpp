#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cfloat>
#include <set>
#include "tools/types.hpp"
#include <limits>
#include <algorithm>
#include "voro++.hh"
using namespace voro;//FIXME: Do not use namespace, instead use voro:: for all functions of voro++

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


// Calculates the distance squared between two points
inline double Model::calc_distance2(Size point1,Size point2)
{
  return (geometry.points.position[point1]-geometry.points.position[point2]).squaredNorm();
}

// Calculates the normalized inner product between two directions (base to point1) and (base to point2)
// Better known as cosine similarity
inline double Model::calc_cosine(Size base, Size point1,Size point2)
{
  double length1=std::sqrt(calc_distance2(base,point1));
  double length2=std::sqrt(calc_distance2(base,point2));
  double innerprod=(geometry.points.position[point1]-geometry.points.position[base])
              .dot((geometry.points.position[point2]-geometry.points.position[base]));
  return innerprod/(length1*length2);
}

// Determines whether a point p2 lies further than p1 from point base when projected on (p1 base).
// note: should raname to 'lies far enough'
// Maybe TODO:make sure we first divide by length1, then apply the inner product//but for this we need a division operator for Vector3D
inline bool Model::lies_further_than(Size base, Size point1, Size point2)
{
  double length_squared=calc_distance2(base,point1);//length (p1 base) for normalization purposes
  double innerprod=(geometry.points.position[point1]-geometry.points.position[base])
              .dot((geometry.points.position[point2]-geometry.points.position[base]));
  double norm=innerprod/length_squared;
  return (norm>0.5);
}

// Calculates the relative difference of the density of point1 with respect to point2
//    @Returns: abs((d-d_other)/(d+d_other)) such that the smallest values will be assigned
inline double Model :: calc_diff_abundance_with_point(Size point1, Size point2)
{
    double abundance_self=chemistry.species.abundance[point1][1];
    double abundance_other=chemistry.species.abundance[point2][1];
    //std::cout<< "abundance self: " << abundance_self << std::endl;
    //std::cout<< "abundance other: " << abundance_other << std::endl;
    double temp_diff=std::abs((abundance_self-abundance_other)/(abundance_self+abundance_other));
    // std::cout<<"temp_diff:"<<temp_diff<<std::endl;
    return temp_diff;
}

// /// Returns function that determines whether two points are similar enough, according to the tolerance
// inline std::function<bool(Size,Size)> Model::points_are_similar(double tolerance)
// {
//   std::function<bool(Size,Size)> fun_to_return = [this,tolerance](Size p1, Size p2)
//   {
//     // std::cout<<"diff_abundance:"<<calc_diff_abundance_with_point(p1,p2)<<std::endl;
//     // std::cout<<"tolerance: "<<tolerance<<std::endl;
//     // std::cout<<"result: "<<(calc_diff_abundance_with_point(p1,p2)<tolerance)<<std::endl;
//     return (calc_diff_abundance_with_point(p1,p2)<tolerance);
//   };
//   return fun_to_return;
// }

/// Returns whether two points are similar enough, according to the tolerance
inline bool Model::points_are_similar(Size point1, Size point2, double tolerance)
{
    return (calc_diff_abundance_with_point(point1,point2)<tolerance);
}


// Coarsen the mesh,
// i.e. add another layer of coarsening.
inline void Model::coarsen(double tol)
{
    // auto mult=geometry.points.multiscale;
    geometry.points.multiscale.mask     .push_back(vector<bool>          (geometry.points.multiscale.mask     .back()));//should be deep copy
    geometry.points.multiscale.neighbors.push_back(vector<std::set<Size>>(geometry.points.multiscale.neighbors.back()));//should be deep copy
    std::set<Size> points_coarsened_around;
    geometry.points.multiscale.curr_coarsening_lvl=geometry.points.multiscale.get_max_coars_lvl();

    vector<Size> points_to_process=geometry.points.multiscale.get_current_points_in_grid();
    // vector<Size> next_points_to_process;
    // Size max_nb_neighbors=10;
    //
    // while(!points_to_process.empty())
    // {
      for (Size p: points_to_process)
      {
        if (can_be_coarsened(p, points_coarsened_around, tol))
        {//we want to prioritize coarsening around points with not too much neighbors
          // if (geometry.points.multiscale.get_nb_neighbors(p)<=max_nb_neighbors)
          // {
            coarsen_around_point(p);
          // std::cout << "Deleted around point: " << p <<std::endl;
            points_coarsened_around.insert(p);
          // }
          // else
          // {//if it has too much neighbors, try again in a next iteration
          //   next_points_to_process.push_back(p);
          // }
        }
      }
      // max_nb_neighbors*=2;
      // points_to_process=next_points_to_process;
      // next_points_to_process=vector<Size>();
    // }
}


// Returns whether the mesh at a point (p) can be coarsened.
inline bool Model::can_be_coarsened (const Size p, std::set<Size>& points_coarsened_around, double tol)
{
    // auto mult=geometry.points.multiscale;
    if (!geometry.points.multiscale.mask.back()[p]||!geometry.not_on_boundary(p)) {
      // std::cout<<"point on boundary or mask is set wrong?"<<std::endl;
      return false;}//if point no longer in grid, do not coarsen
    //if the point lies on the boundary, do not waste time trying to coarsen around it

    for (const Size n : geometry.points.multiscale.neighbors.back()[p])
    {
        // Do not coarsen if a neighbor was already coarsend at this level,
        // this avoids creating large holes in the mesh.
        // TODO: replace this....
        // if (!mask.back()[n]) {return false;}
        if (points_coarsened_around.find(n)!=points_coarsened_around.end()) {
          // std::cout<<"neighbor was already coarsened in this level"<<std::endl;
          return false;}

        // Do not coarsen if the required coarsening criterion does not hold.
        if(!points_are_similar(p,n,tol)) {
          // std::cout<<"points are too different"<<std::endl;
          return false;}
    }

    // Otherwise coarsen.
    return true;
}


inline void Model::coarsen_around_point (const Size p)
{
    // auto mult=geometry.points.multiscale;
    // New neighbors for p.
    // std::set<Size> new_neighbors;
    // Boundary neighbors need to be treated differently
    std::set<Size> boundary_neighbors;

    // Identify all neighbors with the current point (p),
    // i.e. remove its neighbors from the mesh by masking,
    // TODO: figure out whether it is ok/necessary to empty their neighbors
    //TODO: get the points to remove from somewhere else
    for (const Size n : geometry.points.multiscale.neighbors.back()[p])
    {
      if (geometry.not_on_boundary(n))//boundary points will NEVER get removed
        {
          geometry.points.multiscale.mask.back()[n] = false;
          //maps neighbor to point which deleted it (fallback plan for interpolation)
          geometry.points.multiscale.point_deleted_map.insert(std::pair<Size,Size>(n,p));
        }
    }

    // ..., and replace the neighbors of p (which are removed)
    // in their neighbors by p, using the symmetry of neighbors.
    // Also can contain non-deleted neighbors
    std::set<Size> neighbors_of_neighbors;

    // keeping track of which neighbors are deleted for each affected point
    // needed for reconnection purposes when deleting more than 1 neighbor from a neighbor of neighbor
    std::map<Size,vector<Size>> n_n_deleted_neighbors;
    // // The reverse map of the previous; this holds for each deleted point all n_n which are neighbors
    // std::map<Size,vector<Size>> n_n_of_deleted_points;

    for (const Size n : geometry.points.multiscale.neighbors.back()[p])
    {
      if (geometry.not_on_boundary(n))
      {
        // std::set<Size> curr_neighbors_of_n;//set of neighbors of n
        for (const Size n_n : geometry.points.multiscale.neighbors.back()[n])
        {
          if (geometry.points.multiscale.mask.back()[n_n]&&n_n!=p)//if neighbor of neighbor is still in the grid, i.e. not just a neighbor; also, do never try to add a point as its own neighbor!!!!
          { //for every deleted point, add its neighbors to its neighbors
            // for (Size other_n_n:curr_neighbors_of_n)
            // {
            //   geometry.points.multiscale.neighbors.back()[n_n].insert(other_n_n);
            //   geometry.points.multiscale.neighbors.back()[other_n_n].insert(n_n);
            // }
            // curr_neighbors_of_n.insert(n_n); //This leads to exponential growth of the number of neighbors in one single iteration//TODO:replace with post-processing
            //FIXME: find a way to select neighbors to add without increasing the total amount of neighbors of a n_n
            neighbors_of_neighbors.insert(n_n);
            geometry.points.multiscale.neighbors.back()[n_n].erase(n);//remove n from neighbors of n_n
            n_n_deleted_neighbors[n_n].push_back(n);
            // n_n_of_deleted_points[n].push_back(n_n);
          }
        }
        //for every deleted point, add its neighbors to its neighbors
        // for (const Size n_n:curr_neighbors_of_n)
        // {
        //   for (const Size other_n_n:curr_neighbors_of_n)
        //   {
        //     if (n_n<other_n_n)//add strict order to prevent adding the neighbors twice
        //     {
        //       geometry.points.multiscale.neighbors.back()[n_n].insert(other_n_n);
        //       geometry.points.multiscale.neighbors.back()[other_n_n].insert(n_n);
        //     }
        //   }
        // }
        geometry.points.multiscale.neighbors.back()[n]=std::set<Size>();//and finally also delete every neighbor of the deleted point
      }
      // also has become redundant
      else
      {//also do not forget to actually add our boundary elements as a neighbor of our point
        // new_neighbors.insert(n);
        boundary_neighbors.insert(n);
      }
    }

    //for now just the neighbors of neighbors of neighbors
    std::set<Size> container_points;
    for (const Size aff_point: neighbors_of_neighbors)
    {
      for (const Size n_n_n : geometry.points.multiscale.neighbors.back()[aff_point])
      {
        if (geometry.points.multiscale.mask.back()[n_n_n]&&n_n_n!=p)//add if neighbor of affected point is still in the grid, i.e. not just a neighbor; also, do never try to add a point as its own neighbor!!!!
        {
          container_points.insert(n_n_n);
        }
      }
    }
    //also add non-deleted neighbors to container_points
    for (Size bound_neigh:boundary_neighbors)
    {
      container_points.insert(bound_neigh);
    }

    //now also contains neighbors of neighbors and their neighbors
    container_points.insert(neighbors_of_neighbors.begin(),neighbors_of_neighbors.end());
    // std::merge(neighbors_of_neighbors, neighbors_of_neighbors.end(),neighbors_of_n_n,neighbors_of_n_n.end(),container_points.begin())
    // std::cout<<"Number of points in neighbors of neighbors: "<<neighbors_of_neighbors.size()<<std::endl;
    // std::cout<<"Number of points in container (-1): "<<container_points.size()<<std::endl;

    //kind of redundant, i guess
    // // Also deleting the neighbors (not on boundary) from the on boundary neighbors
    // for (const Size n:geometry.points.multiscale.neighbors.back()[p])
    // {
    //   if (geometry.not_on_boundary(n))
    //   {
    //     for (const Size b:boundary_neighbors)
    //     {
    //       geometry.points.multiscale.neighbors.back()[b].erase(n);
    //     }
    //   }
    // }


    // old method
    // // std::cout << "Size neighbors_of_neighbors: " << neighbors_of_neighbors.size() << std::endl;
    // for (const Size n_n:neighbors_of_neighbors)
    // {
    //   // Replace the removed points by p
    //   geometry.points.multiscale.neighbors.back()[n_n].insert(p);
    //
    //   // And add the others as new neighbors.
    //   new_neighbors.insert(n_n);
    // }

    //clear the old neighbors of point p
    geometry.points.multiscale.neighbors.back()[p]=std::set<Size>();

    //contains the non-deleted neighbors and the neighbors of deleted neighbors
    std::vector<Size> affected_points;
    affected_points.reserve(1+neighbors_of_neighbors.size());
    affected_points.push_back(p);
    // std::copy(neighbors_of_neighbors.begin(), neighbors_of_neighbors.end(), std::back_inserter(affected_points));
    std::copy(container_points.begin(), container_points.end(), std::back_inserter(affected_points));
    // std::set_union (neighbors_of_neighbors.begin(), neighbors_of_neighbors.end(), boundary_neighbors.begin(), boundary_neighbors.end(),
    //                 std::inserter(affected_points, affected_points.begin()));

    //correctly setting up a rectangular boundary for all affected points
    double x_min=std::numeric_limits<double>::max();
    double x_max=std::numeric_limits<double>::min();
    double y_min=std::numeric_limits<double>::max();
    double y_max=std::numeric_limits<double>::min();
    double z_min=std::numeric_limits<double>::max();
    double z_max=std::numeric_limits<double>::min();

    for (Size aff_point: affected_points)
    {
      Vector3D pos=geometry.points.position[aff_point];
      if (x_min>pos.x())
      {x_min=pos.x();}
      if (x_max<pos.x())
      {x_max=pos.x();}
      if (y_min>pos.y())
      {y_min=pos.y();}
      if (y_max<pos.y())
      {y_max=pos.y();}
      if (z_min>pos.z())
      {z_min=pos.z();}
      if (z_max<pos.z())
      {z_max=pos.z();}
    }
    //also making sure that all points lie strictly inside the boundary
    double deltax=x_max-x_min;
    double deltay=y_max-y_min;
    double deltaz=z_max-z_min;
    //TODO/FIXME: also make sure that no overflow can ever happen!!
    //+-1 is a temporary fix for also allowing 1D simulations
    x_max=x_max+0.001*deltax+1.0;
    x_min=x_min-0.001*deltax-1.0;
    y_max=y_max+0.001*deltay+1.0;
    y_min=y_min-0.001*deltay-1.0;
    z_max=z_max+0.001*deltaz+1.0;
    z_min=z_min-0.001*deltaz-1.0;
    std::cout<<"coarsening around point: "<<p<<std::endl;
    std::cout<<"container xlims: "<<x_min<<","<<x_max<<std::endl;
    std::cout<<"container ylims: "<<y_min<<","<<y_max<<std::endl;
    std::cout<<"container zlims: "<<z_min<<","<<z_max<<std::endl;

    container con(x_min,x_max,y_min,y_max,z_min,z_max,8,8,8,
                  			false,false,false,8);

    particle_order p_order;
    voronoicell_neighbor cell;

  	// add particles into the container, starting with point p
  	for(Size aff_point:affected_points)
    { Vector3D pos=geometry.points.position[aff_point];
      con.put(p_order,static_cast<int>(aff_point),pos.x(),pos.y(),pos.z());
      std::cout<<"aff point: "<<aff_point<<std::endl;
      std::cout<<"(x,y,z)=("<<pos.x()<<","<<pos.y()<<","<<pos.z()<<")"<<std::endl;
    }

    c_loop_order l_order(con,p_order);

    //back to the beginning
    l_order.start();
    // std::cout<<"currently computing neighbors around point: "<<l_order.pid()<<std::endl;
    //for point p, just add all found neighbors to its neighbors list
    con.compute_cell(cell,l_order);
    std::vector<int> found_neighbors;
    cell.neighbors(found_neighbors);
    for (int new_neighbor:found_neighbors)
    {
      //note: voro++ return negative values for the boundaries//not that it should play a role here
      if (new_neighbor>=0)
      { //and finally adding new neighbors to point p
        geometry.points.multiscale.neighbors.back()[p].insert(static_cast<Size>(new_neighbor));
        geometry.points.multiscale.neighbors.back()[static_cast<Size>(new_neighbor)].insert(p);
      }
    }
    //use inc() for incrementing (returns false when not finding a next one)
    if(l_order.inc())//this if-statement should always return true
    {
      do {
      Size current_point=static_cast<Size>(l_order.pid());
      if (neighbors_of_neighbors.find(current_point)==neighbors_of_neighbors.end())
      { //we only care about finding the neighbors of the neighbors of neighbors
        continue;
      }
      con.compute_cell(cell,l_order);
      std::vector<int> found_neighbors;
      cell.neighbors(found_neighbors);

      // std::cout<<"computed partial neighbors around: "<<l_order.pid()<<std::endl;
      //note:negative values refer to the boundaries
      for (int found_neighbor:found_neighbors)
      {
        if (found_neighbor>=0)//voro++ denotes the walls (which we do not need) by negative values
        {
          Size size_found_neighbor=static_cast<Size>(found_neighbor);
          // if (size_found_neighbor>current_point)//in order to not try to add neighbors twice
          // { //check if already a neighbor (thus no need to add again)
          if (geometry.points.multiscale.neighbors.back()[current_point].find(size_found_neighbor) ==
              geometry.points.multiscale.neighbors.back()[current_point].end())
          {
              // for (Size deleted_neighbor:n_n_deleted_neighbors[current_point])
              // { // normally, one should also put the 'neighbors of the neighbors of the neighbors' in the container to accurately find the neighbors of the 'neighbors of the neighbors'.
                //but as we do not change the n_n_n, we can instead check if the new proposed neighbors somehow replace the old deleted neighbors
                //without this check, we might accidentally also include too much neighbors.
                // if (lies_further_than(current_point, deleted_neighbor, size_found_neighbor))
                // {
                  //and finally adding new neighbors to point p
            geometry.points.multiscale.neighbors.back()[current_point].insert(size_found_neighbor);
            geometry.points.multiscale.neighbors.back()[size_found_neighbor].insert(current_point);
                  // break;
                // }
              // }
          }
          // }
        }
      }

      } while(l_order.inc());
    }

    // old method
    // // Finally also check whether neighbors of neighbors need some extra neighbors
    // for (Size n_n:neighbors_of_neighbors)
    // {
    //   Size nb_deleted_neighbors=n_n_deleted_neighbors[n_n].size();
    //   // if (nb_deleted_neighbors>1)//we do not want the amount of neighbors around points to decrease
    //   {
    //     // Determine the deleted neighbor closest to p (when looking at directions)
    //     double maxcos=-1;//maximum cosine similarity
    //     Size maxindex=0;//corresponding index in n_n_deleted_neighbors[n_n]
    //     for (Size idx=0;idx<n_n_deleted_neighbors[n_n].size();idx++)
    //     {
    //       double cossim=calc_cosine(n_n,p,n_n_deleted_neighbors[n_n][idx]);
    //       if (cossim>maxcos)
    //       {
    //         maxcos=cossim;
    //         maxindex=idx;
    //       }
    //     }
    //     // Testing commenting this out
    //     // n_n_deleted_neighbors[n_n].erase(n_n_deleted_neighbors[n_n].begin() + maxindex);
    //     // For all other deleted points, calculate the closest n_n (in direction; calculate using the cosine)
    //     for (Size deleted_neighbor: n_n_deleted_neighbors[n_n])
    //     {
    //       double maxcos=-1;//maximum cosine similarity
    //       Size opt_n_n=parameters.npoints();//corresponding optimal neighbor of neighbor
    //       // vector<Size> n_n_of_deleted_neighbor=n_n_of_deleted_points[deleted_neighbor];
    //       // Size vectorsize=neighbors_of_neighbors.size();
    //       Size vectorsize=affected_points.size();
    //       // if (vectorsize==1)
    //       // {
    //       //   std::cout<<"Wait, this deleted point only had two non-coarsened neighbors. Skipping for now."<<std::endl;
    //       //   break;//TODO?? throw error or resolve this any other way.
    //       // }
    //       // for (Size temp_n_n:neighbors_of_neighbors)
    //       for (Size temp_n_n:affected_points)
    //       {
    //         if (temp_n_n!=n_n)//do not try to add a point as its own neighbor//also division by zero would ensue
    //         {
    //           double cossim=calc_cosine(n_n,p,temp_n_n);
    //           if (cossim>maxcos)
    //           {
    //             maxcos=cossim;
    //             opt_n_n=temp_n_n;
    //           }
    //         }
    //       }
    //       // Size point_to_insert=neighbors_of_neighbors[maxindex];
    //       geometry.points.multiscale.neighbors.back()[n_n].insert(opt_n_n);
    //       geometry.points.multiscale.neighbors.back()[opt_n_n].insert(n_n);
    //     }
    //
    //   }
    // }
    //
    // // Set the new nearest neighbors.
    // geometry.points.multiscale.neighbors.back()[p] = new_neighbors;
}


/// Initializes the multigrid
///   @Parameter[in]: min_nb_points: stop coarsening when number of points of finest layers is equal or less than this
///   @Parameter[in]: max_coars_lvl: The maximum coarsening level we allow
///   @Parameter[in]: tol: the tolerance level for which we still consider points similar enough
inline int Model::setup_multigrid(Size min_nb_points, Size max_coars_lvl, double tol)
{
  //set number of off-diagonal elements in lambda matrix 0; needed because we will interpolate these values
  parameters.n_off_diag=0;
  // std::cout<<"Tolerance: "<<tol<<std::endl;
  // std::function<bool(Size,Size)> fun_to_del=points_are_similar(tol);//function that says if two points are similar enough
  // geometry.points.multiscale.set_not_on_boundary_fun([&](Size p){return geometry.not_on_boundary(p);});//function that says whether a point lies on the boundary
  // geometry.points.multiscale.set_comparison_fun(fun_to_del);
  //first, we coarsen the grid until we either have too few points left or have too many coarsening levels
  while((geometry.points.multiscale.get_max_coars_lvl()<max_coars_lvl)&&
  (geometry.points.multiscale.get_total_points(geometry.points.multiscale.get_max_coars_lvl())>min_nb_points))
  {std::cout<<"coarsening layer"<<std::endl;
  // std::cout<<"max coarsening lvl: "<<geometry.points.multiscale.get_max_coars_lvl()<<std::endl;
  // std::cout<<"curr coarsening lvl: "<<geometry.points.multiscale.get_curr_coars_lvl()<<std::endl;
    coarsen(tol);
    std::cout<<"coarsened layer; number points remaining: "<<geometry.points.multiscale.get_total_points(geometry.points.multiscale.get_max_coars_lvl())<<std::endl;
  }
  std::cout<<"finished coarsening"<<std::endl;
  return (0);
}

/// Computes second order feautrier using multigrid //never mind, this is useless
inline int Model::compute_feautrier_order_2_multigrid()
{
  //set curr coarsening level to max
  geometry.points.multiscale.set_curr_coars_lvl(geometry.points.multiscale.get_max_coars_lvl());
  //solve and interpolate J for all coarser grids
  while(geometry.points.multiscale.get_curr_coars_lvl()>0)
  {
    compute_radiation_field_feautrier_order_2();
    std::cout<<"Probably crashes at interpolation"<<std::endl;
    //for all frequencies, interpolate J
    interpolate_matrix_local(geometry.points.multiscale.get_curr_coars_lvl(),radiation.J);
    std::cout<<"Did not crash at interpolation"<<std::endl;
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

/// FIXME: either delete me or apply same fixes as in interpolate_matrix_local!!!!
///DO NOT USE IN THIS STATE!!!!!
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


// /// NOTE: with some better data structure, we can remove this duplication
/// Interpolates the level populations (linearized matrix) from the coarser level to the first finer level
///   @Parameter [in] coarser_lvl: the coarsening level of the coarser grid
///   @Parameter [in/out] to_interpolate: the paracabs matrix with values at the points of the coarser grid (and some irrelevant values inbetween), has length equal to the total number of points
// template <Real T>//T should be double, long or long double
inline void Model::interpolate_levelpops_local(Size coarser_lvl)
{
  //TODO: only use large vector and some masks
  //we cannot interpolate from lvl 0 to some level -..., so just return
  if (coarser_lvl==0)
  {return;}
  std::cout<<"Starting the interpolation"<<std::endl;
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
    //TODO: also use neighbors of neighbors as better interpolation
    // Note: using the current multigrid creation method, the number of neighbors in the coarse grid is almost always at least 1
    if (neighbors_coarser_grid.size()==0)//this shouldn't be possible, but you never know what goes wrong
    {
      std::cout<<"No neighbors for the current point! Using point which deleted it instead as neighbor."<<std::endl;
      std::cout<<"Current point: "<<diff_point<<std::endl;
      std::cout<<"Point which replaces it: "<<geometry.points.multiscale.point_deleted_map.at(diff_point)<<std::endl;
      neighbors_coarser_grid.push_back(geometry.points.multiscale.point_deleted_map.at(diff_point));
      std::cout<<"Succesfully replaced neighbors"<<std::endl;
    }
    ///commented out for now until we solve the issue with having way too much neighbors if we apply this...
    // In the case we do not really have enough points for a good interpolation, just add more
    if (neighbors_coarser_grid.size()<MIN_INTERPOLATION_POINTS)
    {//because this estimate is just beyond terrible: also add neighbors of neighbors
      // for (Size freqidx=0; freqidx<parameters.nfreqs(); freqidx++)
      // {
      std::set<Size> temp_neighbors_coarser_grid;
      for (Size neighbor_coarse: neighbors_coarser_grid)
      {
        for (Size neighbor_of_coarse_neighbor: geometry.points.multiscale.get_neighbors(neighbor_coarse,coarser_lvl))
        { //just make sure that we do not accidentally insert a point we already have (in neighbors_coarser_grid)
          if (std::find(neighbors_coarser_grid.begin(), neighbors_coarser_grid.end(), neighbor_of_coarse_neighbor) == neighbors_coarser_grid.end())
          {
            temp_neighbors_coarser_grid.insert(neighbor_of_coarse_neighbor);
          }
        }
        // }
      // std::cout<<"Just taking the value of the closest point"<<std::endl;
      // std::cout<<"Interpolation values"<<to_interpolate(neighbors_coarser_grid[0],freqidx)<<std::endl;
      // to_interpolate[diff_point,freqidx]=to_interpolate(neighbors_coarser_grid[0],freqidx);
      }
      //Finally add our new points too
      std::copy(temp_neighbors_coarser_grid.begin(), temp_neighbors_coarser_grid.end(), std::back_inserter(neighbors_coarser_grid));
    }
    if (neighbors_coarser_grid.size()>MAX_INTERPOLATION_POINTS)
    {//now this just becomes: 1) to expensive to calculate and 2)difference between distances might become too large
      /// This gives too much spam, commenting for now
      // std::cout<<"far too many neighbors to interpolate: "<<neighbors_coarser_grid.size()<<std::endl;
      // std::cout<<"using the first 50 neighbors instead to interpolate"<<std::endl;
      Size maxnbneighbors=MAX_INTERPOLATION_POINTS;
      Size nbneighbors=neighbors_coarser_grid.size();
      vector<double> distances2;//stores the distances2 of the n closest points
      distances2.resize(maxnbneighbors);
      std::fill(distances2.begin(), distances2.end(), std::numeric_limits<double>::max());
      vector<Size> closest_points;//stores the n closest points
      double maxdist=std::numeric_limits<double>::max();//maximum distance of the n closest points
      Size maxindex=0;//corresponding index
      closest_points.resize(maxnbneighbors);
      std::fill(closest_points.begin(), closest_points.end(), parameters.npoints());
      for (Size idx=0;idx<nbneighbors;idx++)
      {
        double tempdistance=(geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[diff_point]).squaredNorm();
        if (tempdistance<maxdist)
        {//replace the curr max distance
          closest_points[maxindex]=neighbors_coarser_grid[idx];
          distances2[maxindex]=tempdistance;
          //get new max distance and corresponding max index
          auto tempindex = std::max_element(distances2.begin(), distances2.end());
          maxindex=tempindex - distances2.begin();
          maxdist=distances2[maxindex];
        }
      }
      // vector<double> distances2=(geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[diff_point]).squaredNorm());

      // vector<Size> subvector = indices;
      //Maybe TODO: find a better way to choose the best neighbors, maybe just order them on distance
      neighbors_coarser_grid=closest_points;
    }
    // cout<<"number neighbors to interpolate from: "<<neighbors_coarser_grid.size()<<std::endl;
    // else
    // {//use rbf estimate
    Size nb_neighbors_coarser_grid=neighbors_coarser_grid.size();

    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> rbf_mat(nb_neighbors_coarser_grid, nb_neighbors_coarser_grid);
    Eigen::Matrix<double,1,Eigen::Dynamic> distance_with_neighbors(1,nb_neighbors_coarser_grid);

    for (Size idx=0; idx<nb_neighbors_coarser_grid; idx++)
    {
      // std::cout<<"calc dist with point"<<std::endl;
      distance_with_neighbors(idx)=std::sqrt((geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[diff_point]).squaredNorm());
      rbf_mat(idx,idx)=0;//distance between point and itself is zero
      for (Size idx2=0; idx2<idx; idx2++)
      {
        //calculate radius
        double radius=std::sqrt((geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[neighbors_coarser_grid[idx2]]).squaredNorm());
        rbf_mat(idx,idx2)=radius;
        rbf_mat(idx2,idx)=radius;
      }
    }
    double mindist=distance_with_neighbors.minCoeff();
    mindist=mindist*5.0;
    // double maxdist=rbf_mat.maxCoeff()*5.0;//arbitrary number 5 to make the max_dist larger
    rbf_mat=rbf_mat/mindist;
    rbf_mat=rbf_mat.unaryExpr(std::ptr_fun(rbf_local<double>));
    distance_with_neighbors=distance_with_neighbors/mindist;
    distance_with_neighbors=distance_with_neighbors.unaryExpr(std::ptr_fun(rbf_local<double>));
    // std::cout<<"determinant: "<<distance_with_neighbors.determinant()<<std::endl;
    //going with more robust LDLT decomposition in the hope that this reduces the number of nans generated
    Eigen::LDLT<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> ldltdec(rbf_mat);
    //now that we have our llt decomposition, calculate the interpolated value
    for (Size specidx=0; specidx<parameters.nlspecs(); specidx++)
    {
      for (Size levidx=0; levidx<lines.lineProducingSpecies[specidx].linedata.nlev; levidx++)
      {
        Eigen::Vector<double,Eigen::Dynamic> right_hand_side(nb_neighbors_coarser_grid);
        for (Size idx=0; idx<nb_neighbors_coarser_grid; idx++)
        {
          // std::cout<<"Interpolation values"<<to_interpolate(idx,freqidx)<<std::endl;
          // interpolating the logarithm to avoid any semblance of negative numbers as result
          right_hand_side(idx)=static_cast<double>(lines.lineProducingSpecies[specidx].get_level_pop(neighbors_coarser_grid[idx],levidx));
        }
        Eigen::Vector<double,Eigen::Dynamic> weights=ldltdec.solve(right_hand_side);
        double interpolated_value=(distance_with_neighbors*weights)(0,0);
        // std::cout<<"interpolated value: "<<interpolated_value<<std::endl;
        // sometimes, this procedure can give negative values (which we do not want),
        // FIXME: at least add a warning when this happens
        lines.lineProducingSpecies[specidx].set_level_pop(diff_point,levidx,abs(static_cast<Real>(interpolated_value)));
        if (std::isnan(interpolated_value))
        {
          std::cout<<"Something went wrong during interpolating: nan value occuring"<<std::endl;
          throw std::runtime_error("Nan encountered during interpolation");
        }
      }
    }
    // }

  }
}


/// FIXME: use doubles for radius instead of T !!!// or just think about just enforcing doubles instead...
/// Interpolates the paracabs matrix from the coarser level to the first finer level
///   @Parameter [in] coarser_lvl: the coarsening level of the coarser grid
///   @Parameter [in/out] to_interpolate: the paracabs matrix with values at the points of the coarser grid (and some irrelevant values inbetween), has length equal to the total number of points
template <typename T>//T should be double, long or long double
inline void Model::interpolate_matrix_local(Size coarser_lvl, Matrix<T> &to_interpolate)
{//TODO: only use large vector and some masks
  //we cannot interpolate from lvl 0 to some level -..., so just return
  if (coarser_lvl==0)
  {return;}
  std::cout<<"Starting the interpolation"<<std::endl;
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
    //TODO: also use neighbors of neighbors as better interpolation
    // Note: using the current multigrid creation method, the number of neighbors in the coarse grid is almost always at least 1
    if (neighbors_coarser_grid.size()==0)//this shouldn't be possible, but you never know what goes wrong
    {
      std::cout<<"No neighbors for the current point! Using point which deleted it instead as neighbor."<<std::endl;
      std::cout<<"Current point: "<<diff_point<<std::endl;
      std::cout<<"Point which replaces it: "<<geometry.points.multiscale.point_deleted_map.at(diff_point)<<std::endl;
      neighbors_coarser_grid.push_back(geometry.points.multiscale.point_deleted_map.at(diff_point));
      std::cout<<"Succesfully replaced neighbors"<<std::endl;
    }
    ///commented out for now until we solve the issue with having way too much neighbors if we apply this...
    // In the case we do not really have enough points for a good interpolation, just add more
    if (neighbors_coarser_grid.size()<MIN_INTERPOLATION_POINTS)
    {//because this estimate is just beyond terrible: also add neighbors of neighbors
      // for (Size freqidx=0; freqidx<parameters.nfreqs(); freqidx++)
      // {
      std::set<Size> temp_neighbors_coarser_grid;
      for (Size neighbor_coarse: neighbors_coarser_grid)
      {
        for (Size neighbor_of_coarse_neighbor: geometry.points.multiscale.get_neighbors(neighbor_coarse,coarser_lvl))
        { //just make sure that we do not accidentally insert a point we already have (in neighbors_coarser_grid)
          if (std::find(neighbors_coarser_grid.begin(), neighbors_coarser_grid.end(), neighbor_of_coarse_neighbor) == neighbors_coarser_grid.end())
          {
            temp_neighbors_coarser_grid.insert(neighbor_of_coarse_neighbor);
          }
        }
        // }
      // std::cout<<"Just taking the value of the closest point"<<std::endl;
      // std::cout<<"Interpolation values"<<to_interpolate(neighbors_coarser_grid[0],freqidx)<<std::endl;
      // to_interpolate[diff_point,freqidx]=to_interpolate(neighbors_coarser_grid[0],freqidx);
      }
      //Finally add our new points too
      std::copy(temp_neighbors_coarser_grid.begin(), temp_neighbors_coarser_grid.end(), std::back_inserter(neighbors_coarser_grid));
    }
    if (neighbors_coarser_grid.size()>MAX_INTERPOLATION_POINTS)
    {//now this just becomes: 1) to expensive to calculate and 2)difference between distances might become too large
      /// This gives too much spam, commenting for now
      // std::cout<<"far too many neighbors to interpolate: "<<neighbors_coarser_grid.size()<<std::endl;
      // std::cout<<"using the first 50 neighbors instead to interpolate"<<std::endl;
      Size maxnbneighbors=MAX_INTERPOLATION_POINTS;
      Size nbneighbors=neighbors_coarser_grid.size();
      vector<double> distances2;//stores the distances2 of the n closest points
      distances2.resize(maxnbneighbors);
      std::fill(distances2.begin(), distances2.end(), std::numeric_limits<double>::max());
      vector<Size> closest_points;//stores the n closest points
      double maxdist=std::numeric_limits<double>::max();//maximum distance of the n closest points
      Size maxindex=0;//corresponding index
      closest_points.resize(maxnbneighbors);
      std::fill(closest_points.begin(), closest_points.end(), parameters.npoints());
      for (Size idx=0;idx<nbneighbors;idx++)
      {
        double tempdistance=(geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[diff_point]).squaredNorm();
        if (tempdistance<maxdist)
        {//replace the curr max distance
          closest_points[maxindex]=neighbors_coarser_grid[idx];
          distances2[maxindex]=tempdistance;
          //get new max distance and corresponding max index
          auto tempindex = std::max_element(distances2.begin(), distances2.end());
          maxindex=tempindex - distances2.begin();
          maxdist=distances2[maxindex];
        }
      }
      // vector<double> distances2=(geometry.points.position[neighbors_coarser_grid[idx]]-geometry.points.position[diff_point]).squaredNorm());

      // vector<Size> subvector = indices;
      //Maybe TODO: find a better way to choose the best neighbors, maybe just order them on distance
      neighbors_coarser_grid=closest_points;
    }
    // cout<<"number neighbors to interpolate from: "<<neighbors_coarser_grid.size()<<std::endl;
    // else
    // {//use rbf estimate
    Size nb_neighbors_coarser_grid=neighbors_coarser_grid.size();

    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> rbf_mat(nb_neighbors_coarser_grid, nb_neighbors_coarser_grid);
    Eigen::Matrix<T,1,Eigen::Dynamic> distance_with_neighbors(1,nb_neighbors_coarser_grid);

    for (Size idx=0; idx<nb_neighbors_coarser_grid; idx++)
    {
      // std::cout<<"calc dist with point"<<std::endl;
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
    T mindist=distance_with_neighbors.minCoeff();
    mindist=mindist*5.0;
    // T maxdist=rbf_mat.maxCoeff()*5.0;//arbitrary number 5 to make the max_dist larger
    rbf_mat=rbf_mat/mindist;
    rbf_mat=rbf_mat.unaryExpr(std::ptr_fun(rbf_local<T>));
    distance_with_neighbors=distance_with_neighbors/mindist;
    distance_with_neighbors=distance_with_neighbors.unaryExpr(std::ptr_fun(rbf_local<T>));
    // std::cout<<"determinant: "<<distance_with_neighbors.determinant()<<std::endl;
    //going with more robust LDLT decomposition in the hope that this reduces the number of nans generated
    Eigen::LDLT<Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>> ldltdec(rbf_mat);
    //now that we have our llt decomposition, calculate the interpolated value
    for (Size freqidx=0; freqidx<parameters.nfreqs(); freqidx++)
    {
      Eigen::Vector<T,Eigen::Dynamic> right_hand_side(nb_neighbors_coarser_grid);
      for (Size idx=0; idx<nb_neighbors_coarser_grid; idx++)
      {
        // std::cout<<"Interpolation values"<<to_interpolate(idx,freqidx)<<std::endl;
        // interpolating the logarithm to avoid any semblance of negative numbers as result
        right_hand_side(idx)=to_interpolate(neighbors_coarser_grid[idx],freqidx);
      }
      Eigen::Vector<T,Eigen::Dynamic> weights=ldltdec.solve(right_hand_side);
      T interpolated_value=(distance_with_neighbors*weights)(0,0);
      // std::cout<<"interpolated value: "<<interpolated_value<<std::endl;
      // sometimes, this procedure can give negative values (which we do not want),
      to_interpolate(diff_point,freqidx)=abs(interpolated_value);
      if (std::isnan(interpolated_value))
      {
        std::cout<<"Something went wrong during interpolating: nan value occuring"<<std::endl;
        throw std::runtime_error("Nan encountered during interpolation");
      }
    }
    // }

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
