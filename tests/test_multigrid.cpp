#include <iostream>
using std::cout;
using std::endl;

#include "io/python/io_python.hpp"
// #include "model/model.cpp"
#include "solver/solver.hpp"
#include <set>
#include "tools/timer.hpp"


//__global__
//void kernel ()
//{
//    return;
//}

//inline void test (Geometry &geometry)
//{
//    geometry.points.position.vec[0].print();
//    geometry.points.position.vec[1].print();
//
//    Geometry* geometry_copy = (Geometry*) pc::accelerator::malloc(sizeof(Geometry));
//    pc::accelerator::memcpy_to_accelerator (geometry_copy, &geometry, sizeof(Geometry));
//    kernel<<<1,1>>> (geometry_copy);
//
//    points.position.vec[0].print();
//    points.position.copy_ptr_to_vec ();
//    points.position.vec[0].print();
//}


int main (int argc, char **argv)
{
    cout << "Running test_multigrid..." << endl;

    /// Store model name
    const string modelName = argv[1];

    cout << modelName << endl;

    pc::accelerator::list_accelerators();

    IoPython io ("hdf5", modelName);


    cout << "sizeof Model    = " << sizeof(Model)    << endl;
    cout << "sizeof geometry = " << sizeof(Geometry) << endl;
    cout << "sizeof points   = " << sizeof(Points)   << endl;
    cout << "sizeof Vector3D = " << sizeof(Vector3D) << endl;

    cout << "sizeof Real = " <<sizeof(Real) << endl;

    cout << "n threads = " << paracabs::multi_threading::n_threads_avail() << endl;
    // paracabs::multi_threading::set_n_threads_avail(6);//set nb threads to one for debugging output
    cout << "n threads = " << paracabs::multi_threading::n_threads_avail() << endl;


    Model model;
    model.read (io);
    Parameters parameters=model.parameters;
    model.compute_spectral_discretisation();
    model.compute_LTE_level_populations();
    model.compute_inverse_line_widths();
    cout << "model read"<<endl;




    // cout <<"solving without multigrid"<<endl;
    // model.compute_level_populations(true,100);
    //NOTE TO SELF: do NOT every try to use a ridiculous amount of coarsening: if only boundary points are left, the interpolation part will probably be hell
    //mgImplementation choices: 1:"NaiveMG", 2:"VCycle", 3:"WCycle"
    model.setup_multigrid(10, 2, 0.1, 2, 20);
    cout << "setup multigrid" << endl;

    std::cout<<"checking symmetry of neighbors"<<std::endl;
    vector<Size> current_points_in_grid=model.geometry.points.multiscale.get_current_points_in_grid();
    for (Size point_to_check: current_points_in_grid)
    {
      std::set<Size> neighbors_of_point=model.geometry.points.multiscale.get_neighbors(point_to_check);
      for (Size neighbor:neighbors_of_point)
      {
        std::set<Size> neighbors_of_neighbor=model.geometry.points.multiscale.get_neighbors(neighbor);
        if (neighbors_of_neighbor.find(point_to_check)==neighbors_of_neighbor.end())
        {
          std::cout<<"Point "<<point_to_check<<"is not a neighbor of "<<neighbor<<std::endl;
        }
      }
    }

    std::cout<<"Checking whether a point contains itself as neighbor"<<std::endl;
    for (Size point_to_check: current_points_in_grid)
    {
      std::set<Size> neighbors_of_point=model.geometry.points.multiscale.get_neighbors(point_to_check);
      // for (Size neighbor:neighbors_of_point)
      // {
      //   std::set<Size> neighbors_of_neighbor=model.geometry.points.multiscale.get_neighbors(neighbor);
      if (neighbors_of_point.find(point_to_check)!=neighbors_of_point.end())
      {
        std::cout<<"Point "<<point_to_check<<"is a neighbor of itself."<<std::endl;
      }
      // }
    }

    model.writing_populations_to_disk=true;

    model.compute_level_populations_multigrid(true);
    // vector<Size> current_points_in_grid=model.geometry.points.multiscale.get_current_points_in_grid();
    // for (Size idx=0;idx<current_points_in_grid.size();idx++)
    // {
    //   Size nb_neighbors=model.geometry.points.multiscale.get_nb_neighbors(current_points_in_grid[idx]);
    //   cout<<nb_neighbors<<std::endl;
    //
    //   // std::set<Size>neighbors=model.geometry.points.multiscale.get_neighbors(point);
    //   // for (Size n:neighbors){
    //   //   cout<<n<<endl;
    //   // }
    // }

    // Size nb_boundary_points=0;
    // for (Size p=0;p<parameters.npoints();p++)
    // {
    //   //check whether all boundary points have neighbors
    //   if (!model.geometry.not_on_boundary(p))
    //   {
    //     nb_boundary_points++;
    //     if (model.geometry.points.multiscale.get_nb_neighbors(p)==0)
    //     {
    //       std::cout<<"Boundary point: "<<p<<" does not have any neighbors"<<std::endl;
    //     }
    //   }
    // }
    // cout<<"done checking boundary points"<<endl;
    // cout<<"number boundary points left: "<<nb_boundary_points<<endl;

    // for (Size p=0;p<parameters.npoints();p++)
    // {
    //   std::cout<<"point: "<<p<<" has number neighbors: "<<model.geometry.points.multiscale.get_nb_neighbors(p)<<std::endl;
    //   for (Size n: model.geometry.points.multiscale.get_neighbors(p,0))
    //   {
    //     std::cout<<"Distance: "<<std::sqrt((model.geometry.points.position[p]-model.geometry.points.position[n]).squaredNorm())<<std::endl;
    //   }
    // }


    // std::set<Size> test_neighbors=model.geometry.points.multiscale.get_neighbors(152970,0);
    // cout<<"neighbors of point 152970 in finest grid:"<<endl;
    // for (Size fine_neighbor:test_neighbors)
    // {
    //   cout<<fine_neighbor<<" is part of coarse grid?: "<<model.geometry.points.multiscale.get_mask(2)[fine_neighbor]<<endl;
    // }
    //
    // std::cout<<"Is point 152970 still in coarse grid? "<<(model.geometry.points.multiscale.get_mask(1)[152970])<<std::endl;
    //
    // test_neighbors=model.geometry.points.multiscale.get_neighbors(152970,1);
    // cout<<"neighbors of point 152970 in coarser grid:"<<endl;
    // for (Size coarse_neighbor:test_neighbors)
    // {
    //   cout<<coarse_neighbor<<" is part of coarser grid."<<endl;
    // }
    //
    // std::cout<<"Is point 152970 still in coarsest grid? "<<(model.geometry.points.multiscale.get_mask(2)[152970])<<std::endl;
    // test_neighbors=model.geometry.points.multiscale.get_neighbors(152970,2);
    // cout<<"neighbors of point 152970 in coarsest grid:"<<endl;
    // for (Size coarse_neighbor:test_neighbors)
    // {
    //   cout<<coarse_neighbor<<" is part of coarsest grid."<<endl;
    // }
    //
    //
    //
    // std::set<Size> test_neighbors2=model.geometry.points.multiscale.get_neighbors(152833,0);
    // cout<<"neighbors of point 152833 in finest grid:"<<endl;
    // for (Size fine_neighbor:test_neighbors2)
    // {
    //   cout<<fine_neighbor<<" is part of coarse grid?: "<<model.geometry.points.multiscale.get_mask(2)[fine_neighbor]<<endl;
    // }
    //
    // std::cout<<"Is point 152833 still in coarse grid? "<<(model.geometry.points.multiscale.get_mask(1)[152833])<<std::endl;
    //
    // test_neighbors2=model.geometry.points.multiscale.get_neighbors(152833,1);
    // cout<<"neighbors of point 152833 in coarser grid:"<<endl;
    // for (Size coarse_neighbor:test_neighbors2)
    // {
    //   cout<<coarse_neighbor<<" is part of coarser grid."<<endl;
    // }
    //
    // std::cout<<"Is point 152833 still in coarsest grid? "<<(model.geometry.points.multiscale.get_mask(2)[152833])<<std::endl;
    // test_neighbors2=model.geometry.points.multiscale.get_neighbors(152833,2);
    // cout<<"neighbors of point 152833 in coarsest grid:"<<endl;
    // for (Size coarse_neighbor:test_neighbors2)
    // {
    //   cout<<coarse_neighbor<<" is part of coarsest grid."<<endl;
    // }



    // model.compute_level_populations (true,100);

    // model.compute_level_populations_multigrid(true, 100);


    // auto fun_to_del=model.points_are_similar(0.1);
    // model.geometry.points.multiscale.set_not_on_boundary_fun([&](Size p){return model.geometry.not_on_boundary(p);});
    // model.geometry.points.multiscale.set_comparison_fun(fun_to_del);
    // model.geometry.points.multiscale.coarsen();
    //
    // model.interpolate_matrix_local(1,model.radiation.J);
    // cout << "interpolation worked" << endl;


    // model.compute_spectral_discretisation();
    // model.compute_LTE_level_populations();
    // model.compute_inverse_line_widths();
    // //just testing to see how long it takes
    // Timer timer("solver: 2nd order Feautrier");
    // timer.start();
    // model.compute_radiation_field_feautrier_order_2();
    // timer.stop();
    //
    // //
    // auto fun_to_del=model.points_are_similar(0.1);
    // model.geometry.points.multiscale.set_not_on_boundary_fun([&](Size p){return model.geometry.not_on_boundary(p);});
    // model.geometry.points.multiscale.set_comparison_fun(fun_to_del);
    // model.geometry.points.multiscale.coarsen();
    //
    // Timer timer2("solver: multigrid implementation");
    // timer2.start();
    // model.compute_radiation_field_feautrier_order_2();
    // timer2.stop();


    //TODO also add interpolation and another hime the computation

    // vector<Size> old_neighbors=model.geometry.points.curr_neighbors.get_neighbors(0);
    // std::set<Size> neighbors=model.geometry.points.multiscale.get_neighbors(0,0);
    //
    //
    // for (Size nb:neighbors)
    // {
    // cout << nb << endl;
    // }
    // cout << "compare with" << endl;
    //
    // for (Size nb:old_neighbors)
    // {
    // cout << nb << endl;
    // }
    // cout << "now testing coarsening" << endl;
    // cout << "no points deleted" << endl;
    // //lambda expression that always returns false, just to see whether it crashes (or not)
    // auto alwaysfalse = [](Size p1, Size p2)
    // {
    //     return false;
    // };
    // cout << "here" << endl;
    // model.geometry.points.multiscale.set_comparison_fun(alwaysfalse);
    // model.geometry.points.multiscale.set_not_on_boundary_fun([&](Size p){return model.geometry.not_on_boundary(p);});
    // cout << "here too" << endl;
    // model.geometry.points.multiscale.coarsen();
    // cout << "max points deleted" << endl;
    // auto alwaystrue = [](Size p1, Size p2)
    // {
    //   return true;
    // };
    // model.geometry.points.multiscale.set_comparison_fun(alwaystrue);
    //
    // model.geometry.points.multiscale.coarsen();
    // vector<bool> curr_mask=model.geometry.points.multiscale.mask[2];
    // cout << "nb of points remaining after coarsening: " << std::count (curr_mask.begin(), curr_mask.end(), true) << endl;
    // //TODO: change function name to get_MAX_...
    // cout << model.geometry.points.multiscale.get_curr_coars_lvl() << endl;
    // //because of the way how we coarsen, point 0 should still lie in the grid
    // //TODO add output for debugging....
    // std::set<Size> neighbors_after_del=model.geometry.points.multiscale.get_neighbors(0,2);
    // for (Size nb:neighbors_after_del)
    // {
    // cout << nb << endl;
    // }
    // cout << "now printing neighbors of neighbors; should all contain the original point" << endl;
    // for (Size nb:neighbors_after_del)
    // {
    //   cout << "point: " << nb << endl;
    //   std::set<Size> temp_neighbors=model.geometry.points.multiscale.get_neighbors(nb,2);
    //   for (Size nnb: temp_neighbors)
    //   {
    //     cout << nnb << endl;
    //   }
    // }
    // //also try to coarsen again
    // cout << "Another coarsening" << endl;
    // model.geometry.points.multiscale.coarsen();
    // neighbors_after_del=model.geometry.points.multiscale.get_neighbors(0,3);
    // for (Size nb:neighbors_after_del)
    // {
    // cout << nb << endl;
    // }
    // curr_mask=model.geometry.points.multiscale.mask[3];
    // cout << "nb of points remaining after second coarsening: " << std::count (curr_mask.begin(), curr_mask.end(), true) << endl;
    //
    //
    //
    // //now we try to interpolate a vector of ones
    // std::vector<double> ones(parameters.npoints(), 1);
    // cout<<"now trying to interpolate"<<endl;
    // model.interpolate_vector_local(2,ones);//note: in this test, the first coarsing does nothing; this leads to throwing bad_alloc
    // cout<<"end interpolation"<<endl;
    // for (Size idx=0; idx<parameters.npoints(); idx++)
    // {cout<<ones[idx]<<std::endl;}



    //trying to delete 0.2 percent of the points in the grid
    //cout << "no of points to delete = " << int(sizeof(Points)*0.01) << endl;
    // model.coarsen_grid(0.1);
    // cout << "done with deleting points" << endl;
    // std::vector<double> vector1(68994, 0.0);
    // model.interpolate_vector(1, 0, vector1);




//    Solver solver (10000, 100);
//    solver.trace (model);
//
//    for (Size i = 0; i < 10; i++)
//    {
//        cout << model.geometry.lengths[i] << endl;
//    }

    // Size1 lengths = model.geometry.get_ray_lengths ();
//    Size1 lengths = model.geometry.get_ray_lengths_gpu (512, 512);

    // for (Size i = 0; i < 100; i++)
    // {
//        cout << model.geometry.lengths[i] << endl;
        // cout << lengths[i] << endl;
    // }

    cout << "Done." << endl;

    return (0);
}
