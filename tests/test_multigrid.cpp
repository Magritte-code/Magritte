#include <iostream>
using std::cout;
using std::endl;

#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"
#include <set>


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

    cout << "n threads = " << paracabs::multi_threading::n_threads_avail() << endl;


    Model model;
    model.read (io);
    vector<Size> old_neighbors=model.geometry.points.curr_neighbors.get_neighbors(9001);
    std::set<Size> neighbors=model.geometry.points.multiscale.get_neighbors(9001,0);
    for (Size nb:neighbors)
    {
    cout << nb << endl;
    }
    cout << "compare with" << endl;
    for (Size nb:old_neighbors)
    {
    cout << nb << endl;
    }

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

//     Size1 lengths = model.geometry.get_ray_lengths ();
// //    Size1 lengths = model.geometry.get_ray_lengths_gpu (512, 512);
//
//     for (Size i = 0; i < 100; i++)
//     {
// //        cout << model.geometry.lengths[i] << endl;
//         cout << lengths[i] << endl;
//     }
//
//     cout << "Done." << endl;
//
//     return (0);
}
