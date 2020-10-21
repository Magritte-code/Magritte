#include <iostream>
using std::cout;
using std::endl;

#include "io/python/io_python.hpp"
#include "model/model.hpp"
#include "solver/solver.hpp"
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
    /// Store model name
    const string modelName = argv[1];

    cout << "Running test_0th_short..." << endl;
    cout << "-------------------------" << endl;
    cout << "Model name: " << modelName << endl;
    cout << "n threads = " << paracabs::multi_threading::n_threads_avail() << endl;


    Model model (modelName);

    model.compute_spectral_discretisation ();
    model.compute_LTE_level_populations   ();
    model.compute_inverse_line_widths     ();

    model.compute_radiation_field_0th_short_characteristics ();
    // solver.trace (model);

    // Solver solver (model.parameters.npoints(),
                   // model.parameters.nfreqs () );
    // solver.solve  (model);

//
//    for (Size i = 0; i < 10; i++)
//    {
//        cout << model.geometry.lengths[i] << endl;
//    }

    // Size1 lengths = model.geometry.get_ray_lengths ();
//    Size1 lengths = model.geometry.get_ray_lengths_gpu (512, 512);

    // cout << "r = 0" << endl;

    // for (Size p = 0; p < model.parameters.npoints(); p++)
    // {
        // model.radiation.print(0, p);
    // }

    // cout << "r = 1" << endl;

    // for (Size p = 0; p < model.parameters.npoints(); p++)
    // {
        // model.radiation.print(1, p);
    // }

    cout << "Done." << endl;

    return (0);
}
