#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "gtest/gtest.h"
#include "paracabs.hpp"
using namespace paracabs;


/// Test the accelerator list
/////////////////////////////
TEST (acceleratpr, accelerator_list)
{
    accelerator::list_accelerators();
}



/// Test the accelerated_for loop
//////////////////////////////
TEST (multi_threading, threaded_for)
{
    const size_t n = 100;

    accelerator::nblocks  = 1;
    accelerator::nthreads = 1;

    datatypes::Vector<double> a (n, 0.5);
    datatypes::Vector<double> b (n, 1.5);
    datatypes::Vector<double> c (n, 0.0);

    accelerator::list_accelerators();

    accelerated_for_outside_class (i, n,
    {
        c[i] = a[i] + b[i];
    })

    accelerator::synchronize();

    c.copy_vec_to_ptr();

    for (size_t i = 0; i < n; i++) {cout << c[i] << endl;}
}


/// Main function for GoogleTest
////////////////////////////////
int main (int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
