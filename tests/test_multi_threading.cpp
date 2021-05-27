#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "gtest/gtest.h"
#include "paracabs.hpp"
using namespace paracabs;


/// Test the setter for the available threads
/////////////////////////////////////////////
TEST (multi_threading, set_n_threads_avail)
{
    multi_threading::set_n_threads_avail (1);
    EXPECT_EQ (1, multi_threading::n_threads_avail());

    multi_threading::set_n_threads_avail (2);
    EXPECT_EQ (2, multi_threading::n_threads_avail());

    multi_threading::set_n_threads_avail (4);
    EXPECT_EQ (4, multi_threading::n_threads_avail());
}


/// Test the start and stop values
//////////////////////////////////
TEST (multi_threading, start_stop)
{
    const size_t n_threads =    4;
    const size_t n         = 1234;

    multi_threading::set_n_threads_avail (n_threads);

    vector<size_t> starts (n_threads);
    vector<size_t> stops  (n_threads);

    PRAGMA_PARALLEL
    {
        starts[multi_threading::thread_id()] = multi_threading::start(n);
        stops [multi_threading::thread_id()] = multi_threading::stop (n);
    }

    EXPECT_EQ (0, starts[          0]);
    EXPECT_EQ (n, stops [n_threads-1]);

    for (size_t t = 1; t < n_threads; t++)
    {
        EXPECT_EQ (starts[t], stops[t-1]);
    }
}


/// Test the threaded_for loop
//////////////////////////////
TEST (multi_threading, threaded_for)
{
    const size_t n_threads =    4;
    const size_t n         = 1000;

    vector<double> a (n, 0.5);
    vector<double> b (n, 1.5);
    vector<double> c (n, 0.0);

    multi_threading::set_n_threads_avail (n_threads);

    threaded_for (i, n,
    {
        c[i] = a[i] + b[i];
    })

    double total = 0.0;
    for (const double x : c) {total += x;}

    EXPECT_EQ (total, 2.0*n);
}


/// Test ThreadPrivate variable
///////////////////////////////
TEST (multi_threading, ThreadPrivate)
{
    const size_t n_threads =    4;
    const size_t n         = 1000;

    multi_threading::ThreadPrivate<vector<double>> a = vector<double> (n);

    multi_threading::set_n_threads_avail (n_threads);

    PRAGMA_PARALLEL
    {
        for (size_t i = 0; i < n; i++)
        {
            a()[i] = multi_threading::thread_id();
        }
    }

    for (size_t t = 0; t < n_threads; t++)
    {
        for (size_t i = 0; i < n; i++)
        {
            EXPECT_EQ (t, a(t)[i]);
        }
    }
}


/// Main function for GoogleTest
////////////////////////////////
int main (int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
