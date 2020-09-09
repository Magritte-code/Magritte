#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "paracabs.hpp"
using namespace paracabs;


int main ()
{
    cout << "Paracabs test multi threading." << endl;

    const unsigned int n_threads       = multi_threading::n_threads();
    const unsigned int thread_id       = multi_threading::thread_id();
    const unsigned int n_threads_avail = multi_threading::n_threads_avail();

    cout << "n_threads        = " << n_threads       << endl;
    cout << "thread_id        = " << thread_id       << endl;
    cout << "n_threads_avail  = " << n_threads_avail << endl;

    const size_t total = 10;

    const size_t start = multi_threading::start(total);
    const size_t stop  = multi_threading::stop (total);

    cout << "total = " << total << endl;
    cout << "start = " << start << endl;
    cout << "stop  = " << stop  << endl;

    vector <double> vec (total, 1.0);

    threaded_for (i, total,{
        vec[i] += 1.0;
    })

    PRAGMA_PARALLEL
    {
        cout << "thread_id = " << multi_threading::thread_id () << ", "
             <<     "start = " << multi_threading::start(total) << ", "
             <<      "stop = " << multi_threading::stop (total) << endl;
    }

    cout << "Done." << endl;

    return (0);
}