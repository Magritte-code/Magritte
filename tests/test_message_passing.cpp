#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;

#include "Paracabs/message_passing/message_passing.hpp"


int main ()
{
    cout << "Paracabs test message passing." << endl;

    message_passing::initialize();


    const unsigned int size = message_passing::comm_size();
    const unsigned int rank = message_passing::comm_rank();

    cout << "size  = " << size  << endl;
    cout << "rank  = " << rank  << endl;


    const size_t total = 10;

    const size_t start = message_passing::start(total);
    const size_t stop  = message_passing::stop (total);

    cout << "total = " << total << endl;
    cout << "start = " << start << endl;
    cout << "stop  = " << stop  << endl;

    vector <double> vec (total, 1.0);

    distributed_for (i, i_loc, total,
    {
        vec[i] += 1.0;
    })

    if (message_passing::comm_rank() == 0)
    {
//        for (size_t w = 0; w < message_passing::comm_size(); w++)
        {
            for (size_t i = 0; i < total; i++)
            {
                cout << "vec[" << i << "] = " << vec[i] << endl;
            }
        }
    }

    if (message_passing::comm_rank() == 1)
    {
//        for (size_t w = 0; w < message_passing::comm_size(); w++)
        {
            for (size_t i = 0; i < total; i++)
            {
                cout << "             ";
                cout << "vec[" << i << "] = " << vec[i] << endl;
            }
        }
    }

    message_passing::finalize();


    cout << "Done." << endl;

    return (0);
}