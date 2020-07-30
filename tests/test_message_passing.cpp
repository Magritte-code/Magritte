#include <iostream>
using std::cout;
using std::endl;

#include "Paracabs/message_passing/message_passing.hpp"


int main ()
{
    cout << "Paracabs test message passing." << endl;

    message_passing::initialize();

    const unsigned int size = message_passing::comm_size();
    const unsigned int rank = message_passing::comm_rank();

    const size_t total = 1000000;

    const size_t start = message_passing::start(total);
    const size_t stop  = message_passing::stop (total);

    message_passing::finalize();

    cout << "size  = " << size  << endl;
    cout << "rank  = " << rank  << endl;
    cout << "start = " << start << endl;
    cout << "stop  = " << stop  << endl;

    cout << "Done." << endl;

    return (0);
}