#pragma once


#include "Paracabs/configure.hpp"


#if PARACABS_USE_MESSAGE_PASSING && PARACABS_USE_MPI

    #include "message_passing_mpi.hpp"

#else

    namespace message_passing
    {
        inline void initialize () {}
        inline void   finalize () {}

        inline unsigned int comm_size () {return 1;}
        inline unsigned int comm_rank () {return 0;}

        inline size_t start  (const size_t total) {return 0;    }
        inline size_t stop   (const size_t total) {return total;}
        inline size_t length (const size_t total) {return total;}
    }

#endif


namespace message_passing
{
    #define distributed_for (i, total) for (size_t i = start(total); i < stop(total); i++)
}
