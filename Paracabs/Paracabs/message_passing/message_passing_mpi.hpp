#pragma once


#include <mpi.h>


namespace paracabs
{
    namespace message_passing
    {
        inline void initialize ()
        {
            MPI_Init (NULL, NULL);
        }


        inline void finalize ()
        {
            MPI_Finalize ();
        }


        inline unsigned int comm_size ()
        {
            int size;
            MPI_Comm_size (MPI_COMM_WORLD, &size);

            return size;
        }


        inline unsigned int comm_rank ()
        {
            int rank;
            MPI_Comm_rank (MPI_COMM_WORLD, &rank);

            return rank;
        }


        inline size_t start (const size_t total)
        {
            return (comm_rank() * total) / comm_size();
        }


        inline size_t stop (const size_t total)
        {
            return ((comm_rank()+1) * total) / comm_size();
        }


        inline size_t length (const size_t total)
        {
            return stop(total) - start(total);
        }
    }
}
