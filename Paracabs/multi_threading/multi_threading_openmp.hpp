#pragma once


#include <omp.h>


namespace multi_threading
{
    inline unsigned int n_threads ()
    {
        return omp_get_num_threads();
    }


    inline unsigned int thread_id ()
    {
        return omp_get_thread_num();
    }


    inline size_t start (const size_t total)
    {
        return (thread_id() * total) / n_threads();
    }


    inline size_t stop (const size_t total)
    {
        return ((thread_id()+1) * total) / n_threads();
    }


    inline unsigned int n_threads_avail ()
    {
        int n;

        #pragma omp parallel
        {
            n = n_threads();
        }

        return n;
    }
}
