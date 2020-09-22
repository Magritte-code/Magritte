#pragma once


#include "configure.hpp"


#if PARACABS_USE_MULTI_THREADING && PARACABS_USE_OPENMP

    #include "multi_threading_openmp.hpp"

    #define PRAGMA_PARALLEL     _Pragma("omp parallel     default (shared)")
    #define PRAGMA_PARALLEL_FOR _Pragma("omp parallel for default (shared)")

#else

    namespace paracabs
    {
        namespace multi_threading
        {
            inline unsigned int n_threads () {return 1;}
            inline unsigned int thread_id () {return 0;}

            inline size_t start (const size_t total) {return 0;    }
            inline size_t stop  (const size_t total) {return total;}

            inline unsigned int n_threads_avail () {return 1;}
        }
    }

    #define PRAGMA_PARALLEL
    #define PRAGMA_PARALLEL_FOR

#endif


#define threaded_for(i, total, ...)                             \
    PRAGMA_PARALLEL                                             \
    for (size_t i = paracabs::multi_threading::start(total);    \
                i < paracabs::multi_threading::stop (total);    \
                i++                                         )   \
    {                                                           \
        __VA_ARGS__;                                            \
    }