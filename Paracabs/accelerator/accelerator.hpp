#pragma once


#include "configure.hpp"


#if PARACABS_USE_ACCELERATOR && PARACABS_USE_CUDA

    #include "accelerator_cuda.cuh"

#else

    namespace accelerator
    {
    }

#endif
