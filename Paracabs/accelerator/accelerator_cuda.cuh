#pragma once


//#ifdef __CUDACC__


#include "cuda_runtime.h"
#include "device_launch_parameters.h"


namespace accelerator
{

    struct MemoryManagement
    {
        inline void* allocate_memory (const size_t size)
        {
            void* ptr;
            cudaMalloc (&ptr, size);
            return ptr;
        }


        inline void free_memory (void* ptr)
        {
            cudaFree (ptr);
        }
    };

}

//#endif