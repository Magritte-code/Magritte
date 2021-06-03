#pragma once

#include <iostream>
#include <cstring>
#include <string>
using std::string;

#include "configure.hpp"


#if   PARACABS_USE_ACCELERATOR && PARACABS_USE_CUDA

    #include "accelerator_cuda.hpp"

#elif PARACABS_USE_ACCELERATOR && PARACABS_USE_SYCL

    #include "accelerator_sycl.hpp"

#else

    #define accel

    #include "multi_threading/multi_threading.hpp"

    namespace paracabs
    {
        namespace accelerator
        {
            ///  Meyers' singleton for nblocks
            //////////////////////////////////
            inline size_t& nblocks()
            {
                static size_t value = 1;
                return value;
            }

            ///  Meyers' singleton for nthreads
            ///////////////////////////////////
            inline size_t& nthreads()
            {
                static size_t value = 1;
                return value;
            }

            ///  Getter for the number of available GPUs
            ///    @returns number of available GPUs
            ////////////////////////////////////////////
            inline unsigned int nGPUs ()
            {
                return 0;
            }

            ///  Getter for the name of a GPU
            ///    @param[in] i : number of the GPU
            ///    @returns name of GPU with number i
            //////////////////////////////////////////
            inline string get_gpu_name (const int i)
            {
                return "";
            }

            ///  Lists the available accelerators
            /////////////////////////////////////
            inline void list_accelerators ()
            {
                for (unsigned int i = 0; i < nGPUs (); i++)
                {
                    std::cout << get_gpu_name (i) << std::endl;
                }
            }

            ///  Barrier for accelerator threads
            ////////////////////////////////////
            inline void synchronize ()
            {
                return;
            }

            ///  Allocate memory on the device
            ///    @param[in] num : number of bytes to allocate
            ///    @returns pointer to the allocated block of memory
            ////////////////////////////////////////////////////////
            inline void* malloc (const size_t num)
            {
                return std::malloc (num);
            }

            ///  Free memory on the device
            ///    @param[in] ptr : pointer to the block to free
            ////////////////////////////////////////////////////
            inline void free (void* ptr)
            {
                std::free (ptr);
            }

            ///  Copy memory from the host to the accelerator
            ///    @param[in] dst  : pointer to destination
            ///    @param[in] src  : pointer to source
            ///    @param[in] size : size of the memory block
            /////////////////////////////////////////////////
            inline void memcpy_to_accelerator (void* dst, const void* src, const size_t size)
            {
                std::memcpy (dst, src, size);
            }

            ///  Copy memory from the host to the accelerator
            ///    @param[in] dst  : pointer to destination
            ///    @param[in] src  : pointer to source
            ///    @param[in] size : size of the memory block
            /////////////////////////////////////////////////
            inline void memcpy_from_accelerator (void* dst, const void* src, const size_t size)
            {
                std::memcpy (dst, src, size);
            }
            
            /// Accelerator thread functionality: fall back to host
            ///////////////////////////////////////////////////////
            using AcceleratorThreads = paracabs::multi_threading::HostThreads;
        }
    }


#define accelerated_for(i, total, ... )    \
{                                          \
    threaded_for(i, total, __VA_ARGS__);   \
}


#define accelerated_for_outside_class(i, total, ... )   \
{                                                       \
    threaded_for(i, total, __VA_ARGS__);                \
}


#endif

