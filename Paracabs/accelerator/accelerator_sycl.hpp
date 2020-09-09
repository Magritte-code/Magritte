#pragma once

#include <iostream>
#include <string>
using std::string;

#include <CL/sycl.hpp>

#define PARACABS_DEBUG false

#if PARACABS_DEBUG
    #define handle_sycl_error(body)   \
    {                                 \
        body                          \
    }
#else
    #define handle_sycl_error(body) body
#endif


#define accel


namespace paracabs
{
    namespace accelerator
    {
        ///  Getter for the number of available GPUs
        ///    @returns number of available GPUs
        ////////////////////////////////////////////
        inline unsigned int nGPUs ()
        {
            int n;
            handle_cuda_error (cudaGetDeviceCount (&n));
            return n;
        }


        ///  Getter for the name of a GPU
        ///    @param[in] i : number of the GPU
        ///    @returns name of GPU with number i
        //////////////////////////////////////////
        inline string get_gpu_name (const int i)
        {
            cudaDeviceProp prop;
            handle_cuda_error (cudaGetDeviceProperties (&prop, i));
            return prop.name;
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
            handle_cuda_error (cudaDeviceSynchronize ());
        }


        ///  Allocate memory on the device
        ///    @param[in] num : number of bytes to allocate
        ///    @returns pointer to the allocated block of memory
        ////////////////////////////////////////////////////////
        inline void* malloc (const size_t num)
        {
            void* ptr;
            handle_cuda_error (cudaMalloc (&ptr, num));
            handle_cuda_error (cudaDeviceSynchronize ());
            return ptr;
        }


        ///  Free memory on the device
        ///    @param[in] ptr : pointer to the block to free
        ////////////////////////////////////////////////////
        inline void free (void* ptr)
        {
            handle_cuda_error (cudaDeviceSynchronize ());
            handle_cuda_error (cudaFree (ptr));
        }


//        inline void memcpy (void* dst, const void* src, const size_t size)
//        {
//            handle_cuda_error (cudaMemcpy (dst, src, size, cudaMemcpyDeviceToDevice));
//        }


        inline void memcpy_to_accelerator (void* dst, const void* src, const size_t size)
        {
            handle_cuda_error (cudaMemcpy (dst, src, size, cudaMemcpyHostToDevice));
        }


        inline void memcpy_from_accelerator (void* dst, const void* src, const size_t size)
        {
            handle_cuda_error (cudaMemcpy (dst, src, size, cudaMemcpyDeviceToHost));
        }
    }
}


/// Assumes to be used within a class!

#define accelerated_for(i, total, nblocks, nthreads, ... )   \
{							                                 \
    copyContextAccelerator() = true;                         \
    auto lambda = [=, *this] __device__ (size_t i) mutable   \
    {		                                                 \
        __VA_ARGS__;							             \
    };									                     \
    apply_lambda <<<nblocks, nthreads>>> (total, lambda);    \
    copyContextAccelerator() = false;                        \
}


#define accelerated_for_outside_class(i, total, nblocks, nthreads, ... )   \
{							                                               \
    copyContextAccelerator() = true;                                       \
    auto lambda = [=] __device__ (size_t i) mutable                        \
    {		                                                               \
        __VA_ARGS__;							                           \
    };									                                   \
    apply_lambda <<<nblocks, nthreads>>> (total, lambda);                  \
    copyContextAccelerator() = false;                                      \
}


template <typename Lambda>
__global__ void apply_lambda (const size_t stop, Lambda lambda)
{
    const size_t start  = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t stride =  gridDim.x * blockDim.x;

    for (size_t i = start; i < stop; i += stride)
    {
        lambda (i);
    }
}
