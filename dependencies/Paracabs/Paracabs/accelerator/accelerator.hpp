#pragma once


#include "configure.hpp"


#if PARACABS_USE_ACCELERATOR && PARACABS_USE_CUDA

    #include "accelerator_cuda.hpp"

#else

    #include "multi_threading/multi_threading.hpp"

    namespace paracabs
    {
        namespace accelerator
        {
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
                free (ptr);
            }


//        inline void memcpy (void* dst, const void* src, const size_t size)
//        {
//            handle_cuda_error (cudaMemcpy (dst, src, size, cudaMemcpyDeviceToDevice));
//        }


            inline void memcpy_to_accelerator (void* dst, const void* src, const size_t size)
            {
                std::memcpy (dst, src, size);
            }


            inline void memcpy_from_accelerator (void* dst, const void* src, const size_t size)
            {
                std::memcpy (dst, src, size);
            }
        }
    }


    #define accel

    #define accelerated_for(i, total, nblocks, nthreads, ... ) \
        threaded_for(i, total, __VA_ARGS__)

#endif

