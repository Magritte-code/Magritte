#pragma once

#include <iostream>
#include <string>
using std::string;

#include <CL/sycl.hpp>
namespace sycl = cl::sycl;

#define PARACABS_DEBUG true

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

        extern sycl::queue acceleratorQueue;

        ///  Getter for the number of available GPUs
        ///    @returns number of available GPUs
        ////////////////////////////////////////////
//        inline unsigned int nGPUs ()
//        {
//            int n;
//            handle_cuda_error (cudaGetDeviceCount (&n));
//            return n;
//        }


        ///  Getter for the name of a GPU
        ///    @param[in] i : number of the GPU
        ///    @returns name of GPU with number i
        //////////////////////////////////////////
//        inline string get_gpu_name (const int i)
//        {
//            cudaDeviceProp prop;
//            handle_cuda_error (cudaGetDeviceProperties (&prop, i));
//            return prop.name;
//        }


        ///  Lists the available accelerators
        /////////////////////////////////////
        inline void list_accelerators ()
        {
	    // getting the list of all supported sycl platforms
	    auto platform_list = sycl::platform::get_platforms();
	    // looping over platforms
	    for (const auto &platform : platform_list)
	    {
	        // getting the list of devices from the platform
	        auto device_list = platform.get_devices();
	        // looping over devices
	        for (const auto &device : device_list)
		{
		    cout << device.get_info <sycl::info::device::name>() << endl;
	        }
	    }
        }


        ///  Barrier for accelerator threads
        ////////////////////////////////////
        inline void synchronize ()
        {
            acceleratorQueue.wait();
        }


        ///  Allocate memory on the device
        ///    @param[in] num : number of bytes to allocate
        ///    @returns pointer to the allocated block of memory
        ////////////////////////////////////////////////////////
        inline void* malloc (const size_t num)
        {
            return sycl::malloc_device (num, acceleratorQueue);
        }


        ///  Free memory on the device
        ///    @param[in] ptr : pointer to the block to free
        ////////////////////////////////////////////////////
        inline void free (void* ptr)
        {
            acceleratorQueue.wait();
            sycl::free (ptr, acceleratorQueue);
        }


        inline void memcpy (void* dst, const void* src, const size_t size)
        {
//            handle_cuda_error (cudaMemcpy (dst, src, size, cudaMemcpyDeviceToDevice));
        }


        inline void memcpy_to_accelerator (void* dst, const void* src, const size_t size)
        {
            acceleratorQueue.memcpy (dst, src, size);
            acceleratorQueue.wait();
        }


        inline void memcpy_from_accelerator (void* dst, const void* src, const size_t size)
        {
            acceleratorQueue.memcpy (dst, src, size);
            acceleratorQueue.wait();
        }
    }
}


/// Assumes to be used within a class!

#define accelerated_for(i, total, nblocks, nthreads, ... )   \
{                                                            \
    copyContextAccelerator() = true;                         \
    paracabs::accelerator::acceleratorQueue.submit (         \
        [&](sycl::handler &cgh)                              \
        {                                                    \
            cgh.parallel_for<class kernel> (                 \
                sycl::range<1> {total},                      \
                [=] (sycl::id<1> i) mutable                  \
                {                                            \
                    __VA_ARGS__                              \
                }                                            \
            );                                               \
        }                                                    \
    );                                                       \
    copyContextAccelerator() = false;                        \
}


#define accelerated_for_outside_class(i, total, nblocks, nthreads, ... )   \
//{                                                                          \
//    copyContextAccelerator() = true;                                       \
//    auto lambda = [=] __device__ (size_t i) mutable                        \
//    {                                                                      \
//        __VA_ARGS__;                                                       \
//    };                                                                     \
//    apply_lambda <<<nblocks, nthreads>>> (total, lambda);                  \
//    copyContextAccelerator() = false;                                      \
//}

