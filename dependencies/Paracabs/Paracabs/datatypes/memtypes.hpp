#pragma once


#include <cstring> // std::memcpy
#include <iostream>
using std::cout;
using std::endl;

#include "accelerator/accelerator.hpp"


namespace paracabs
{
    namespace datatypes
    {
        struct MemTypeDefault
        {
            ///  Memory allocation
            ///    @param[in] num : nuber of bytes to allocate
            ///    @returns pointer to allocated block of memory
            ////////////////////////////////////////////////////
            inline void* malloc (const size_t num)
            {
                cout << "--- Called MemTypeDefault malloc" << endl;
                return std::malloc (num);
            }

            ///  Memory deallocation
            ///    @param[in] ptr : pointer to memory block to free
            ///////////////////////////////////////////////////////
            inline void free (void* ptr)
            {
                cout << "--- Called MemTypeDefault free" << endl;
                std::free (ptr);
            }

//            inline void memcpy (void* dst, const void* src, const size_t num)
//            {
//                std::memcpy (dst, src, num);
//            }

            /// Overloaded new operator
            ///    @param[in] num : nuber of bytes to allocate
            ///    @returns pointer to allocated block of memory
            ////////////////////////////////////////////////////
            void* operator new (size_t num)
            {
                cout << "--- Called MemTypeDefault new" << endl;
                return std::malloc (num);
            }

            /// Overloaded delete operator
            ///    @param[in] ptr : pointer to memory block to free
            ///////////////////////////////////////////////////////
            void operator delete (void* ptr)
            {
                cout << "--- Called MemTypeDefault delete" << endl;
                std::free (ptr);
            }
        };


        #if PARACABS_USE_ACCELERATOR

        struct MemTypeAccelerator
        {
            ///  Memory allocation
            ///    @param[in] num : nuber of bytes to allocate
            ///    @returns pointer to allocated block of memory
            ////////////////////////////////////////////////////
            inline void* malloc (const size_t num)
            {
                cout << "--- Called MemTypeAccelerator malloc" << endl;
                return paracabs::accelerator::malloc (num);
            }

            ///  Memory deallocation
            ///    @param[in] ptr : pointer to memory block to free
            ///////////////////////////////////////////////////////
            inline void free (void* ptr)
            {
                cout << "--- Called MemTypeAccelerator free" << endl;
                paracabs::accelerator::free (ptr);
            }

            /// Overloaded new operator
            ///    @param[in] num : nuber of bytes to allocate
            ///    @returns pointer to allocated block of memory
            ////////////////////////////////////////////////////
            void* operator new (size_t num)
            {
                cout << "--- Called MemTypeAccelerator new" << endl;
                return paracabs::accelerator::malloc (num);
            }

            /// Overloaded delete operator
            ///    @param[in] ptr : pointer to memory block to free
            ///////////////////////////////////////////////////////
            void operator delete (void* ptr)
            {
                cout << "--- Called MemTypeAccelerator delete" << endl;
                paracabs::accelerator::free (ptr);
            }

//            inline void memcpy (void* dst, const void* src, const size_t num)
//            {
//                paracabs::accelerator::memcpy (dst, src, num);
//            }

            inline void memcpy_to_accelerator (void* dst, const void* src, const size_t num)
            {
                paracabs::accelerator::memcpy_to_accelerator (dst, src, num);
            }

            inline void memcpy_from_accelerator (void* dst, const void* src, const size_t num)
            {
                paracabs::accelerator::memcpy_from_accelerator (dst, src, num);
            }
        };

        #else

//        using MemtypeAccelerator = MemTypeDefault;
        typedef MemTypeDefault MemTypeAccelerator;

        #endif
    }
}