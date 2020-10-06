#pragma once


#include "configure.hpp"
#include "memtypes.hpp"
#include "vector3d.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "array.hpp"


namespace paracabs
{
    namespace datatypes
    {
//    template <typename type>
//    inline void array1d <type, MemTypeDefault> :: copy_to (const array1d <type, MemTypeAccelerator>& array)
//    {
//        assert (size == array.size);
//        array.copy_to_accelerator (array.data, data, size*sizeof(type));
//    }


//        template <typename type>
//        inline void array1d<type, MemTypeDefault> :: copy_to (array1d <type, MemTypeAccelerator>& arr)
//        {
//            arr.memcpy_to_accelerator (arr.data, data, size*sizeof(type));
//        }


//        ///  Copy data from array 1 to array 2
//        ///    @param[in]     arr1 : array1 (MemTypeDefault)
//        ///    @param[in,out] arr2 : array2 (MemTypeAccelerator)
//        ////////////////////////////////////////////////////////
//        template <typename type1, typename type2>
//        inline void my_copy (const Array <type1, MemTypeDefault>*     arr1,
//                                   Array <type2, MemTypeAccelerator>* arr2 )
//        {
//            arr2->memcpy_to_accelerator (arr2->data, arr1->data, arr1->size*sizeof(type1));
//        }
//
//        template <typename type>
//        inline void my_copy (const vector<type, paracabs::allocator<type, MemTypeDefault>>     vec1,
//                                   vector<type, paracabs::allocator<type, MemTypeAccelerator>> vec2 )
//        {
//            assert (vec1.size() == vec2.size());
//
//            for (size_t i = 0; i < vec1.size(); i++)
//            {
//                paracabs::accelerator::memcpy_to_accelerator (&vec2[i], &vec1[i], sizeof(type));
//            }
//        }

//        template <typename type>
//        inline void my_copy (const array1d <type, MemTypeAccelerator>& arr1,
//                                   array1d <type, MemTypeDefault>&     arr2 )
//        {
//            arr1.memcpy_to_accelerator (arr2.data, arr1.data, arr1.size*sizeof(type));
//        }


//    template <typename type, class MemType>
//    struct varray1d : public MemType
//    {
//        type*   data;    ///< variable array data
//        size_t  size;    ///< variable array size
//        size_t* sizes;   ///< sizes of constituent objects
//
//    };
    }
}
