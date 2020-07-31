#pragma once


#include "configure.hpp"
#include "datatypes_vector3d.hpp"


//enum MemoryMode {shared, distributed, accelerator};





namespace datatypes
{
    struct MemoryManagement
    {
        inline void* allocate_memory (const size_t size) const
        {
            return malloc (size);
        }


        inline void free_memory (void* ptr) const
        {
            free (ptr);
        }
    };


    template <typename type, class MemoryManagement>
    struct array1d : private MemoryManagement
    {
        type* data;

        const size_t size;

        array1d (const size_t s) : size (s)
        {
            data = (type*) this->allocate_memory (size*sizeof(type));
        }

        ~array1d ()
        {
            this->free_memory (data);
        }
    };


}