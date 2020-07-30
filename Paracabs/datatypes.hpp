#pragma once


#include "configure.hpp"
#include "simd/simd.hpp"
#include "accelerator/accelerator.hpp"
#include "message_passing/message_passing.hpp"
#include "multi_threading/multi_threading.hpp"


// Array types
template <typename type> struct Array1d
{
    type *data;

    const size_t size;

    Array1d (const size_t s) : size (s)
    {
        data = new type[size];
    }

    ~Array1d ()
    {
        delete[] data;
    }
};