#pragma once

inline bool& copyContextAccelerator()
{
    static bool value = false;
    return value;
}

#include "configure.hpp"
#include "simd/simd.hpp"
#include "datatypes/datatypes.hpp"
#include "accelerator/accelerator.hpp"
#include "message_passing/message_passing.hpp"
#include "multi_threading/multi_threading.hpp"
