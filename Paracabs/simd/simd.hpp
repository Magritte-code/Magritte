#pragma once

#if PARACABS_USE_SIMD

    // SIMD types
    typedef double vdouble __attribute__ (( vector_size (4*sizeof(double)) ));
    typedef float  vfloat  __attribute__ (( vector_size (4*sizeof(float )) ));

#else

    typedef double vdouble;
    typedef float  vfloat;

#endif