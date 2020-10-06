#pragma once

#include <assert.h>
#include "tools/types.hpp"


template <typename type>
inline Size max (const vector<type>& a, const Size n, const Size i, const Size j, const Size k)
{
    Size m = i;

    if ( (j < n) && (a[j] > a[m]) ) {m = j;}
    if ( (k < n) && (a[k] > a[m]) ) {m = k;}

    return m;
}


template <typename type>
inline void downheap (vector<type>& a, Size1& b, Size n, Size i)
{
    while (1)
    {
        Size j = max (a, n, i, 2*i+1, 2*i+2);

        if (j == i) {break;}

        double temp1 = a[i];
        Size   temp2 = b[i];

        a[i] = a[j];
        a[j] = temp1;

        b[i] = b[j];
        b[j] = temp2;

        i = j;
    }
}


template <typename type>
inline void heapsort (vector<type>& a, Size1& b)
{
    // Get vector length
    const Size n = a.size();

    // Assert that both vectors have the same size
    assert (n == b.size());

    for (long i = (n-2)/2; i >=0 ; i--) // Warning: long is required here, couting down!!!
    {
        downheap (a, b, n, i);
    }

    for (Size i = 0; i < n; i++)
    {
        double temp1 = a[n-i-1];
        Size   temp2 = b[n-i-1];

        a[n-i-1] = a[0];
        a[0]     = temp1;

        b[n-i-1] = b[0];
        b[0]     = temp2;

        downheap (a, b, n-i-1, 0);
    }
}
