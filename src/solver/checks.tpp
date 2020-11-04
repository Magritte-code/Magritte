#include <algorithm>


inline void Solver :: setup_T ()
{
    const Size N = n_tot[w/nfreqs_red];

    test_T = MatrixXd::Zero (N, N);

    test_T(0,0) =  1.0 + 2.0 / dtau[I(0,w)] + 2.0 / (dtau[I(0,w)]*dtau[I(0,w)]) ;
    test_T(0,1) = -C[I(0,w)];

    for (Size n = 1; n < N-1; n++)
    {
        test_T(n,n-1) = -A[I(n,w)];
        test_T(n,n  ) = 1.0 + A[I(n,w)] + C[I(n,w)];
        test_T(n,n+1) = -C[I(n,w)];
    }

    test_T(N-1,N-2) = -A[I(N-1,w)];
    test_T(N-1,N-1) =  1.0 + 2.0 / dtau[I(N-2,w)] + 2.0 / (dtau[I(N-2,w)]*dtau[I(N-2,w)]) ;
}




template <typename Real, typename DataLayout>
inline void Solver<Real, DataLayout> :: setup_L (const Size w)
{
    const Size N = n_tot[w/nfreqs_red];

    test_L = MatrixXd::Zero (N, N);

    for (Size n = 0; n < N; n++)
    {
        test_L(n,n) = L_diag[I(n,w)];
    }

    for (Size m = 0; (m < n_off_diag) && (m < N-1); m++)
    {
        for (Size n = 0; n < N-m-1; n++)
        {
            test_L(n,n+m+1) = L_upper[M(m,I(n+m+1,w))];
            test_L(n+m+1,n) = L_lower[M(m,I(n    ,w))];
        }
    }
}




template <typename Real>
inline Real abs_rel_diff (const MatrixXd A, const MatrixXd B, const Size i, const Size j)
{
//    if (isnan(A(i,j)) || isnan(B(i,j)))
//    {
//        return 1.0e+99;
//    }

    return fabs (2.0 * (A(i,j) - B(i,j)) / (A(i,j) + B(i,j)));
}




template <typename Real, typename DataLayout>
inline Real Solver<Real, DataLayout> :: check_L (const Size w)
{
    const Size N = n_tot[w/nfreqs_red];

    setup_T (w);
    setup_L (w);

    // Define the inverse of T
    const MatrixXd test_T_inverse = test_T.inverse();

    // Check precision
    Real prec = 0.0;

    for (Size n = 0; n < N; n++)
    {
        prec = std::max(prec, abs_rel_diff<Real>(test_L, test_T_inverse, n, n));
    }

    for (Size m = 0; (m < n_off_diag) && (m < N-1); m++)
    {
        for (Size n = 0; n < N-m-1; n++)
        {
            prec = std::max(prec, abs_rel_diff<Real>(test_L, test_T_inverse, n, n+m+1));
            prec = std::max(prec, abs_rel_diff<Real>(test_L, test_T_inverse, n+m+1, n));
        }
    }

    if (prec > 1.0e-3)
    {
        cout         << "- test: n_off_diag = " << n_off_diag << endl;
        cout << endl << "- test: T ="      << endl << test_T         << endl;
        cout << endl << "- test: T^(-1) =" << endl << test_T_inverse << endl;
        cout << endl << "- test: L ="      << endl << test_L         << endl;
        cout << endl << "- test: precision = " << prec               << endl;
    }

    return prec;
}
