#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <cmath>
using std::sqrt;
#include <Eigen/Dense>

#include "gtest/gtest.h"
#include "model/model.hpp"
#include "solver/solver.hpp"


MatrixXr setup_T (Solver& solver)
{
    const Size N = solver.n_tot_();
    cout << "N = " << N << endl;

    Vector<Real>& A = solver.A_();
    Vector<Real>& C = solver.C_();

    cout << "solver.A_()[10] = " << solver.A_()[10] << endl;
    cout << "       A   [10] = " << A[10] << endl;
    cout << "solver.C_()[10] = " << solver.C_()[10] << endl;
    cout << "       C   [10] = " << C[10] << endl;

    const Size first = solver.first_();
    const Size last  = solver.last_ ();
    cout << "first = " << first << endl;
    cout << "last  = " << last  << endl;

    MatrixXr T = MatrixXr::Zero(N, N);

    T(0,0) = 1.0 + C[first] + 2.0*sqrt(0.5*C[first]);
    T(0,1) = -C[first];

    for (Size n = 1; n < N-1; n++)
    {
        cout << "A[first+n] = " << A[first+n] << endl;
        T(n,n-1) = -A[first+n];
        T(n,n  ) = 1.0 + A[first+n] + C[first+n];
        T(n,n+1) = -C[first+n];
    }

    T(N-1,N-2) = -A[last];
    T(N-1,N-1) = 1.0 + A[last] + 2.0*sqrt(0.5*A[last]);

    cout << "T = " << endl;
    cout << T      << endl;

    return T;
}


MatrixXr setup_L (Solver& solver)
{
    const Size N = solver.n_tot_();

    Vector<Real>& L_diag  = solver.L_diag_ ();
    Matrix<Real>& L_upper = solver.L_upper_();
    Matrix<Real>& L_lower = solver.L_lower_();

    const Size first = solver.first_();
    const Size last  = solver.last_ ();

    MatrixXr L = MatrixXr::Zero (N, N);

    for (Size n = 0; n < N; n++)
    {

        L(n,n) = L_diag[first+n];
    }

    for (Size m = 0; (m < solver.n_off_diag) && (m < N-1); m++)
    {
        for (Size n = 0; n < N-m-1; n++)
        {
            L(n,n+m+1) = L_upper(m,first+n+m+1);
            L(n+m+1,n) = L_lower(m,first+n    );
        }
    }

    cout << "L = " << endl;
    cout << L      << endl;

    return L;
}


// MatrixXr extract_L (Model& model)
// {
//     // Lambda lambda = model.lines.lineProducingSpecies[0].lambda;
//
//
//
//
//     const Size N = solver.n_tot_();
//
//     Vector<Real>& L_diag  = solver.L_diag_ ();
//     Matrix<Real>& L_upper = solver.L_upper_();
//     Matrix<Real>& L_lower = solver.L_lower_();
//
//     const Size first = solver.first_();
//     const Size last  = solver.last_ ();
//
//     MatrixXr L = MatrixXr::Zero (N, N);
//
//     for (Size n = 0; n < N; n++)
//     {
//
//         L(n,n) = L_diag[first+n];
//     }
//
//     for (Size m = 0; (m < solver.n_off_diag) && (m < N-1); m++)
//     {
//         for (Size n = 0; n < N-m-1; n++)
//         {
//             L(n,n+m+1) = L_upper(m,first+n+m+1);
//             L(n+m+1,n) = L_lower(m,first+n    );
//         }
//     }
//
//     cout << "L = " << endl;
//     cout << L      << endl;
//
//     return L;
// }




TEST (solver_lambda, lambda)
{
    const string modelFile = magritte_folder + "/tests/models/density_distribution_VZa_1D.hdf5";

    Model model = Model (modelFile);

    const Size length_max = 4*model.parameters->npoints() + 1;
    const Size  width_max =   model.parameters->nfreqs ();

    model.parameters->n_off_diag = model.parameters->npoints();

    Solver solver;
    solver.setup <CoMoving>               (model);
    solver.solve_feautrier_order_2 <None> (model);

    MatrixXr T = setup_T (solver);
    MatrixXr L = setup_L (solver);

    MatrixXr T_inverse = T.inverse();

    cout << "det(T) = " << T.determinant() << endl;

    cout << "L = " << endl;
    cout <<  L     << endl;

    cout << "T_inverse = " << endl;
    cout <<  T_inverse     << endl;


}


int main (int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
