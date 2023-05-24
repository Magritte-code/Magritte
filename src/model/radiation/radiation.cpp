#include "radiation.hpp"

#include "tools/constants.hpp"

#include <iomanip>
#include <iostream>

const string prefix = "radiation/";

///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Radiation ::read(const Io& io) {
    cout << "Reading radiation..." << endl;

    frequencies.read(io);

    if (parameters->use_scattering()) {
        cout << "Using scattering, make sure you have enough memory!" << endl;
    } else {
        cout << "Not using scattering!" << endl;
    }

    // Size and initialize I, u, v
    if (parameters->store_intensities) {
        // Define the local (i.e. in this process) lengths of the distributed arrays
        const Size nrays_local  = pc::message_passing::length(parameters->nrays());
        const Size hnrays_local = pc::message_passing::length(parameters->hnrays());

        I.resize(nrays_local, parameters->npoints(), parameters->nfreqs());
        u.resize(hnrays_local, parameters->npoints(), parameters->nfreqs());
        v.resize(hnrays_local, parameters->npoints(), parameters->nfreqs());
        J.resize(parameters->npoints(), parameters->nfreqs());
    }

    // if (parameters->use_scattering())
    // {
    //     U.resize (parameters->nrays_red());
    //     V.resize (parameters->nrays_red());
    //
    //     for (Size r = 0; r < parameters->nrays_red(); r++)
    //     {
    //         U[r].resize (npoints, nfreqs_red);
    //         V[r].resize (npoints, nfreqs_red);
    //     }
    // }
}

///  write: write out data structure
///    @param[in] io: io object
/////////////////////////////////
void Radiation ::write(const Io& io) const {
    cout << "Writing radiation..." << endl;

    frequencies.write(io);
}

///  initialize: initialize vector with zero's
//////////////////////////////////////////////
void initialize(Real1& vec) {
    threaded_for(i, vec.size(), { vec[i] = 0.0; })
}

///  Initialize J (angular ean intensity)
/////////////////////////////////////////
void Radiation ::initialize_J() {
    for (Size p = 0; p < parameters->npoints(); p++) {
        for (Size f = 0; f < parameters->nfreqs(); f++) {
            J(p, f) = 0.0;
        }
    }
}

/// calc_J: integrate mean intensity and "flux over all directions
//////////////////////////////////////////////////////////////////
void Radiation ::MPI_reduce_J()
#if (MPI_PARALLEL)
{
    int ierr = MPI_Allreduce(MPI_IN_PLACE, // pointer to data to be reduced -> here in place
        J.dat,                             // pointer to data to be received
        J.size(),                          // size of data to be received
        MPI_TYPE_REAL,                     // type of reduced data
        MPI_SUM,                           // reduction operation
        MPI_COMM_WORLD);
    assert(ierr == 0);
}
#else
{
    return;
}
#endif

/// calc_U_and_V: integrate scattering quantities over all directions
/////////////////////////////////////////////////////////////////////
void Radiation ::calc_U_and_V()
#if (MPI_PARALLEL)
    {
        //    Real1 U_local (parameters->npoints()*parameters->nfreqs());
        //    Real1 V_local (parameters->npoints()*parameters->nfreqs());
        //
        //    for (Size w = 0; w < pc::message_passing::comm_size(); w++)
        //    {
        //        const Size start = ( w   *(parameters->hnrays())) /
        //        pc::message_passing::comm_size(); const Size stop  =
        //        ((w+1)*(parameters->hnrays())) / pc::message_passing::comm_size();
        //
        //        for (Size r1 = start; r1 < stop; r1++)
        //        {
        //            const Size R1 = r1 - start;
        //
        //            initialize (U_local);
        //            initialize (V_local);
        //
        //            distributed_for (r2, R2, parameters->hnrays(),
        //            {
        //                threaded_for (p, parameters->npoints(),
        //                {
        //                    for (Size f = 0; f < parameters->nfreqs(); f++)
        //              	    {
        //                        U_local[index(p,f)] += u[R2][index(p,f)]; // *
        //                        scattering.phase[r1][r2][f]; V_local[index(p,f)] +=
        //                        v[R2][index(p,f)]; // * scattering.phase[r1][r2][f];
        //                    }
        //                })
        //            }) // end of r2 loop over raypairs2
        //
        //            int ierr_u = MPI_Reduce (
        //                           U_local.data(),     // pointer to the data to be
        //                           reduced U[R1].data(),       // pointer to the
        //                           data to be received U_local.size(),     // size
        //                           of the data to be received MPI_DOUBLE,         //
        //                           type of the reduced data MPI_SUM,            //
        //                           reduction operation w,                  // rank
        //                           of root to which we reduce MPI_COMM_WORLD);
        //            assert (ierr_u == 0);
        //
        //            int ierr_v = MPI_Reduce (
        //                           V_local.data(),     // pointer to the data to be
        //                           reduced V[R1].data(),       // pointer to the
        //                           data to be received V_local.size(),     // size
        //                           of the data to be received MPI_DOUBLE,         //
        //                           type of the reduced data MPI_SUM,            //
        //                           reduction operation w,                  // rank
        //                           of root to which we reduce MPI_COMM_WORLD);
        //            assert (ierr_v == 0);
        //        }
        //    }
    }
#else
{
    //    Real1 U_local (parameters->npoints()*parameters->nfreqs());
    //    Real1 V_local (parameters->npoints()*parameters->nfreqs());

    //    for (Size r1 = 0; r1 < parameters->hnrays(); r1++)
    //    {
    //        initialize (U_local);
    //        initialize (V_local);

    //        for (Size r2 = 0; r2 < parameters->hnrays(); r2++)
    //        {
    //            threaded_for (p, parameters->npoints(),
    //            {
    //                for (Size f = 0; f < parameters->nfreqs(); f++)
    //                {
    //                    U_local[index(p,f)] += u[r2][index(p,f)]; // *
    //                    scattering.phase[r1][r2][f]; V_local[index(p,f)] +=
    //                    v[r2][index(p,f)]; // * scattering.phase[r1][r2][f];
    //                }
    //            })
    //        }

    //        U[r1] = U_local;
    //        V[r1] = V_local;
    //    }
}
#endif
