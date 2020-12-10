#include <iostream>
#include <iomanip>

#include "radiation.hpp"
#include "tools/constants.hpp"

const string prefix = "radiation/";


///  read: read in data structure
///    @param[in] io: io object
/////////////////////////////////
void Radiation :: read (const Io& io)
{
    cout << "Reading radiation..." << endl;

    frequencies.read (io);

    if (parameters.use_scattering())
    {
        cout << "Using scattering, make sure you have enough memory!" << endl;
    }
    else
    {
        cout << "Not using scattering!" << endl;
    }

    parameters.set_nrays_red (paracabs::message_passing::length (parameters.hnrays()));




    // Size and initialize I_bdy, u, v, U and V

    I_bdy.resize (parameters.nrays_red());

    for (Size r = 0; r < parameters.nrays_red(); r++)
    {
        I_bdy[r].resize (parameters.nboundary());

        for (Size p = 0; p < parameters.nboundary(); p++)
        {
            I_bdy[r][p].resize (parameters.nfreqs());
        }
    }


    I.resize (parameters.nrays(),  parameters.npoints(), parameters.nfreqs());
    u.resize (parameters.hnrays(), parameters.npoints(), parameters.nfreqs());
    v.resize (parameters.hnrays(), parameters.npoints(), parameters.nfreqs());
    J.resize (                     parameters.npoints(), parameters.nfreqs());

    // for (Size r = 0; r < parameters.nrays(); r++)
    // {
        // I[r].resize (parameters.npoints(), parameters.nfreqs());
    // }

    if (parameters.use_scattering())
    {

        // u.resize (parameters.nrays_red());
        // v.resize (parameters.nrays_red());

        // U.resize (parameters.nrays_red());
        // V.resize (parameters.nrays_red());

        for (Size r = 0; r < parameters.nrays_red(); r++)
        {
            // u[r].resize (parameters.npoints(), parameters.nfreqs());
//            v[r].resize (npoints, nfreqs_red);

//            U[r].resize (npoints, nfreqs_red);
//            V[r].resize (npoints, nfreqs_red);
        }
    }
}


///  write: write out data structure
///    @param[in] io: io object
/////////////////////////////////
void Radiation :: write (const Io &io) const
{
    cout << "Writing radiation..." << endl;

    frequencies.write (io);


//    Double2 JJ (ncells, Double1 (nfreqs));
//
//    OMP_PARALLEL_FOR (p, ncells)
//    {
//        for (size_t f = 0; f < nfreqs; f++)
//        {
//            JJ[p][f] = get_J (p,f);
//        }
//    }
//
//    io.write_array (prefix+"J", JJ);

//    io.write_list (prefix+"J", J);


//    Double2 uu (nrays_red, Double1 (ncells*nfreqs));
//
//    for (size_t R = 0; R < nrays_red; R++)
//    {
//        OMP_PARALLEL_FOR (p, ncells)
//        {
//            for (size_t f = 0; f < nfreqs; f++)
//            {
//                uu[R][index(p,f)] = get_u (R,p,f);
//            }
//        }
//    }
//
//    io.write_array (prefix+"u", uu);


    // io.write_array (prefix+"u", u);


//    Double2 I_bdy_buffer (I_bdy.size(), Double1 (I_bdy[0].size()*I_bdy[0][0].size()));
//
//    for (long r = 0; r < nrays_red; r++)
//    {
//        OMP_PARALLEL_FOR (p, ncells)
//        {
//            for (size_t f = 0; f < nfreqs; f++)
//            {
//                I_bdy_buffer[r][p + f] = get_J (p,f);
//            }
//        }
//    }


  // Print all frequencies (nu)
//# if (GRID_SIMD)
//
//    Double3 u_expanded (ncells, Double1 (ncells*nfreqs));
//
//
//    OMP_PARALLEL_FOR (p, ncells)
//    {
//      long index = 0;
//
//      for (long f = 0; f < nfreqs_red; f++)
//      {
//        for (int lane = 0; lane < n_simd_lanes; lane++)
//        {
//          nu_expanded[p][index] = nu[p][f].getlane (lane);
//          index++;
//        }
//      }
//    }
//
//    io.write_array (prefix+"nu", nu_expanded);
//
//# else
//
//    io.write_array (prefix+"nu", nu);
//
//# endif
//  if (MPI_comm_rank () == 0)
//  {
//    const string file_name_J = output_folder + "J" + tag + ".txt";
//    const string file_name_G = output_folder + "G" + tag + ".txt";
//
//    ofstream outputFile_J (file_name_J);
//    ofstream outputFile_G (file_name_G);
//
//    outputFile_J << scientific << setprecision(16);
//    outputFile_G << scientific << setprecision(16);
//
//
//    for (long p = 0; p < ncells; p++)
//    {
//      for (int f = 0; f < nfreq_red; f++)
//      {
//#       if (GRID_SIMD)
//          for (int lane = 0; lane < n_simd_lanes; lane++)
//          {
//            outputFile_J << J[index(p,f)].getlane(lane) << "\t";
//            outputFile_G << G[index(p,f)].getlane(lane) << "\t";
//          }
//#       else
//          outputFile_J << J[index(p,f)] << "\t";
//          outputFile_G << G[index(p,f)] << "\t";
//#       endif
//      }
//
//      outputFile_J << endl;
//      outputFile_G << endl;
//    }
//
//    outputFile_J.close ();
//    outputFile_G.close ();
//
//
//    const string file_name_bc = output_folder + "bc" + tag + ".txt";
//
//    ofstream outputFile_bc (file_name_bc);
//
//    //for (long r = 0; r < nrays_red; r++)
//    //{
//      long r = 0;
//      for (long b = 0; b < nboundary; b++)
//      {
//        for (long f = 0; f < nfreq_red; f++)
//        {
//#         if (GRID_SIMD)
//            for (int lane = 0; lane < n_simd_lanes; lane++)
//            {
//              outputFile_bc << boundary_intensity[r][b][f].getlane(lane) << "\t";
//            }
//#         else
//            outputFile_bc << boundary_intensity[r][b][f] << "\t";
//#         endif
//        }
//          outputFile_bc << endl;
//      }
//    //}
//
//    outputFile_bc.close ();
//
//  }
//
//
}




///  initialize: initialize vector with zero's
//////////////////////////////////////////////

void initialize (Real1 &vec)
{
    threaded_for (i, vec.size(),
    {
        vec[i] = 0.0;
    })
}


void Radiation :: initialize_J ()
{
    for (Size p = 0; p < parameters.npoints(); p++)
    {
        for (Size f = 0; f < parameters.nfreqs(); f++)
        {
            J(p,f) = 0.0;
        }
    }
}


/// calc_J: integrate mean intensity and "flux over all directions
//////////////////////////////////////////////////////////////////
void Radiation :: MPI_reduce_J ()

#if (MPI_PARALLEL)

{
    int ierr = MPI_Allreduce (
                  MPI_IN_PLACE,      // pointer to data to be reduced -> here in place
                  J.data(),          // pointer to data to be received
                  J.size(),          // size of data to be received
                  MPI_DOUBLE,        // type of reduced data
                  MPI_SUM,           // reduction operation
                  MPI_COMM_WORLD);
    assert (ierr == 0);
}

#else

{
    return;
}

#endif




/// calc_U_and_V: integrate scattering quantities over all directions
/////////////////////////////////////////////////////////////////////
void Radiation :: calc_U_and_V ()

#if (MPI_PARALLEL)

{
//    Real1 U_local (parameters.npoints()*parameters.nfreqs());
//    Real1 V_local (parameters.npoints()*parameters.nfreqs());
//
//    for (Size w = 0; w < pc::message_passing::comm_size(); w++)
//    {
//        const Size start = ( w   *(parameters.hnrays())) / pc::message_passing::comm_size();
//        const Size stop  = ((w+1)*(parameters.hnrays())) / pc::message_passing::comm_size();
//
//        for (Size r1 = start; r1 < stop; r1++)
//        {
//            const Size R1 = r1 - start;
//
//            initialize (U_local);
//            initialize (V_local);
//
//            distributed_for (r2, R2, parameters.hnrays(),
//            {
//                threaded_for (p, parameters.npoints(),
//                {
//                    for (Size f = 0; f < parameters.nfreqs(); f++)
//              	    {
//                        U_local[index(p,f)] += u[R2][index(p,f)]; // * scattering.phase[r1][r2][f];
//                        V_local[index(p,f)] += v[R2][index(p,f)]; // * scattering.phase[r1][r2][f];
//                    }
//                })
//            }) // end of r2 loop over raypairs2
//
//            int ierr_u = MPI_Reduce (
//                           U_local.data(),     // pointer to the data to be reduced
//                           U[R1].data(),       // pointer to the data to be received
//                           U_local.size(),     // size of the data to be received
//                           MPI_DOUBLE,         // type of the reduced data
//                           MPI_SUM,            // reduction operation
//                           w,                  // rank of root to which we reduce
//                           MPI_COMM_WORLD);
//            assert (ierr_u == 0);
//
//            int ierr_v = MPI_Reduce (
//                           V_local.data(),     // pointer to the data to be reduced
//                           V[R1].data(),       // pointer to the data to be received
//                           V_local.size(),     // size of the data to be received
//                           MPI_DOUBLE,         // type of the reduced data
//                           MPI_SUM,            // reduction operation
//                           w,                  // rank of root to which we reduce
//                           MPI_COMM_WORLD);
//            assert (ierr_v == 0);
//        }
//    }
}

#else

{
//    Real1 U_local (parameters.npoints()*parameters.nfreqs());
//    Real1 V_local (parameters.npoints()*parameters.nfreqs());

//    for (Size r1 = 0; r1 < parameters.hnrays(); r1++)
//    {
//        initialize (U_local);
//        initialize (V_local);

//        for (Size r2 = 0; r2 < parameters.hnrays(); r2++)
//        {
//            threaded_for (p, parameters.npoints(),
//            {
//                for (Size f = 0; f < parameters.nfreqs(); f++)
//                {
//                    U_local[index(p,f)] += u[r2][index(p,f)]; // * scattering.phase[r1][r2][f];
//                    V_local[index(p,f)] += v[r2][index(p,f)]; // * scattering.phase[r1][r2][f];
//                }
//            })
//        }

//        U[r1] = U_local;
//        V[r1] = V_local;
//    }
}

#endif
