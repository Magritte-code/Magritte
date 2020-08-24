// Magritte: Multidimensional Accelerated General-purpose Radiative Transfer
//
// Developed by: Frederik De Ceuster - University College London & KU Leuven
// _________________________________________________________________________


#include <iostream>
#include <iomanip>

#include "radiation.hpp"
#include "Tools/constants.hpp"
#include "Tools/Parallel/wrap_omp.hpp"
#include "Tools/Parallel/wrap_mpi.hpp"
#include "Tools/Parallel/wrap_Grid.hpp"
#include "Tools/logger.hpp"

const string Radiation::prefix = "Radiation/";


///  read: read in data structure
///    @param[in] io: io object
///    @param[in] parameters: model parameters object
/////////////////////////////////////////////////////

void Radiation :: read (const Io &io, Parameters &parameters)
{
    cout << "Reading radiation..." << endl;

    frequencies.read (io, parameters);

    ncells     = parameters.ncells     ();
    nrays      = parameters.nrays      ();
    nfreqs     = parameters.nfreqs     ();
    nfreqs_red = parameters.nfreqs_red ();
    nboundary  = parameters.nboundary  ();


    use_scattering = parameters.use_scattering ();


    if (use_scattering) {cout << "Using scattering, make sure you have enough memory!" << endl;}
    else                {cout << "Not using scattering!"                               << endl;}

    nrays_red = MPI_length (nrays/2);

    parameters.set_nrays_red (nrays_red);


    J.resize (ncells*nfreqs_red);


    // Size and initialize I_bdy, u, v, U and V

    I_bdy.resize (nrays_red);

    for (size_t r = 0; r < nrays_red; r++)
    {
        I_bdy[r].resize (nboundary);

        for (size_t p = 0; p < nboundary; p++)
        {
            I_bdy[r][p].resize (nfreqs_red);
        }
    }


  if (use_scattering)
  {
    u.resize (nrays_red);
//    v.resize (nrays_red);

//    U.resize (nrays_red);
//    V.resize (nrays_red);

    for (size_t r = 0; r < nrays_red; r++)
    {
      u[r].resize (ncells*nfreqs_red);
//      v[r].resize (ncells*nfreqs_red);

//      U[r].resize (ncells*nfreqs_red);
//      V[r].resize (ncells*nfreqs_red);
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

    cout << "Writing u's..." << endl;

    io.write_array (prefix+"u", (Double2) u);

    cout << "u's written..." << endl;


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

int initialize (vReal1 &vec)
{
    OMP_PARALLEL_FOR (i, vec.size())
    {
        vec[i] = 0.0;
    }

    return (0);
}


int Radiation :: initialize_J ()
{
    initialize (J);

    return (0);
}


/// mpi_vector_sum: custom reduction operation using vReals for MPI_Reduce
//////////////////////////////////////////////////////////////////////////

#if (MPI_PARALLEL)

void mpi_vector_sum (vReal *in, vReal *inout, int *len, MPI_Datatype *datatype)
{
    for (int i = 0; i < *len; i++)
    {
        inout[i] = in[i] + inout[i];
    }
}

#endif




/// calc_J: integrate mean intensity and "flux over all directions
//////////////////////////////////////////////////////////////////

int Radiation :: MPI_reduce_J ()

#if (MPI_PARALLEL)

{
    MPI_Datatype MPI_VREAL;
    MPI_Type_contiguous (n_simd_lanes, MPI_DOUBLE, &MPI_VREAL);
    MPI_Type_commit (&MPI_VREAL);

    MPI_Op MPI_VSUM;
    MPI_Op_create ((MPI_User_function*) mpi_vector_sum, true, &MPI_VSUM);

    int ierr = MPI_Allreduce (
                  MPI_IN_PLACE,      // pointer to data to be reduced -> here in place
                  J.data(),          // pointer to data to be received
                  J.size(),          // size of data to be received
                  MPI_VREAL,         // type of reduced data
                  MPI_VSUM,          // reduction operation
                  MPI_COMM_WORLD);
    assert (ierr == 0);

    MPI_Type_free (&MPI_VREAL);
    MPI_Op_free   (&MPI_VSUM);

    return (0);
}

#else

{
    return (0);
}

#endif




/// calc_U_and_V: integrate scattering quantities over all directions
/////////////////////////////////////////////////////////////////////

int Radiation :: calc_U_and_V ()

#if (MPI_PARALLEL)

{
    vReal1 U_local (ncells*nfreqs_red);
    vReal1 V_local (ncells*nfreqs_red);


    MPI_Datatype MPI_VREAL;
    MPI_Type_contiguous (n_simd_lanes, MPI_DOUBLE, &MPI_VREAL);
    MPI_Type_commit (&MPI_VREAL);

    MPI_Op MPI_VSUM;
    MPI_Op_create ( (MPI_User_function*) mpi_vector_sum, true, &MPI_VSUM);


    for (int w = 0; w < MPI_comm_size(); w++)
    {
        const size_t start = ( w   *(nrays/2)) / MPI_comm_size();
        const size_t stop  = ((w+1)*(nrays/2)) / MPI_comm_size();

        for (size_t r1 = start; r1 < stop; r1++)
        {
            const size_t R1 = r1 - start;

            initialize (U_local);
            initialize (V_local);

            MPI_PARALLEL_FOR (r2, nrays/2)
            {
                const size_t R2 = r2 - MPI_start (nrays/2);

                OMP_PARALLEL_FOR (p, ncells)
                {
                    for (size_t f = 0; f < nfreqs_red; f++)
              	    {
                        U_local[index(p,f)] += u[R2][index(p,f)] ;//* scattering.phase[r1][r2][f];
                        V_local[index(p,f)] += v[R2][index(p,f)] ;//* scattering.phase[r1][r2][f];
                    }
                }
            } // end of r2 loop over raypairs2

            int ierr_u = MPI_Reduce (
                           U_local.data(),     // pointer to the data to be reduced
                           U[R1].data(),       // pointer to the data to be received
                           ncells*nfreqs_red,  // size of the data to be received
                           MPI_VREAL,          // type of the reduced data
                           MPI_VSUM,           // reduction operation
                           w,                  // rank of root to which we reduce
                           MPI_COMM_WORLD);
            assert (ierr_u == 0);

            int ierr_v = MPI_Reduce (
                           V_local.data(),     // pointer to the data to be reduced
                           V[R1].data(),       // pointer to the data to be received
                           ncells*nfreqs_red,  // size of the data to be received
                           MPI_VREAL,          // type of the reduced data
                           MPI_VSUM,           // reduction operation
                           w,                  // rank of root to which we reduce
                           MPI_COMM_WORLD);
            assert (ierr_v == 0);
        }
    }

    MPI_Type_free (&MPI_VREAL);
    MPI_Op_free   (&MPI_VSUM);

    return (0);
}

#else

{
    vReal1 U_local (ncells*nfreqs_red);
    vReal1 V_local (ncells*nfreqs_red);

    for (size_t r1 = 0; r1 < nrays/2; r1++)
    {
        initialize (U_local);
        initialize (V_local);

        for (size_t r2 = 0; r2 < nrays/2; r2++)
        {
            OMP_PARALLEL_FOR (p, ncells)
            {
                for (size_t f = 0; f < nfreqs_red; f++)
                {
                    U_local[index(p,f)] += u[r2][index(p,f)] ;//* scattering.phase[r1][r2][f];
                    V_local[index(p,f)] += v[r2][index(p,f)] ;//* scattering.phase[r1][r2][f];
                }
            }
        }

        U[r1] = U_local;
        V[r1] = V_local;
    }

    return (0);
}

#endif
