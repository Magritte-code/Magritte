inline void Lambda :: initialize (const Size nrad_new)
{
    nrad = nrad_new;

    Lss.reserve (parameters.npoints() * nrad);
    nrs.reserve (parameters.npoints() * nrad);

    size.resize (parameters.npoints() * nrad);

    Ls.resize (parameters.npoints());
    nr.resize (parameters.npoints());

    for (Size p = 0; p < parameters.npoints(); p++)
    {
        Ls[p].resize (nrad);
        nr[p].resize (nrad);

        for (Size k = 0; k < nrad; k++)
        {
            Ls[p][k].reserve (5);
            nr[p][k].reserve (5);

            size [index_first(p,k)] = 0;
        }
    }
}



inline void Lambda :: clear ()
{
    threaded_for (p, parameters.npoints(),
    {
        for (Size k = 0; k < nrad; k++)
        {
            Ls[p][k].clear();
            nr[p][k].clear();
        }
    })
}


/// Index of the first element belonging to p and k
///    @param[in] p : index of the receiving cell
///    @param[in] k : index of the line transition
///////////////////////////////////////////////////
inline Size Lambda :: index_first (const Size p, const Size k) const
{
    return k + nrad*p;
}


/// Index of the last element belonging to p and k
///    @param[in] p : index of the receiving cell
///    @param[in] k : index of the line transition
//////////////////////////////////////////////////
inline Size Lambda :: index_last (const Size p, const Size k) const
{
    return index_first(p,k) + get_size(p,k) - 1;
}


///  Getter for ALO element
///    @param[in] p      : index of the receiving cell
///    @param[in] k      : index of the line transition
///    @param[in] index  : index of the ALO element
///////////////////////////////////////////////////////
inline Real Lambda :: get_Ls (const Size p, const Size k, const Size index) const
{
    return Ls[p][k][index];
}


///  Getter for cell index corresponding of ALO element
///    @param[in] p      : index of the receiving cell
///    @param[in] k      : index of the line transition
///    @param[in] index  : index of the ALO element
///////////////////////////////////////////////////////
inline Size Lambda :: get_nr (const Size p, const Size k, const Size index) const
{
    return nr[p][k][index];
}


///  Getter for cell index corresponding of ALO element
///    @param[in] p      : index of the receiving cell
///    @param[in] k      : index of the line transition
///    @param[in] index  : index of the ALO element
///////////////////////////////////////////////////////
inline Size Lambda :: get_size (const Size p, const Size k) const
{
    return nr[p][k].size();
}


///  Setter for an ALO element
///    @param[in] p      : index of the receiving cell
///    @param[in] k      : index of the line transition
///    @param[in] nr_new : index of the emitting cell
///    @param[in] Ls_new : new element of the ALO
///////////////////////////////////////////////////////
inline void Lambda :: add_element (const Size p, const Size k, const Size nr_new, const Real Ls_new)
{
    for (Size index = 0; index < nr[p][k].size(); index++)
    {
        if (nr[p][k][index] == nr_new)
        {
            Ls[p][k][index] += Ls_new;
            return;
        }
    }

    Ls[p][k].push_back (Ls_new);
    nr[p][k].push_back (nr_new);
}




inline void Lambda :: linearize_data ()
{
    Size size_total = 0;

#   pragma omp parallel for reduction (+: size_total)
    for (Size p = 0; p < parameters.npoints(); p++)
    {
        for (Size k = 0; k < nrad; k++)
        {
            const Size index = index_first (p,k);

            size[index] = nr[p][k].size();
            size_total += size[index];
        }
    }

    Lss.resize (size_total);
    nrs.resize (size_total);

    Size index = 0;

    for (Size p = 0; p < parameters.npoints(); p++)
    {
        for (Size k = 0; k < nrad; k++)
        {
            for (Size m = 0; m < size[index_first(p,k)]; m++)
            {
                Lss[index] = Ls[p][k][m];
                nrs[index] = nr[p][k][m];

                index++;
            }
        }
    }
}




inline void Lambda :: MPI_gather ()

#if (MPI_PARALLEL)

{
    // Linearize the Lambda operator data
    linearize_data ();

    int size_total = Lss.size();


    // Gather the lengths of the linearized vectors of each process
    Int1 buffer_lengths (pc::message_passing::comm_size(), 0);
    Int1 displacements  (pc::message_passing::comm_size(), 0);


    int ierr_l = MPI_Allgather (
                     &size_total,             // pointer to data to be send
                     1,                       // number of elements in the send buffer
                     MPI_INT,                 // type of the send data
                     buffer_lengths.data(),   // pointer to the data to be received
                     1,                       // number of elements in receive buffer
                     MPI_INT,                 // type of the received data
                     MPI_COMM_WORLD        );
    assert (ierr_l == 0);


    for (int w = 1; w < pc::message_passing::comm_size(); w++)
    {
        displacements[w] = buffer_lengths[w-1];
    }

    Real1 Lss_total;
    Size1 nrs_total;
    Size1 szs_total;

    Size total_buffer_length = 0;

    for (Size length : buffer_lengths) {total_buffer_length += length;}

    Lss_total.resize (total_buffer_length);
    nrs_total.resize (total_buffer_length);
    szs_total.resize (pc::message_passing::comm_size()*parameters.npoints()*nrad);


    int ierr_ls = MPI_Allgatherv (
                      Lss.data(),              // pointer to data to be send
                      size_total,              // number of elements in the send buffer
                      MPI_DOUBLE,              // type of the send data
                      Lss_total.data(),        // pointer to the data to be received
                      buffer_lengths.data(),   // list of numbers of elements in receive buffer
                      displacements.data(),    // displacements between data blocks
  	                  MPI_DOUBLE,              // type of the received data
                      MPI_COMM_WORLD);
    assert (ierr_ls == 0);


    int ierr_nr = MPI_Allgatherv (
                      nrs.data(),              // pointer to data to be send
                      size_total,              // number of elements in the send buffer
                      MPI_LONG,                // type of the send data
                      nrs_total.data(),        // pointer to the data to be received
                      buffer_lengths.data(),   // list of numbers of elements in receive buffer
                      displacements.data(),    // displacements between data blocks
  	                  MPI_LONG,                // type of the received data
                      MPI_COMM_WORLD);
    assert (ierr_nr == 0);


    int ierr_sz = MPI_Allgather (
                      size.data(),               // pointer to data to be send
                      parameters.npoints()*nrad, // number of elements in the send buffer
                      MPI_LONG,                  // type of the send data
                      szs_total.data(),          // pointer to the data to be received
                      parameters.npoints()*nrad, // number of elements in receive buffer
                      MPI_LONG,                  // type of the received data
                      MPI_COMM_WORLD);
    assert (ierr_sz == 0);


    clear();


    Size index_1 = 0;
    Size index_2 = 0;

    for (Size w = 0; w < pc::message_passing::comm_size(); w++)
    {
        for (Size p = 0; p < parameters.npoints(); p++)
        {
            for (Size k = 0; k < nrad; k++)
            {
                for (Size m = 0; m < szs_total[index_1]; m++)
                {
                    add_element (p, k, nrs_total[index_2], Lss_total[index_2]);
                    index_2++;
                }
                index_1++;
            }
        }
    }
}

#else

{
    return;
}

#endif
