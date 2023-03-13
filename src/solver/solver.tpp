#include <set>

template <Frame frame>
inline void Solver :: setup (Model& model)
{
    const Size length = 2 * get_ray_lengths_max <frame> (model) + 1;
    const Size  width = model.parameters->nfreqs();
    const Size  n_o_d = model.parameters->n_off_diag;

    model.set_dshift_max();//err, probably belongs somewhere else, but we need to compute the max shift for each point

    setup (length, width, n_o_d);
}


inline void Solver :: setup (const Size l, const Size w, const Size n_o_d)
{
    length     = l;
    centre     = l/2;
    width      = w;
    n_off_diag = n_o_d;

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        dZ_          (i).resize (length);
        nr_          (i).resize (length);
        shift_       (i).resize (length);

        eta_c_       (i).resize (width);
        eta_n_       (i).resize (width);

        chi_c_       (i).resize (width);
        chi_n_       (i).resize (width);

        source_c_    (i).resize (width);
        source_n_    (i).resize (width);

        inverse_chi_ (i).resize (length);

        tau_         (i).resize (width);

        Su_          (i).resize (length);
        Sv_          (i).resize (length);

        A_           (i).resize (length);
        C_           (i).resize (length);
        inverse_A_   (i).resize (length);
        inverse_C_   (i).resize (length);

        FF_          (i).resize (length);
        FI_          (i).resize (length);
        GG_          (i).resize (length);
        GI_          (i).resize (length);
        GP_          (i).resize (length);

        L_diag_      (i).resize (length);

        L_upper_     (i).resize (n_off_diag, length);
        L_lower_     (i).resize (n_off_diag, length);
    }
}


//l=max length of ray, w=nfreqs, n_o_d is number of off-diagonals; will be ignored
// inline void Solver :: setup_comoving (Model& model, const Size l, const Size w)
inline void Solver :: setup_comoving (Model& model)
{
    length     = 2 * get_ray_lengths_max <Rest> (model) + 1;
    centre     = length/2;
    width      = model.parameters->nfreqs();

    model.set_dshift_max();//err, probably belongs somewhere else, but we need to compute the max shift for each point

    std::cout<<"setup length: "<<length<<std::endl;
    std::cout<<"setup width: "<<width<<std::endl;


    points_to_trace_ray_through.resize(model.parameters->hnrays());
    for (Size i = 0; i < model.parameters->hnrays(); i ++)
    {
        points_to_trace_ray_through[i].resize(model.parameters->npoints());
    }

    //For determining which ray lies closest to each point
    n_rays_through_point.resize(model.parameters->hnrays(), model.parameters->npoints());
    min_ray_distsqr     .resize(model.parameters->hnrays(), model.parameters->npoints());
    closest_ray         .resize(model.parameters->hnrays(), model.parameters->npoints());

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        //general ray tracing variables
        dZ_          (i).resize (length);
        nr_          (i).resize (length);
        shift_       (i).resize (length);

        start_indices_(i).resize(length, width);
        for (Size j=0; j<length; j++)
        {
            for (Size k=0; k<width; k++)
            {
                start_indices_(i)(j,k).resize(2);
            }
        }

        line_quad_discdir_overlap_(i).resize(model.parameters->nlines());
        left_bound_(i).resize(model.parameters->nlines());
        right_bound_(i).resize(model.parameters->nlines());
        left_bound_index_(i).resize(model.parameters->nlines());
        right_bound_index_(i).resize(model.parameters->nlines());

        line_count_(i).resize(model.parameters->nlines());
        quad_range_weight_(i).resize(model.parameters->nquads());

        intensities_(i).resize(length, width);
        delta_tau_(i).resize(length, width);
        S_curr_(i).resize(length, width);
        S_next_(i).resize(length, width);

        dIdnu_coef1_curr_(i).resize(length, width);
        dIdnu_coef2_curr_(i).resize(length, width);
        dIdnu_coef3_curr_(i).resize(length, width);

        dIdnu_index1_curr_(i).resize(length, width);
        dIdnu_index2_curr_(i).resize(length, width);
        dIdnu_index3_curr_(i).resize(length, width);

        dIdnu_coef1_next_(i).resize(length, width);
        dIdnu_coef2_next_(i).resize(length, width);
        dIdnu_coef3_next_(i).resize(length, width);

        dIdnu_index1_next_(i).resize(length, width);
        dIdnu_index2_next_(i).resize(length, width);
        dIdnu_index3_next_(i).resize(length, width);

    }

    //Finally also trace the rays in advance to prune the unnecessary ones.
    get_static_rays_to_trace(model);

}


//l=max length of ray, w=nfreqs, n_o_d is number of off-diagonals; will be ignored
// Approximated version of comoving solver requires less setup stuff
inline void Solver :: setup_comoving_local_approx (Model& model)
{
    length     = 2 * get_ray_lengths_max <Rest> (model) + 1;
    centre     = length/2;
    width      = model.parameters->nfreqs();

    model.set_dshift_max();//err, probably belongs somewhere else, but we need to compute the max shift for each point

    std::cout<<"setup length: "<<length<<std::endl;
    std::cout<<"setup width: "<<width<<std::endl;

    //Setup for allowing a ray to compute multiple values
    points_to_trace_ray_through.resize(model.parameters->hnrays());
    for (Size i = 0; i < model.parameters->hnrays(); i ++)
    {
        points_to_trace_ray_through[i].resize(model.parameters->npoints());
    }

    //for determining which ray lies closest to each point
    n_rays_through_point.resize(model.parameters->hnrays(), model.parameters->npoints());
    min_ray_distsqr     .resize(model.parameters->hnrays(), model.parameters->npoints());
    closest_ray         .resize(model.parameters->hnrays(), model.parameters->npoints());

    for (Size i = 0; i < pc::multi_threading::n_threads_avail(); i++)
    {
        //general ray tracing variables
        dZ_          (i).resize (length);
        nr_          (i).resize (length);
        shift_       (i).resize (length);

        //For simplicity, we define other variables for the comoving approx solver, as we store less (in full solver, we use matrices, here vectors suffice)
        cma_computed_intensities_(i).resize(width);
        cma_start_intensities_(i).resize(width);
        cma_start_frequencies_(i).resize(width);
        cma_start_frequency_index_(i).resize(width);
        cma_end_frequencies_(i).resize(width);
        // cma_dIdnu_curr_(i).resize(width); explicit, computed intensity derivative with respect to frequency; not used
        cma_chi_curr_(i).resize(width);
        cma_chi_next_(i).resize(width);
        cma_compute_curr_opacity_(i).resize(width);
        cma_compute_next_opacity_(i).resize(width);
        cma_S_curr_(i).resize(width);
        cma_S_next_(i).resize(width);

        cma_dIdnu_coef1_next_(i).resize(width);
        cma_dIdnu_coef2_next_(i).resize(width);
        cma_dIdnu_coef3_next_(i).resize(width);
        cma_dIdnu_expl_(i).resize(width);
        cma_delta_tau_(i).resize(width);

        // The implicit frequency derivative will be computed on the fly
        // cma_dIdnu_index1_next_(i).resize(width);
        // cma_dIdnu_index2_next_(i).resize(width);
        // cma_dIdnu_index3_next_(i).resize(width);
    }

    //Finally also trace the rays in advance to prune the unnecessary ones.
    get_static_rays_to_trace(model);

}



///  Getter for the maximum allowed shift value determined by the smallest line
///    @param[in] o : number of point under consideration
///    @retrun maximum allowed shift value determined by the smallest line
///////////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: get_dshift_max (
    const Model& model,
    const Size   o     )
{
    Real dshift_max = std::numeric_limits<Real>::max();

    for (const LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        const Real inverse_mass   = lspec.linedata.inverse_mass;
        const Real new_dshift_max = model.parameters->max_width_fraction
                                    * model.thermodynamics.profile_width (inverse_mass, o);

        if (dshift_max > new_dshift_max)
        {
            dshift_max = new_dshift_max;
        }
    }

    return dshift_max;
}


template <Frame frame>
inline void Solver :: get_ray_lengths (Model& model)
{
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        accelerated_for (o, model.parameters->npoints(),
        {
            const Real dshift_max = get_dshift_max (model, o);

            model.geometry.lengths(rr,o) =
                model.geometry.get_ray_length <frame> (o, rr, dshift_max)
              + model.geometry.get_ray_length <frame> (o, ar, dshift_max);
        })

        pc::accelerator::synchronize();
    }

    model.geometry.lengths.copy_ptr_to_vec();
}


template <Frame frame>
inline Size Solver :: get_ray_lengths_max (Model& model)
{
    get_ray_lengths <frame> (model);

    Geometry& geo = model.geometry;

    geo.lengths_max = *std::max_element(geo.lengths.vec.begin(),
                                        geo.lengths.vec.end()   );

    return geo.lengths_max;
}


///  Traces all points on a given ray and checks whether the ray lies closest (among currently traced rays) to the given ray
///  Should be called in both directions, but as tracing antipodal rays should be symmetric, we only need to keep track of one direction
accel inline void Solver :: trace_ray_points (
    const Geometry& geometry,
    const Size      o,//origin point of ray
    const Size      rdir,//ray direction to trace ∈ [0, nrays-1]
    const Size      rsav,//ray index save direction ∈ [0,hnrays-1]
    const Size      rayidx)//for indexing the rays (counts from 0 to total number of rays traced through the domain
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = geometry.get_next (o, rdir, o, Z, dZ);

    if (geometry.valid_point(nxt))
    {
        // get distance and check if closest ray
        Real dist2 = geometry.get_dist2_ray_point(o, nxt, rdir);
        //If it is the first time we encounter this point, or this is the closest ray: assign this ray to compute the stuff (J, lambda) of the point
        if (n_rays_through_point(rsav,nxt)==0||dist2<min_ray_distsqr(rsav,nxt))
        {
            min_ray_distsqr(rsav,nxt)=dist2;
            closest_ray(rsav,nxt)=rayidx;
        }

        n_rays_through_point(rsav,nxt)++;

        Size crt = o;

        while (geometry.not_on_boundary(nxt))
        {
            crt =       nxt;
            nxt = geometry.get_next(o, rdir, nxt, Z, dZ);

            // get distance and check if closest ray
            Real dist2 = geometry.get_dist2_ray_point(o, nxt, rdir);

            if (n_rays_through_point(rsav,nxt)==0||dist2<min_ray_distsqr(rsav,nxt))
            {
                min_ray_distsqr(rsav,nxt)=dist2;
                closest_ray(rsav,nxt)=rayidx;
            }

            n_rays_through_point(rsav,nxt)++;
        }
    }
}




// For all directions, determines a ray covering of the points
inline void Solver :: get_static_rays_to_trace (Model& model)
{
    accelerated_for (rr, model.parameters->hnrays(),
    {
        Size n_rays_to_trace=0;
        const Size ar=model.geometry.rays.antipod[rr];
        //To make sure that we order the rays starting from the largest elements
        for (Size pointidx=0; pointidx<model.parameters->npoints(); pointidx++)
        {
            const Size o=model.geometry.sorted_position_indices[pointidx];

            //DEBUG: commenting this out results in tracing through all points
            if (n_rays_through_point(rr,o)>0)
            {continue;}

            //seemingly no ray has been traced through this point, so we must trace a ray through it
            points_to_trace_ray_through[rr][n_rays_to_trace]=o;
            //tracing rays is symmetric, so only keep for the first half of the ray directions

            //trace ray through point
            n_rays_through_point(rr,o)++;
            //antipod has exactly same number of rays through the point, so do not save

            // For generality, I assign a ray index to each ray of a given direction rsav
            // Also the ray direction is identified with the forward dir
            closest_ray(rr,o)=n_rays_to_trace;
            min_ray_distsqr(rr,o)=0.0;

            //now trace ray through rest of model
            trace_ray_points(model.geometry, o, rr, rr, n_rays_to_trace);
            trace_ray_points(model.geometry, o, ar, rr, n_rays_to_trace);

            n_rays_to_trace++;
        }
        // n_points_to_trace_ray_through[rr]=n_rays_to_trace;
        points_to_trace_ray_through[rr].resize(n_rays_to_trace);//and now the correct size, instead of parameters->npoints()

    })

    //debug print stuff
    for (Size rr=0; rr<model.parameters->hnrays(); rr++)
    {
        std::cout<<"rr: "<<rr<<" size points_to_trace_ray_through: "<<points_to_trace_ray_through[rr].size()<<std::endl;
        for (Size idx=0; idx<points_to_trace_ray_through[rr].size(); idx++)
        {
            // std::cout<<"point: "<<points_to_trace_ray_through[rr][idx]<<std::endl;
        }
        // std::cout<<"number of rays per point"<<std::endl;
        for (Size p=0; p<model.parameters->npoints(); p++)
        {
            // std::cout<<"closest ray: "<<closest_ray(rr,p)<<std::endl;
            // std::cout<<"point: "<<p<<"#: "<<n_rays_through_point(rr, p)<<std::endl;
        }
    }
}


///  Matches the sorted frequency indices such that the frequency difference is minimal
///  Internal function for the non-approximate comoving solver
// uses the doppler shifts (curr_shift, next_shift)
// ERR, should just iterate over all frequencies, ignoring those different splits. Boundary conditions
// is_upward_disc denotes whether we want an upward discretization
// nextpointonrayindex should denote how far the next point lies on the ray
// currpointonrayindex should do the same, but for the previous point
inline void Solver :: match_frequency_indices(Model& model, const Size nextpoint, const Size currpoint, const Real next_shift, const Real curr_shift, Size nextpointonrayindex, Size currpointonrayindex, bool is_upward_disc)
{
    //if upward discretization, we start from the uppermost part
    if (is_upward_disc)
    {
        //Starting from the highest frequency
        //+1 to all indices due to using unsigned ints when looping down (overflow otherwise)
        Size curr_freq_idx=model.parameters->nfreqs()-1+1;
        Size next_freq_idx=model.parameters->nfreqs()-1+1;
        const Size min_freq_idx=0+1;
        //assumes at least a single frequency will not be a boundary condition (is reasonable if we limit the doppler shift)
        //NOTE: oob freqs will be dealt with later on
        while (curr_freq_idx>=min_freq_idx && next_freq_idx>=min_freq_idx)
        {
            //check if the frequency at the previous point index is higher than this frequency (in static frame)
            if (model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx-1)*curr_shift>
                model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx-1)*next_shift)
            {
                curr_freq_idx--;
            }
            else//freq matching index must be determined if it is just higher than the freq at the previous point
            {
                start_indices_()(nextpointonrayindex,next_freq_idx-1)[0]=currpointonrayindex;//stores point on ray index
                start_indices_()(nextpointonrayindex,next_freq_idx-1)[1]=curr_freq_idx-1;//stores freq
                //for the sake of boundary conditions, these start_indices_ might be overwritten later
                next_freq_idx--;
            }
        }
        //match remaining points as best as possible (with the outermost point)
        while (next_freq_idx>=min_freq_idx)
        {
            start_indices_()(nextpointonrayindex,next_freq_idx-1)[0]=currpointonrayindex;//stores point
            start_indices_()(nextpointonrayindex,next_freq_idx-1)[1]=min_freq_idx-1;//stores freq
            next_freq_idx--;
        }
    }
    else
    {
        //Starting from the lowest frequency
        Size curr_freq_idx = 0;
        Size next_freq_idx = 0;
        const Size max_freq_idx = model.parameters->nfreqs()-1;
        //assumes at least a single frequency will not be a boundary condition (is reasonable if we limit the doppler shift)
        while (curr_freq_idx<=max_freq_idx && next_freq_idx<=max_freq_idx)
        {
            //check if a previous frequency index exists which is higher than this (in static frame)
            if (model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx)*curr_shift<
                model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx)*next_shift)
            {
                curr_freq_idx++;
            }
            else//freq matching index must be determined if it is just higher than the freq at the previous point
            {
                start_indices_()(nextpointonrayindex,next_freq_idx)[0]=currpointonrayindex;//stores point
                start_indices_()(nextpointonrayindex,next_freq_idx)[1]=curr_freq_idx;//stores freq
                //for the sake of boundary conditions, these start_indices_ might be overwritten later
                next_freq_idx++;
            }
        }
        //match remaining points as best as possible (with the outermost point)
        while (next_freq_idx<=max_freq_idx)
        {
            start_indices_()(nextpointonrayindex,next_freq_idx)[0]=currpointonrayindex;//stores point
            start_indices_()(nextpointonrayindex,next_freq_idx)[1]=max_freq_idx;//stores freq
            next_freq_idx++;
        }
    }
}


///  Computes which lines are overlapping at the given point (at the next point)
///  Internal function for (non-approximate) comoving solver
///  Is necessary for determining the boundary conditions (we require to know which lines overlap to determine those)
// TODO: USING THE SINGLE LINE APPROX, THIS CAN BE SIMPLIFIED
inline void Solver :: get_overlapping_lines(Model& model, const Size nextpoint, bool is_upward_disc)
{
    //For checking whether the lines overlap, we simply check whether the frequency quadrature overlaps
    //Should be conceptually more simple than using the line centers and widths to figure out whether they overlap
    Vector<unsigned char>& line_quad_discdir_overlap=line_quad_discdir_overlap_();
    for (Size lineidx=0; lineidx<model.parameters->nlines(); lineidx++)
    {
        line_quad_discdir_overlap[lineidx]=false;//initialize to no overlap
    }

    Vector<Size>& line_count=line_count_();//counts the number of quadratures encountered for each line

    if (is_upward_disc)
    {
        for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
        {
            //get corresponding line of the sorted frequency index
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpoint, freqidx);//arbitrary frame, as within a single point, one does not need to care about doppler shifts
            const Size lineidx=model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            //increment number of freqs of line quadrature encountered
            line_count[lineidx]++;
            //if we have counted all relevant quads belonging to a specific line:
            if (line_count[lineidx]==model.parameters->nquads())
            {
                line_count[lineidx]=0;//reset the index (as counting this line is no longer needed)
                //Then if no other lines are currently being counted (i.e. entire count vectors is zero),
                // we conclude that the current continguous range has ended
                // TODO SIMPLIFY: LOGIC USING SET AND TOTAL COUNT
                if (!std::all_of(line_count.vec.begin(), line_count.vec.end(), [](int i) { return i==0; }))
                {
                    //just set curr line index to be overlapping (on the right side)
                    line_quad_discdir_overlap[lineidx]=true;
                }
            }
        }
    }
    else
    {   //long for more easily handling reverse loops
        for (long freqidx=model.parameters->nfreqs()-1; freqidx>=0; freqidx--)
        {
            //get corresponding line of the sorted frequency index
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpoint, freqidx);//arbitrary frame, as within a single point, one does not need to care about doppler shifts
            const Size lineidx=model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            //increment number of freqs of line quadrature encountered
            line_count[lineidx]++;
            //if we have counted all relevant quads belonging to a specific line:
            if (line_count[lineidx]==model.parameters->nquads())
            {
                line_count[lineidx]=0;//reset the index (as counting this line is no longer needed)
                //Then if no other lines are currently being counted (i.e. entire count vectors is zero),
                // we conclude that the current continguous range has ended
                if (!std::all_of(line_count.vec.begin(), line_count.vec.end(), [](int i) { return i==0; }))
                {
                    //just set curr line index to be overlapping (on the right side)
                    line_quad_discdir_overlap[lineidx]=true;
                }
            }
        }
    }
    //TODO: check line_count to be uniformly 0!! We did not reinitialize this variable...
}


///  Computes the continguous ranges spanned by the line quadratures (at the current point)
///  Internal function of the (non-approximate) comoving solver
///  In this version, we exclude the 2 farthest frequency points due to boundary conditions
///  Also requires the shift, as we will use these computed ranges to determine overlap with boundary frequencies
//TODO: USING THE SINGLE LINE APPROX, THIS CAN BE SIMPLIFIED ENORMOUSLY (just use all individual line quadratures)
inline void Solver :: get_line_ranges(Model& model, const Size curr_point, bool is_upward_disc, Real curr_shift)
{
    Vector<Real>& left_bound=left_bound_();//Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound=right_bound_();//Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges=nb_ranges_();//contains the number of ranges
    Vector<Size>& left_bound_index=left_bound_index_();//Specifies the corresponding freq index to the left bounds of the ranges
    Vector<Size>& right_bound_index=right_bound_index_();//Specifies the corresponding freq index to the right bounds of the ranges

    nb_ranges=0;//by default, set to 0.
    //First, we specify which quadratues are counted (every one except the outmost 2 ones depending on discretization direction)
    Vector<Size>& quad_range_weight=quad_range_weight_();
    Size& tot_quad_range_weight=tot_quad_range_weight_();
    tot_quad_range_weight=model.parameters->nquads()-2;//Number of used quadratures is the same for both cases, so we do not need to worry about this

    // TODO: maybe replace this implementation with the following:
    // const Size tot_quad_range_weight = model.parameters->nquads()-2;
    // Size tot_quad_count = 0;
    // std::set<Size> encountered_lines = std::set<Size>();
    //... In this, the set might make it a bit more efficient ...

    //We exclude the farthest 2 frequency points for each line, as these correspond to boundary conditions
    if (is_upward_disc)
    {//then everything except the final two quads may be used for determining the range
        Size quadidx=0;
        while (quadidx<model.parameters->nquads()-2)
        {
            quad_range_weight[quadidx]=1;
            quadidx++;
        }

        while (quadidx<model.parameters->nquads())
        {
            quad_range_weight[quadidx]=0;
            quadidx++;
        }
    }
    else
    {//then everything except the first two quads may be used for determining the range
        Size quadidx=0;
        while (quadidx<2)
        {
            quad_range_weight[quadidx]=0;
            quadidx++;
        }

        while (quadidx<model.parameters->nquads())
        {
            quad_range_weight[quadidx]=1;
            quadidx++;
        }
    }

    //Dummy initialization, as we do not know a priori what the first counted quadrature will be
    Real leftbound=0.0;
    Real rightbound=0.0;
    Size leftboundidx=0;
    Size rightboundidx=0;

    bool leftbound_specified=false;//Specifies whether the left bound for the current range has already been set

    Vector<Size>& line_count=line_count_();//counts the number of quadratures encountered for each line
    //ASSUMED AT THIS POINT TO BE INITIALIZED TO 0.

    for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
    {
        //get corresponding line of the sorted frequency index
        const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(curr_point, freqidx);
        const Size lineidx=model.radiation.frequencies.corresponding_line[unsorted_freqidx];
        const Size quadidx=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
        const Real curr_freq=model.radiation.frequencies.sorted_nu(curr_point, freqidx)*curr_shift;//in static frame

        if (!leftbound_specified&&(quad_range_weight[quadidx]))
        {//specify left bound if not yet done so (and the point is not a boundary point)
            leftbound_specified=true;
            leftbound=curr_freq;
            leftboundidx=freqidx;
        }
        rightbound=curr_freq;
        rightboundidx=freqidx;

        // TODO: implementation using a set
        // //better implementation of checking whether we have counted all useful quadrature points of a (maybe overlapping) line region.
        // encountered_lines.insert(lineidx);
        // tot_quad_count+=quad_range_weight[quadidx];//the boundary points need not be counted,
        // //err, I accidentally forgot that we do not increment the tot_quad_count when approaching bdy points (adding nonspecified bounds makes no sense)
        // if ((leftbound_specified)&&((encountered_lines.size()*tot_quad_range_weight) == tot_quad_count))
        // {
        //     //add range
        //     left_bound[nb_ranges]=leftbound;
        //     right_bound[nb_ranges]=rightbound;
        //     left_bound_index[nb_ranges]=leftboundidx;
        //     right_bound_index[nb_ranges]=rightboundidx;
        //     nb_ranges++;
        //
        //     //set flag for new left bound
        //     leftbound_specified=false;
        // }


        //increment number of freqs encountered if a non-boundary point is encountered (implicit conversion from bool to int)
        line_count[lineidx]+=quad_range_weight[quadidx];
        //if we have counted all relevant quads belonging to a specific line:
        if (line_count[lineidx]==tot_quad_range_weight)
        {
            line_count[lineidx]=0;
            //Then if no other lines are currently being counted (i.e. entire count is zero),
            // we conclude that the current continguous range has ended
            //TODO: make this thing scalable, using a set for the lines encoutered, and some counter for total counts!
            if (std::all_of(line_count.vec.begin(), line_count.vec.end(), [](int i) { return i==0; }))
            {
                //add range
                left_bound[nb_ranges]=leftbound;
                right_bound[nb_ranges]=rightbound;
                left_bound_index[nb_ranges]=leftboundidx;
                right_bound_index[nb_ranges]=rightboundidx;
                nb_ranges++;

                //set flag for new left bound
                leftbound_specified=false;
            }
        }
    }
    //TODO: check line_count to be uniformly 0!! As we do not reinitialize it anywhere (and we assume it to be 0 when entering this and another function)
}



///  Sets all implicit boundary conditions, using the computed frequency ranges for currpoint and overlaps for nextpoint
///  Internal function for the (non-approximate) comoving solver
///  Uses the line ranges computed by get_line_ranges to figure out which frequencies (at the next point) lie out of bounds
inline void Solver :: set_implicit_boundary_frequencies(Model& model, const Size nextpoint, const Size nextpointonrayindex, Real shift_next,
                std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index, bool is_upward_disc)
{
    Size curr_range_index=0;
    Vector<Real>& left_bound=left_bound_();
    Vector<Real>& right_bound=right_bound_();
    Vector<unsigned char>& line_quad_overlap=line_quad_discdir_overlap_();
    Size nb_ranges=nb_ranges_();
    Real left_bound_freq=left_bound[curr_range_index];
    Real right_bound_freq=right_bound[curr_range_index];

    if (is_upward_disc)
    {
        // The only difference between upward and downward direction lies in the location of the forced boundary conditions
        for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
        {
            Real next_freq=model.radiation.frequencies.sorted_nu(nextpoint, freqidx)*shift_next;//in static frame
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpoint, freqidx);
            const Size lineidx=model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            const Size quadidx=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];

            // Nonoverlapping edges of a line quadrature need boundary conditions
            if (quadidx>=model.parameters->nquads()-2&&!line_quad_overlap[lineidx])
            {
                //Thus we add this frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
                continue;
            }
            //ranges are sorted, just like the freqs; therefore we can relatively efficently compare the non-overlap
            // by just advancing the index for the ranges
            while((right_bound_freq<next_freq) && (curr_range_index<nb_ranges-1))
            {
                curr_range_index++;
                left_bound_freq=left_bound[curr_range_index];
                right_bound_freq=right_bound[curr_range_index];
            }
            //The frequency lies outside the range if smaller than left_bound OR larger than right_bound
            if (next_freq<left_bound_freq || next_freq>right_bound_freq)
            {
                //Thus we add this frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
            }
        }
    }
    else
    {
        for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
        {
            Real next_freq=model.radiation.frequencies.sorted_nu(nextpoint, freqidx)*shift_next;//in static frame
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpoint, freqidx);
            const Size lineidx=model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            const Size quadidx=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            // Nonoverlapping edges of a line quadrature need boundary conditions
            if (quadidx<2&&!line_quad_overlap[lineidx])
            {
                //Thus we add the frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
                continue;
            }
            //ranges are sorted, just like the freqs; therefore we can relatively efficently compare the non-overlap
            while((right_bound_freq<next_freq) && (curr_range_index<nb_ranges-1))
            {
                curr_range_index++;
                left_bound_freq=left_bound[curr_range_index];
                right_bound_freq=right_bound[curr_range_index];
            }
            //The frequency lies outside the range if smaller than left_bound OR larger than right_bound
            if (next_freq<left_bound_freq || next_freq>right_bound_freq)
            {
                //Thus we add the frequency as a boundary condition
                multimap_freq_to_bdy_index.emplace(next_freq, std::make_tuple(nextpointonrayindex, freqidx));
            }
        }
    }
}


///  Matches the overlapping boundary conditions, using the currently computed ranges
///  Internal function for the (non-approximate) comoving solver
inline void Solver :: match_overlapping_boundary_conditions(Model& model, const Size currpoint, const Size curr_point_on_ray_index,
                                                            const Real curr_shift, std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index)
{
    Vector<Real>& left_bound=left_bound_();//Specifies the left bounds of the ranges in [Hz]
    Vector<Real>& right_bound=right_bound_();//Specifies the right bounds of the ranges in [Hz]
    Size& nb_ranges=nb_ranges_();//contains the number of ranges
    Vector<Size>& left_bound_index=left_bound_index_();//Specifies the corresponding freq index to the left bounds of the ranges
    Vector<Size>& right_bound_index=right_bound_index_();//Specifies the corresponding freq index to the right bounds of the ranges

    Size curr_range_index=0;
    Size curr_freq_index=left_bound_index[curr_range_index];//contains the current frequency index
    //Assumption: at least a single line exists -> a single range exists
    Real curr_range_min=left_bound[curr_range_index];
    Real curr_range_max=right_bound[curr_range_index];
    Real curr_range_freq=curr_range_min;

    //By just using the builtin lower_bound, one can easily find the boundary conditions overlapping with the current ranges
    std::multimap<Real, std::tuple<Size, Size>>::iterator it=multimap_freq_to_bdy_index.lower_bound(curr_range_min);
    //Using a while loop, as we simulataneously iterate over the boundary conditions and over the line ranges
    while (it!=multimap_freq_to_bdy_index.end()&&curr_range_index<nb_ranges)
    {
        const Real bdy_freq=it->first;//get boundary freq
        const Size bdy_point_on_ray_idx=std::get<0>(it->second);//get corresponding point index on ray
        const Size bdy_freq_idx=std::get<1>(it->second);//get corresponding freq index
        //check if bdy freq is smaller than the right bound
        if (bdy_freq>curr_range_max)
        {
            //The boundary frequency lies further than the current line, so use the next line instead
            curr_range_index++;
            if (curr_range_index==nb_ranges)
            {   //no next line exist, so all remaining boundary frequencies lie outside the last range
                return;
            }
            //update info about ranges
            curr_range_min=left_bound[curr_range_index];
            curr_range_max=right_bound[curr_range_index];
            curr_range_freq=curr_range_min;
            curr_freq_index=left_bound_index[curr_range_index];//stores the freq index of curr point
            it=multimap_freq_to_bdy_index.lower_bound(curr_range_min);//and adjust iterator accordingly
            continue;
        }
        //the bdy frequency now lies within the bounds, so just linearly iterate until we find the correct freqs of the ranges (static frame) to interpolate with

        //Find which freq bounds at curr_point exactly correspond to the bdy freq (enclosing it from the right)
        while (bdy_freq>curr_range_freq&&curr_range_freq<curr_range_max)
        {
            curr_freq_index++;
            curr_range_freq=model.radiation.frequencies.sorted_nu(currpoint, curr_freq_index)*curr_shift;//in static frame
        }

        //In this section: I have implemented a few different options for interpolating the boundary conditions

        ///zeroth order interpolation for the boundary frequency
        // Size left_curr_freq_idx=curr_freq_index-1;
        // Size right_curr_freq_idx=curr_freq_index;
        //
        // //the curr range freq should now be larger/equal to the bdy freq
        // //In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as index for interpolation, so use curr_index as left bound instead
        // if (curr_freq_index==left_bound_index[curr_range_index])
        // {
        //   left_curr_freq_idx=curr_freq_index;
        //   right_curr_freq_idx=curr_freq_index+1;
        // }
        //
        // const Real left_freq = model.radiation.frequencies.sorted_nu(currpoint, left_curr_freq_idx)*curr_shift;
        // const Real right_freq = model.radiation.frequencies.sorted_nu(currpoint, right_curr_freq_idx)*curr_shift;
        // Size zeroth_order_freq = 0;
        //
        // //get closest freq
        // if (bdy_freq-left_freq<0.5*(right_freq-left_freq))
        // {
        //     zeroth_order_freq = left_curr_freq_idx;
        // }
        // else
        // {
        //     zeroth_order_freq = right_curr_freq_idx;
        // }

        ///First order interpolation for the boundary frequency
        Size left_curr_freq_idx=curr_freq_index-1;
        Size right_curr_freq_idx=curr_freq_index;

        //the curr range freq should now be larger/equal to the bdy freq
        //In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as index for interpolation, so use curr_index as left bound instead
        if (curr_freq_index==left_bound_index[curr_range_index])
        {
          left_curr_freq_idx=curr_freq_index;
          right_curr_freq_idx=curr_freq_index+1;
        }

        const Real deltafreq=(model.radiation.frequencies.sorted_nu(currpoint, right_curr_freq_idx)
                              -model.radiation.frequencies.sorted_nu(currpoint, left_curr_freq_idx))*curr_shift;//in static frame

        ///Second order interpolation
        // Size left_curr_freq_idx=curr_freq_index-2;
        // Size middle_curr_freq_idx=curr_freq_index-1;
        // Size right_curr_freq_idx=curr_freq_index;
        // //the curr range freq should now be larger/equal to the bdy freq
        // //In the case that bdy_freq lies on the left bound, doing left_bound-1 makes no sense as index for interpolation, so use curr_index as left bound instead
        // if (curr_freq_index==left_bound_index[curr_range_index]||middle_curr_freq_idx==left_bound_index[curr_range_index])
        // {
        //     left_curr_freq_idx=left_bound_index[curr_range_index];
        //     middle_curr_freq_idx=left_bound_index[curr_range_index]+1;
        //     right_curr_freq_idx=left_bound_index[curr_range_index]+2;
        // }
        //
        // //first compute coefficents for the second order accurate freq derivative for the explicit part
        // const Real dfreqsmall=(model.radiation.frequencies.sorted_nu(currpoint, right_curr_freq_idx)-model.radiation.frequencies.sorted_nu(currpoint, middle_curr_freq_idx))
        //             *curr_shift;//in static frame
        // const Real dfreqlarge=(model.radiation.frequencies.sorted_nu(currpoint, right_curr_freq_idx)-model.radiation.frequencies.sorted_nu(currpoint, left_curr_freq_idx))
        //             *curr_shift;//in static frame


        ///Setting the interpolation coefficients for the boundary condition


        //Set frequency derivative correctly

        //Note: the actual factor with which to multiply is Δτ^2/(1-exp(-Δτ)-Δτ*exp(-Δτ)), which in the limit Δτ→0 corresponds to 2 (wait, I did forget a factor exp(-Δτ) in front due to optical depth itself)
        //So in total, the multiplication factor is exp(-Δτ)*Δτ^2/(1-exp(-Δτ)-Δτ*exp(-Δτ))
        //Set explicit coefficents to TWICE the normal 1/Δν value (as we are only treating the explicit part)

        ///zeroth order interpolation
        // //set starting indices correctly for zeroth order interpolation
        // start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[0]=curr_point_on_ray_index;
        // start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1]=zeroth_order_freq;
        //
        // dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        // dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        // dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        // //as the zeroth order does not use the freq derivative, the indices here do not matter
        // dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
        // dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
        // dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;

        //Set starting index correctly (for first, second order interpolation)
        start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[0]=curr_point_on_ray_index;
        start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1]=left_curr_freq_idx;
        ///first order interpolation
        dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=-2.0/deltafreq;
        dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=2.0/deltafreq;
        dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;

        dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
        dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=right_curr_freq_idx;
        dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=right_curr_freq_idx;
        //as the coefficient for the last on is 0, no effort should be made to renumber that last index

        ///second order interpolation (might be worse due to trying to fit some discontinouous thing with a higher order polynomial)
        // dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) =-2.0*dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
        // dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) =2.0*dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
        // dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx) =-dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)-dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx);//curr point itself
        //
        // dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=left_curr_freq_idx;
        // dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=middle_curr_freq_idx;
        // dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=right_curr_freq_idx;

        //Implicit part coefs are set to 0. (as long as the corresponding indicies are in bounds, i do not particularly care about the exact values)
        dIdnu_coef1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        //exact indices do not matter, as long as they are in bounds
        dIdnu_index1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
        dIdnu_index2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
        dIdnu_index3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;

        //set dtau to approx 0
        delta_tau_()(bdy_point_on_ray_idx, bdy_freq_idx)=model.parameters->comoving_min_dtau;//should be small enough, but not small enough to crash my solver (due to /Δτ^2 necessary)
        //The source function is set to 0, as the optical depth increment Δτ≃0 anyway.
        S_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        S_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        //TODO: In case of large doppler shifts, we might need to rethink setting dtau→0; as this will underestimate the optical depth (and influence of source)
        //Concretely, we would need to use the optical depth of the last segment, and somehow incorporate the interpolated intensity of way earlier (and set the source function correctly)

        //pop value from iterator, as we have succesfully set this boundary condition
        it=multimap_freq_to_bdy_index.erase(it);
    }
}


///  For all remaining boundary points, this sets the boundary condition correctly
///  Interal function of the (non-approximate) comoving solver.
///  Uses the initial_bdy (point) for determining which boundary condition to use
inline void Solver :: set_initial_boundary_conditions(Model& model, const Size inital_bdy, const Real curr_shift, std::multimap<Real, std::tuple<Size, Size>>& multimap_freq_to_bdy_index)
{
    std::multimap<Real, std::tuple<Size, Size>>::iterator it=multimap_freq_to_bdy_index.begin();
    while (it!=multimap_freq_to_bdy_index.end())
    {
        const Real bdy_freq=it->first;//get boundary freq (in static frame)
        const Size bdy_point_on_ray_idx=std::get<0>(it->second);//get corresponding point index on ray
        const Size bdy_freq_idx=std::get<1>(it->second);//get corresponding freq index

        //For using somewhat reasonable frequency indices
        Size curr_freq_idx=start_indices_()(bdy_point_on_ray_idx, bdy_freq_idx)[1];

        //Compute boundary intensity
        Real bdy_intensity=boundary_intensity(model, initial_bdy, bdy_freq);

        //Set the frequency derivative coefficients to 0
        dIdnu_coef1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;

        dIdnu_index1_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=curr_freq_idx;
        dIdnu_index2_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=curr_freq_idx;
        dIdnu_index3_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=curr_freq_idx;

        dIdnu_coef1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;
        dIdnu_coef3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=0.0;

        dIdnu_index1_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
        dIdnu_index2_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;
        dIdnu_index3_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_freq_idx;

        // Set S and dtau such that the resulting intensity will be (approximately) equal to the initial boundary intensity
        delta_tau_()(bdy_point_on_ray_idx, bdy_freq_idx)=50.0;//should be large enough
        S_curr_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_intensity;
        S_next_()(bdy_point_on_ray_idx, bdy_freq_idx)=bdy_intensity;

        //TODO: Hmm, if we were to compute S, Δτ (only using the last part), we might be able to get far better results in case of extreme velocity gradients (and limited quadrature points)
        //However, this does not protect against using far too large doppler shifts... (when going back and forth, might jump into a non-computed part; however then +- bdy intensity applies either way?)

        //pop value from iterator, as this boundary condition is now accounted for
        it=multimap_freq_to_bdy_index.erase(it);
    }
}


//Computes the setup and the boundary conditions for the comoving solver for the forward ray direction
//Assumes we have already traced the ray
//TODO: do MORE specialization for the single line approx (taking inspiration from far too complex deque stuff?)
template<ApproximationType approx>
inline void Solver :: comoving_ray_bdy_setup_forward(Model& model, Size last_interesting_rayposidx)
{
    // std::cout<<"entering forward setup"<<std::endl;
    //Note: this assumes forward ray, not backward ray
    // the only practical differences lie in the shift direction and the traversal direction of the ray
    Matrix<Real>& intensities=intensities_();

    Size first_index=first_();//index of first point on ray
    Real shift_first=2.0-shift_()[first_()];
    // Size rayposidx=last_();//ray position index -> point index through nr_()[rayposidx]
    //  -- START FROM ACTUAL INTERESTING POINT -- (tracing everything is a bit too much)
    Size rayposidx=last_interesting_rayposidx;//ray position index -> point index through nr_()[rayposidx]
    // std::cout<<"setting bdy intensities at rayposidx= "<<first_index<<std::endl;
    //TODO: do we want the bdy frequencies in static (is current option) or comoving frame?
    for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
    {
        intensities(first_index, freqid)=boundary_intensity(model, nr_()[first_index], model.radiation.frequencies.sorted_nu(nr_()[first_index], freqid)*shift_first);
    }

    //from first_+1 to last_ index, trace the ray in reverse direction to treat the boundary conditions

    //dummy initialization; just need some references to Reals for get_eta_and_chi
    //Previously defined vectors are not of the right size for me (O(npointsonray) instead of O(nfreqs))//
    //Wrong: eta_c_, ... are the correct dimension, however I see no reason to use them. (as the values will be directly used)
    //So I just compute per point individually (however the compiler should notice and loop unrolling could still happen)
    Real eta_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real chi_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real eta_curr=0.0;
    Real chi_curr=0.0;
    //TODO: use same data struct as current shortchar method; err, might not be possible due to matching arbitrary frequencies with eachother

    //helper values for computing the freq derivative terms
    Real dfreqsmall=0.0;
    Real dfreqlarge=0.0;

    //POSSIBLE OPTIMIZATION: only query them once, moving ..._curr to ..._next after each rayposidx change
    //But then this should change to a vector holding these values...

    std::multimap<Real, std::tuple<Size, Size>> multimap_freq_to_bdy_index;//boundary conditions can overlap, so multimap it is

    bool computing_curr_opacity=true;
    // std::cout<<"starting iteration"<<std::endl;
    //Yes, for determining boundary conditions, we are going backwards; this seems strange, but is a result of easier treatment of the spatial boundary conditions at the end
    while (rayposidx>first_())
    {
        const Size nextpoint=nr_()[rayposidx];
        const Size currpoint=nr_()[rayposidx-1];
        const Real dZ=dZ_()[rayposidx-1];
        //TODO: for less mistakes, redefine rayposidx, rayposidx-1 using different names

        const Real shift_next=2.0-shift_()[rayposidx];//shift is by default defined backwards for the forward ray direction; so under the non-relativistic approx, we can just do 2-shift to get the shift in the other direction
        const Real shift_curr=2.0-shift_()[rayposidx-1];//this is due to its meaning; it is the transformation from the static to the comoving frame of that specific point

        //For the sake of accurately putting boundary conditions, this should point in the correct direction; otherwise we might accidentally put boundary conditions somewhat closer to the middle of a line
        const bool is_upward_disc=(shift_next>=shift_curr);

        //TODO: maybe refactor loop from this point onwards; should be the same for both forward and backward rays

        // std::cout<<"before matching freq indices"<<std::endl;

        //somewhere, we should also compute all of S, Δτ, dI/dν stuff
        //this first part might be appropriate, as we can as well overwrite all this stuff later for the boundary conditions
        //err, Δτ is not available until we know which freq matches which other freq
        //First, we match all frequency indices as good as possible.
        match_frequency_indices(model, nextpoint, currpoint, shift_next, shift_curr, rayposidx, rayposidx-1, is_upward_disc);//O(nfreqs)

        // std::cout<<"before setting default computation values"<<std::endl;

        bool compute_curr_opacity;
        //Now compute all default stuff, computing bogus for the non-matched indices (but as these correspond to boundary indices, we will overwrite this anyway)
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            const Real nextfreq=model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx);//=comoving frame freq at next point
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            const Real currfreq=model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx);//=comoving frame freq at curr point (different frame)
            //technically, we also need to read the point index from this start_indices_ (start_indices_()(rayposidx, next_freq_idx)[0]), but this should still correspond to currpoint at this moment in time
            // const Size nextlineidx=model.radiation.frequencies.corresponding_line_matrix(nextpoint, next_freq_idx);
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpoint, next_freq_idx);
            const Size nextlineidx=model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            //only useful for single line approx

            Real dtau, Snext, Scurr;
            Real chicurr, chinext;//dummy stuff
            //new method for computing S, dtau
            compute_curr_opacity=computing_curr_opacity;
            compute_curr_opacity=true;//due to some shenanigans with the boundary conditions, I'll probably need to compute all sources, dtau's twice
            //OTHERWISE TODO: add some way to remember previous source function, dtau (if boundary condition implementation is flexible enough)
            //needs doppler shift difference, so definition of direction does not matter; err, line index might also not be 100% correct; it can change easily with a large doppler shift...
            //FIXME: shift to same frame! e.g. shift curr point to frame of next point for computing dtau
            // const Real relshift=shift_curr-shift_next;//TODO: check sign; I think it needs to be the default shift computation
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq*relshift, nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            //only using the old method for computing dtau
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, 1.0, 1.0, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            //using only new method for computing dtau
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, 0.0, 1.0, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            //Possible optimization, only query eta/chi_curr once for each point, as the next value is the same either way when traversing through this ray
            // const Size currlineidx=model.radiation.frequencies.corresponding_line_matrix(currpoint, curr_freq_idx);
            // get_eta_and_chi<approx>(model, currpoint, currlineidx, currfreq, eta_curr, chi_curr);
            S_next_()(rayposidx, next_freq_idx)=Snext;
            S_curr_()(rayposidx, next_freq_idx)=Scurr;
            // const Real dtau = std::max(trap(chi_curr, chi_next, dZ), COMOVING_MIN_DTAU);
            // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
            dtau = std::max(dtau, model.parameters->comoving_min_dtau);
            delta_tau_()(rayposidx, next_freq_idx)=dtau;
            // std::cout<<"setting dtau: "<<delta_tau_()(rayposidx, next_freq_idx)<<std::endl;
            // std::cout<<"chi_next: "<<chi_next<<std::endl;
            // std::cout<<"chi_curr: "<<chi_curr<<std::endl;
            // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        }
        //for the next point, we need to know whether to compute
        computing_curr_opacity=compute_curr_opacity;

        // std::cout<<"setting derivative coeffs"<<std::endl;

        //This part precomputes the frequency derivative coefficients (and corresponding indices), conveniently ignoring the fact that for some outermost line quadrature frequencies, we must do more complicated stuff.
        // However, the complicated stuff is handled by the boundary conditions.
        //The discretization is evidently depends on freq disc direction
        if (is_upward_disc)
        {
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            //FIXME: replace with while loop to deal with nfreqs < 3!
            for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs()-2; next_freq_idx++)
            {
                // std::cout<<"get start indices"<<std::endl;
                Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
                // std::cout<<"getting start_indices_"<<std::endl;
                //Modulo operator seems weird, but for the last few freqs, the corresponding freq might be the last one of currpoint == nfreqs()-1
                //If we want to prevent oob accesses and divides by zero, we need to remap the oob indices (my choice is just using modelo)
                //The exact remapped indices should not matter, as all of this will be overwritten when determining the boundary indices
                Size curr_freq_idxp1=(curr_freq_idx+1)%model.parameters->nfreqs();//index+1
                Size curr_freq_idxp2=(curr_freq_idx+2)%model.parameters->nfreqs();//index+2
                // std::cout<<"freq+2: "<<curr_freq_idxp2<<" nfreqs: "<<model.parameters->nfreqs()<<std::endl;
                // std::cout<<"computing curr coeffs"<<std::endl;
                //first compute coefficents for the second order accurate freq derivative for the explicit part
                dfreqsmall=(model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxp1))
                            *shift_curr;
                dfreqlarge=(model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx))
                            *shift_curr;
                // std::cout<<"setting curr coeffs and indices"<<std::endl;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                // std::cout<<"computing next coeffs"<<std::endl;

                //And now do the same (but simpler; less index management) for the implicit part
                dfreqsmall=(model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx+1))
                            *shift_next;
                dfreqlarge=(model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx))
                            *shift_next;

                // std::cout<<"setting next coeffs and indices"<<std::endl;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                // std::cout<<"setting next indices"<<std::endl;
                // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
                // std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx+2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx+1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
        }
        else
        {
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            //Assumes nfreqs>=3
            for (Size next_freq_idx=2; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
            {
                Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
                //Modulo operator seems weird, but for the last few freqs, the corresponding freq might be the last one of currpoint == nfreqs()-1
                //If we want to prevent oob accesses and divides by zero, we need to remap the oob indices (my choice is just using modelo)
                //The exact remapped indices should not matter, as all of this will be overwritten when determining the boundary indices

                //MODULUS FOR 'negative' (more like overflowing) unsigned integer is in general not that easy to compute
                //However, as i know it only goes ever so slightly negative ('-1,-2'): (mod+i)%mod is sufficient for the computation (mod>=1,2)
                Size curr_freq_idxm1=(curr_freq_idx-1+model.parameters->nfreqs())%model.parameters->nfreqs();//index-1
                Size curr_freq_idxm2=(curr_freq_idx-2+model.parameters->nfreqs())%model.parameters->nfreqs();//index-2
                //first compute coefficents for the second order accurate freq derivative for the explicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side
                dfreqsmall=-(model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxm1)-model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxm2))
                *shift_curr;
                dfreqlarge=-(model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx)-model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idxm2))
                *shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side
                dfreqsmall=-(model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx-1)-model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx-2))
                *shift_next;
                dfreqlarge=-(model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx)-model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx-2))
                *shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx-2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx-1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
        }

        //Note to self: we can just use all default freq derivative stuff, except for the 2 outermost points on the grid.
        // Because I define every computation in terms of the next point index, I will have no problems with non-existant/too far frequencies, as I have made the computed curr regions slightly smaller
        //Note: it might seem we will try to use oob values, but these should be overwritten anyway when setting boundary conditions
        //CHECK THIS FOR A SIMPLE EXAMPLE!!!

        // std::cout<<"getting line ranges"<<std::endl;

        //now start classifying the new boundary points
        //For this we first need to compute the ranges of curr_point (minus some boundary freqs)
        get_line_ranges(model, currpoint, is_upward_disc, shift_curr);//for the current point, we definitely need the exact ranges (fiddled with to ensure enough bdy conditions)

        // std::cout<<"getting overlapping lines"<<std::endl;
        // and compute the overlap between lines on next_point to make sure we do not accidentally add extra boundary conditions where lines might overlap
        get_overlapping_lines(model, nextpoint, is_upward_disc);//technically, we should compute which lines overlap only in this part

        // std::cout<<"setting implicit boundary freqs"<<std::endl;
        //First, we set boundary conditions by comparing the frequencies we have at next_point
        // versus the range of frequencies we have a curr_point; any outside that range will be treated as boundary points
        set_implicit_boundary_frequencies(model, nextpoint, rayposidx, shift_next, multimap_freq_to_bdy_index, is_upward_disc);

        // std::cout<<"matching implicit boundary freqs"<<std::endl;
        // And now we check what boundary conditions need to be evaluated using the intensities at curr_point
        match_overlapping_boundary_conditions(model, currpoint, rayposidx-1, shift_curr, multimap_freq_to_bdy_index);
        rayposidx--;
    }
    // std::cout<<"setting initial bounds"<<std::endl;
    //After going through all points, the remaining boundary frequencies have not matched any frequencies near any line,
    // so we will use initial boundary conditions computed at the first point to put on the ray
    //note: currpoint should be point on boundary of domain; in this way, we can simply add the corresponding boundary conditions
    set_initial_boundary_conditions(model, nr_()[first_index], shift_first, multimap_freq_to_bdy_index);
}



//Computes the setup and the boundary conditions for the comoving solver for the backward ray direction
//Assumes we have already traced the ray
//TODO: do MORE specialization for the single line approx (taking inspiration from far too complex deque stuff?)
template<ApproximationType approx>
inline void Solver :: comoving_ray_bdy_setup_backward(Model& model, Size first_interesting_rayposidx)
{
    // std::cout<<"setup backward"<<std::endl;
    //Note: this assumes backward ray, not forward ray
    Matrix<Real>& intensities=intensities_();

    Size first_index=last_();//index of first point on ray
    Real shift_first=shift_()[first_index];//shift on backward ray is reverse of shift of forward ray (which itself is the reverse operation of mapping from comoving to static frame)
    // Size rayposidx=first_();//ray position index -> point index through nr_()[rayposidx]
    Size rayposidx = first_interesting_rayposidx; //as starting from
    for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
    {
        intensities(first_index, freqid)=boundary_intensity(model, nr_()[first_index], model.radiation.frequencies.sorted_nu(nr_()[first_index], freqid)*shift_first);
    }

    //from first_+1 to last_ index, trace the ray in reverse direction to treat the boundary conditions

    //dummy initialization; just need some references to Reals for get_eta_and_chi
    //Previously defined vectors are not of the right size for me (O(npointsonray) instead of O(nfreqs))
    //So I just compute per point individually (however the compiler should notice and loop unrolling could still happen)
    Real eta_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real chi_next=0.0;//get_eta_and_chi<TODO APPROX>(model, nextpoint, ll?, nextfreq, eta, chi);
    Real eta_curr=0.0;
    Real chi_curr=0.0;

    //helper values for computing the freq derivative terms
    Real dfreqsmall=0.0;
    Real dfreqlarge=0.0;

    //POSSIBLE OPTIMIZATION: only query them once, moving ..._curr to ..._next after each rayposidx change
    //But then this should change to a vector holding these values...

    bool computing_curr_opacity=true;
    std::multimap<Real, std::tuple<Size, Size>> multimap_freq_to_bdy_index;//boundary conditions can overlap, so multimap it is
    //Yes, for determining boundary conditions, we are going backwards; this seems strange, but is a result of easier treatment of the spatial boundary conditions at the end
    while (rayposidx<last_())
    {
        const Size nextpoint=nr_()[rayposidx];
        const Size currpoint=nr_()[rayposidx+1];
        const Real dZ=dZ_()[rayposidx];//dZ is stored somewhat finnicky, in locations [first_(),last_()[ (so not including last_())
        //TODO: for less mistakes, redefine rayposidx, rayposidx+1 using different names

        const Real shift_next=shift_()[rayposidx];//shift is by default defined backwards for the forward ray direction; so under the non-relativistic approx, we can just do 2-shift to get the shift in the other direction
        const Real shift_curr=shift_()[rayposidx+1];

        //For the sake of accurately putting boundary conditions, this should point in the correct direction; otherwise we might accidentally put boundary conditions somewhat closer to the middle of a line
        const bool is_upward_disc=(shift_next>=shift_curr);

        //TODO: maybe refactor loop from this point onwards; should be the same for both forward and backward rays

        //somewhere, we should also compute all of S, Δτ, dI/dν stuff
        //this first part might be appropriate, as we can as well overwrite all this stuff later for the boundary conditions
        //err, Δτ is not available until we know which freq matches which other freq
        //First, we match all frequency indices as good as possible.
        match_frequency_indices(model, nextpoint, currpoint, shift_next, shift_curr, rayposidx, rayposidx+1, is_upward_disc);//O(nfreqs)

        bool compute_curr_opacity;
        //Now compute all default stuff, computing bogus for the non-matched indices (but as these correspond to boundary indices, we will overwrite this anyway)
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            const Real nextfreq=model.radiation.frequencies.sorted_nu(nextpoint, next_freq_idx);//=comoving frame freq
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            const Real currfreq=model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx);//
            //technically, we also need to read the point index from this start_indices_ (start_indices_()(rayposidx, next_freq_idx)[0]), but this should still correspond to currpoint at this moment in time
            // const Size nextlineidx=model.radiation.frequencies.corresponding_line_matrix(nextpoint, next_freq_idx);
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpoint, next_freq_idx);
            const Size nextlineidx=model.radiation.frequencies.corresponding_line[unsorted_freqidx];
            Real dtau, Snext, Scurr;
            Real chicurr, chinext;//dummy stuff
            //new method for computing S, dtau
            compute_curr_opacity=computing_curr_opacity;
            compute_curr_opacity=true;//due to some shenanigans with the boundary conditions, I'll probably need to compute all sources, dtau's twice
            //OTHERWISE TODO: add some way to remember previous source function, dtau (if boundary condition implementation is flexible enough)
            //needs doppler shift difference, so definition of direction does not matter; err, line index might also not be 100% correct; it can change easily with a large doppler shift...
            //FIXME: shift to same frame! e.g. shift curr point to frame of next point for computing dtau
            // const Real relshift=shift_curr-shift_next;//TODO: check sign; I think it needs to be the default shift computation
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq*relshift, nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            //only using the old method for computing dtau
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, 1.0, 1.0, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            //using only new method for computing dtau
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, 0.0, 1.0, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            // compute_source_dtau<approx>(model, currpoint, nextpoint, nextlineidx, currfreq, nextfreq, shift_curr, shift_next, dZ, compute_curr_opacity, dtau, chicurr, chinext, Scurr, Snext);
            //Possible optimization, only query eta/chi_curr once for each point, as the next value is the same either way when traversing through this ray
            // const Size currlineidx=model.radiation.frequencies.corresponding_line_matrix(currpoint, curr_freq_idx);
            // get_eta_and_chi<approx>(model, currpoint, currlineidx, currfreq, eta_curr, chi_curr);
            // const Real Snext=eta_next/chi_next;
            // const Real Scurr=eta_curr/chi_curr;
            S_next_()(rayposidx, next_freq_idx)=Snext;
            S_curr_()(rayposidx, next_freq_idx)=Scurr;
            // std::cout<<"setting curr source: "<<Scurr<<std::endl;
            // const Real dtau = std::max(trap (chi_curr, chi_next, dZ), COMOVING_MIN_DTAU);
            // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
            dtau = std::max(dtau, model.parameters->comoving_min_dtau);
            delta_tau_()(rayposidx, next_freq_idx)=dtau;
            // std::cout<<"setting dtau: "<<delta_tau_()(rayposidx, next_freq_idx)<<std::endl;
            // std::cout<<"chi_next: "<<chi_next<<std::endl;
            // std::cout<<"chi_curr: "<<chi_curr<<std::endl;
            // std::cout<<"dZ: "<<dZ<<std::endl;
            // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        }
        computing_curr_opacity=compute_curr_opacity;

        //This part precomputes the frequency derivative coefficients (and corresponding indices), conveniently ignoring the fact that for some outermost line quadrature frequencies, we must do more complicated stuff.
        // However, the complicated stuff is handled by the boundary conditions.
        //The discretization is evidently depends on freq disc direction
        if (is_upward_disc)
        {
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs()-2; next_freq_idx++)
            {
                Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
                //Modulo operator seems weird, but for the last few freqs, the corresponding freq might be the last one of currpoint == nfreqs()-1
                //If we want to prevent oob accesses and divides by zero, we need to remap the oob indices (my choice is just using modelo)
                //The exact remapped indices should not matter, as all of this will be overwritten when determining the boundary indices
                Size curr_freq_idxp1=(curr_freq_idx+1)%model.parameters->nfreqs();//index+1
                Size curr_freq_idxp2=(curr_freq_idx+2)%model.parameters->nfreqs();//index+2
                //first compute coefficents for the second order accurate freq derivative for the explicit part
                dfreqsmall=(model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idxp1))
                            *shift_curr;
                dfreqlarge=(model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idx))
                            *shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxp1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                dfreqsmall=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx+1))
                            *shift_next;
                dfreqlarge=(model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx))
                            *shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx+2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx+1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
        }
        else
        {
            //Do not forget to exclude the last two frequencies, as these cannot have any freq der computation (should definitely bdy conditions)
            for (Size next_freq_idx=2; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
            {
                Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
                //Modulo operator seems weird, but for the last few freqs, the corresponding freq might be the last one of currpoint == nfreqs()-1
                //If we want to prevent oob accesses and divides by zero, we need to remap the oob indices (my choice is just using modelo)
                //The exact remapped indices should not matter, as all of this will be overwritten when determining the boundary indices

                //MODULUS FOR 'negative' (more like overflowing) unsigned integer is in general not that easy to compute
                //However, as i know it only goes ever so slightly negative ('-1,-2'): (mod+i)%mod is sufficient for the computation (mod>=1,2)
                Size curr_freq_idxm1=(curr_freq_idx-1+model.parameters->nfreqs())%model.parameters->nfreqs();//index-1
                Size curr_freq_idxm2=(curr_freq_idx-2+model.parameters->nfreqs())%model.parameters->nfreqs();//index-2
                //first compute coefficents for the second order accurate freq derivative for the explicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side (switched order of term, making sign clearer; just compare to upward disc version)
                dfreqsmall=-(model.radiation.frequencies.nu(currpoint, curr_freq_idxm1)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                *shift_curr;
                dfreqlarge=-(model.radiation.frequencies.nu(currpoint, curr_freq_idx)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
                *shift_curr;
                dIdnu_coef3_curr_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_curr_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_curr_()(rayposidx, next_freq_idx) =-dIdnu_coef3_curr_()(rayposidx, next_freq_idx)-dIdnu_coef2_curr_()(rayposidx, next_freq_idx);//curr point itself
                dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm2;
                dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm1;
                dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

                //And now do the same (but simpler; less index management) for the implicit part
                //Note: minus signs in front, as we are now computing the first derivative using points on the other side
                dfreqsmall=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx-1)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next;
                dfreqlarge=-(model.radiation.frequencies.nu(nextpoint, next_freq_idx)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
                *shift_next;

                dIdnu_coef3_next_()(rayposidx, next_freq_idx) =-dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
                dIdnu_coef2_next_()(rayposidx, next_freq_idx) =dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
                dIdnu_coef1_next_()(rayposidx, next_freq_idx) =-dIdnu_coef3_next_()(rayposidx, next_freq_idx)-dIdnu_coef2_next_()(rayposidx, next_freq_idx);//curr point itself;
                dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx-2;
                dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx-1;
                dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

            }
        }

        //Note to self: we can just use all default freq derivative stuff, except for the 2 outermost points on the grid.
        // Because I define every computation in terms of the next point index, I will have no problems with non-existant/too far frequencies, as I have made the computed curr regions slightly smaller
        //Note: it might seem we will try to use oob values, but these should be overwritten anyway
        //CHECK THIS FOR A SIMPLE EXAMPLE!!!

        //now start classifying the new boundary points
        //For this we first need to compute the ranges of curr_point (minus some boundary freqs)
        get_line_ranges(model, currpoint, is_upward_disc, shift_curr);//for the current point, we definitely need the exact ranges (fiddled with to ensure enough bdy conditions)

        // std::cout<<"getting overlapping lines"<<std::endl;
        // and compute the overlap between lines on next_point
        get_overlapping_lines(model, nextpoint, is_upward_disc);//technically, we should compute which lines overlap only in this part

        // std::cout<<"setting implicit boundary freqs"<<std::endl;
        //First, we set boundary conditions by comparing the frequencies we have at next_point
        // versus the range of frequencies we have a curr_point; any outside that range will be treated as boundary points
        set_implicit_boundary_frequencies(model, nextpoint, rayposidx, shift_next, multimap_freq_to_bdy_index, is_upward_disc);
        // And now we check what boundary conditions need to be evaluated using the intensities at curr_point
        // std::cout<<"matching implicit boundary freqs"<<std::endl;
        match_overlapping_boundary_conditions(model, currpoint, rayposidx+1, shift_curr, multimap_freq_to_bdy_index);
        rayposidx++;
    }
    //After going through all points, the remaining boundary frequencies have not matched any frequencies near any line,
    // so we will use initial boundary conditions computed at the first point to put on the ray
    //note: nr_()[first_index] should be point on boundary of domain; in this way, we can simply add the corresponding boundary conditions
    set_initial_boundary_conditions(model, nr_()[first_index], shift_first, multimap_freq_to_bdy_index);
}



// Note: single line approximation might be useful for determining the splits. IMPLEMENT FANCY VERSION
//Note: first order accurate in space, second order accurate in frequency
template<ApproximationType approx>
inline void Solver :: solve_comoving_order_2_sparse (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        accelerated_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J(o,k) = 0.0;
            }
        })
    }

    std::cout<<"hnrays: "<<model.parameters->hnrays()<<std::endl;
    std::cout<<"number threads: "<<paracabs::multi_threading::n_threads ()<<std::endl;
    // For each ray, solve the radiative transfer equation in a comoving manner
    accelerated_for (rr, model.parameters->hnrays(),
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        // const Real     wt = model.geometry.rays.weight   [rr] * two;
        // const Vector3D nn = model.geometry.rays.direction[rr];

        std::cout << "--- rr = " << rr << std::endl;
        // std::cout<<"number threads: "<<paracabs::multi_threading::n_threads ()<<std::endl;

        // std::cout<<"thread id: "<<paracabs::multi_threading::thread_id()<<std::endl;

        //for every ray to trace
        const Size n_rays_to_trace=points_to_trace_ray_through[rr].size();
        for (Size rayidx=0; rayidx<n_rays_to_trace; rayidx++)
        {
        // accelerated_for (rayidx, n_rays_to_trace,
        // {
            //get corresponding origins of the rays
            const Size o = points_to_trace_ray_through[rr][rayidx];
            const double dshift_max = get_dshift_max (model, o);
            // std::cout<<"here"<<std::endl;
            //FIXME: using dshift_max computed from a single point only works when the line width does not change too much; however, all other solvers are already using this bad approximation
            //trace and solve over the ray in both directions, incrementing J as necessary
            solve_comoving_order_2_sparse<approx>(model, o, rr, rayidx, dshift_max);//solves both for forward and backward ray
            //complicated to decouple ray tracing from solving, so more logic is moved inside this functions
        // })
        }
    })
}



//As stepping a single step should be exactly the same in both directions on the ray, we might as well refactor it out
// //rayposidx==rayposidx_curr_point+-1
//Ergo rayposidx stand for the ray position index of the next point
//rr direction index necessary for determining whether to add intensity to J//TODO? replace with bool denoting whether to add it (precompute if clause somewhere else)
//TODO: figure out whether rayposidx_currpoint is needed!! As it should be replaced with start_indices_...
inline void Solver :: solve_comoving_single_step (Model& model, const Size rayposidx, const Size rayidx, const Size rr, const bool is_upward_disc, const bool forward_ray)
{
    Vector<Size>& nr=nr_();//stores the exact point indices
    // Vector<unsigned char>& real_pt=real_pt_();//stores whether each point is real point (not added extra for interpolation purposes)
    Vector<double>& shift=shift_();//stores the shifts versus the static frame; warning: for forward rays, this should be changed to 2.0-shift; for backward rays, it is correct
    Matrix<Vector<Size>>& start_indices=start_indices_();//for getting the previous indices associated with computing the intensity for the next point
    Matrix<Real>& delta_tau=delta_tau_();//stores the optical depths
    Matrix<Real>& S_curr=S_curr_();
    Matrix<Real>& S_next=S_next_();
    Matrix<Real>& intensities=intensities_();

    Matrix<Real>& dIdnu_coef1_curr=dIdnu_coef1_curr_();//coefficient of the freq point itself
    Matrix<Real>& dIdnu_coef2_curr=dIdnu_coef2_curr_();
    Matrix<Real>& dIdnu_coef3_curr=dIdnu_coef3_curr_();

    Matrix<Size>& dIdnu_index1_curr=dIdnu_index1_curr_();
    Matrix<Size>& dIdnu_index2_curr=dIdnu_index2_curr_();
    Matrix<Size>& dIdnu_index3_curr=dIdnu_index3_curr_();

    Matrix<Real>& dIdnu_coef1_next=dIdnu_coef1_next_();
    Matrix<Real>& dIdnu_coef2_next=dIdnu_coef2_next_();
    Matrix<Real>& dIdnu_coef3_next=dIdnu_coef3_next_();

    Matrix<Size>& dIdnu_index1_next=dIdnu_index1_next_();
    Matrix<Size>& dIdnu_index2_next=dIdnu_index2_next_();
    Matrix<Size>& dIdnu_index3_next=dIdnu_index3_next_();

    // std::cout<<"after renaming all vars"<<std::endl;
    // std::cout<<"rayposidx=next point on ray idx= "<<rayposidx<<std::endl;

    //rayposidx represents the next position on the ray
    const Size nextpointidx=nr[rayposidx];
    const Real next_shift=(forward_ray) ? 2.0-shift[rayposidx] : shift[rayposidx];//absolute shifts, versus static frame!

    // std::cout<<"after getting shift"<<std::endl;

    //Do the explicit part
    for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
    {
        // std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
        //Get the point indices
        const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx)[0];
        // std::cout<<"curr_point_on_ray_index: "<<curr_point_on_ray_index<<std::endl;
        const Size curr_point_idx=nr[curr_point_on_ray_index];
        const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
        const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame

        // std::cout<<"after computing shift"<<std::endl;

        const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
        const Real dtau=delta_tau(rayposidx, next_freq_idx);
        // std::cout<<"dtau: "<<dtau<<std::endl;
        // std::cout<<"rayposidx: "<<rayposidx<<std::endl;

        // std::cout<<"after computing dtau"<<std::endl;

        const Real expl_term=(-expm1(-dtau)-dtau*exp(-dtau))/dtau;
        // std::cout<<"expl term: "<<expl_term<<std::endl;
        // std::cout<<"curr source: "<<S_curr(rayposidx, next_freq_idx)<<std::endl;

        const Real expl_freq_der=dIdnu_coef1_curr(rayposidx, next_freq_idx)*intensities(curr_point_on_ray_index, dIdnu_index1_curr(rayposidx, next_freq_idx))
        +dIdnu_coef2_curr(rayposidx, next_freq_idx)*intensities(curr_point_on_ray_index, dIdnu_index2_curr(rayposidx, next_freq_idx))
        +dIdnu_coef3_curr(rayposidx, next_freq_idx)*intensities(curr_point_on_ray_index, dIdnu_index3_curr(rayposidx, next_freq_idx));
        // std::cout<<"curr freq der coef 1: "<<dIdnu_coef1_curr(rayposidx, next_freq_idx)<<std::endl;
        // std::cout<<"curr freq der coef 2: "<<dIdnu_coef2_curr(rayposidx, next_freq_idx)<<std::endl;
        // std::cout<<"curr freq der coef 3: "<<dIdnu_coef3_curr(rayposidx, next_freq_idx)<<std::endl;
        // std::cout<<"curr_point_on_ray_index: "<<curr_point_on_ray_index<<std::endl;
        // std::cout<<"respective freq index: "<<dIdnu_index1_curr(rayposidx, next_freq_idx)<<std::endl;
        // std::cout<<"intensity coef 1: "<<intensities(curr_point_on_ray_index, dIdnu_index1_curr(rayposidx, next_freq_idx))<<std::endl;
        // std::cout<<"intensity coef 2: "<<intensities(curr_point_on_ray_index, dIdnu_index2_curr(rayposidx, next_freq_idx))<<std::endl;
        // std::cout<<"intensity coef 3: "<<intensities(curr_point_on_ray_index, dIdnu_index3_curr(rayposidx, next_freq_idx))<<std::endl;
        // std::cout<<"expl_freq_der: "<<expl_freq_der<<std::endl;
        // std::cout<<"full freq der term: "<<expl_term*expl_freq_der*deltanu/dtau<<std::endl;
        //TODO: actually do the explicit part
        // std::cout<<"intensities before explicit part: "<<intensities(curr_point_on_ray_index, curr_freq_idx)<<std::endl;
        // std::cout<<"intensities after optical depth: "<<intensities(curr_point_on_ray_index, curr_freq_idx)*exp(-dtau)<<std::endl;
        intensities(rayposidx, next_freq_idx)=intensities(curr_point_on_ray_index, curr_freq_idx)*exp(-dtau)
        +expl_term*S_curr(rayposidx, next_freq_idx) //source term
        +expl_term*expl_freq_der*deltanu/dtau;//frequency derivative term
        // std::cout<<"explicit part intensities: "<<intensities(rayposidx, next_freq_idx)<<std::endl;
    }

    // std::cout<<"after explicit part"<<std::endl;

    //Implicit part ordering depends on the discretization direction
    if (is_upward_disc)
    {
        // std::cout<<"is_upward_disc"<<std::endl;
        //in case of positive doppler shifts, we need to start from the largest frequencies (where we have extra boundary conditions on the end)
        //reverse loop over unsigned int, so indices +1
        for (Size next_freq_idx=model.parameters->nfreqs(); next_freq_idx>0; next_freq_idx--)
        {
            //TODO: actually do the implicit part
            const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx-1)[0];
            const Size curr_point_idx=nr[curr_point_on_ray_index];
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx-1)[1];
            const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx-1)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real dtau=delta_tau(rayposidx, next_freq_idx-1);
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"next freq idx: "<<next_freq_idx-1<<std::endl;
            // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
            const Real impl_term=(dtau+expm1(-dtau))/dtau;
            // std::cout<<"impl term: "<<impl_term<<std::endl;
            // std::cout<<"next source: "<<S_next(rayposidx, next_freq_idx-1)<<std::endl;
            //only the parts of the other points (to subtract/add)
            const Real impl_freq_der=dIdnu_coef2_next(rayposidx, next_freq_idx-1)*intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx-1))
            +dIdnu_coef3_next(rayposidx, next_freq_idx-1)*intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx-1));
            // std::cout<<"impl freq der: "<<impl_freq_der<<std::endl;

            // std::cout<<"intensities before implicit part: "<<intensities(rayposidx, next_freq_idx-1)<<std::endl;

            intensities(rayposidx, next_freq_idx-1)=(intensities(rayposidx, next_freq_idx-1)
            +impl_term*S_next(rayposidx, next_freq_idx-1)//source term
            +impl_term*impl_freq_der*deltanu/dtau)/(1.0-dIdnu_coef1_next(rayposidx, next_freq_idx-1)*impl_term*deltanu/dtau); //freq derivative term
            // std::cout<<"implicit part intensities: "<<intensities(rayposidx, next_freq_idx-1)<<std::endl;
        }
    }
    else
    {
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            //TODO: actually do the implicit part
            const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx)[0];
            const Size curr_point_idx=nr[curr_point_on_ray_index];
            const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real dtau=delta_tau(rayposidx, next_freq_idx);
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"next freq idx: "<<next_freq_idx<<std::endl;
            const Real impl_term=(dtau+expm1(-dtau))/dtau;
            // std::cout<<"impl_term: "<<impl_term<<std::endl;
            //only the parts of the other points (to subtract/add)
            const Real impl_freq_der=dIdnu_coef2_next(rayposidx, next_freq_idx)*intensities(rayposidx, dIdnu_index2_next(rayposidx, next_freq_idx))
            +dIdnu_coef3_next(rayposidx, next_freq_idx)*intensities(rayposidx, dIdnu_index3_next(rayposidx, next_freq_idx));
            // std::cout<<"impl freq der: "<<impl_freq_der<<std::endl;
            // std::cout<<"full freq der term: "<<impl_term*impl_freq_der*deltanu/dtau<<std::endl;
            // std::cout<<"1.0-dIdnu coef times...: "<<1.0-dIdnu_coef1_next(rayposidx, next_freq_idx)*deltanu/dtau<<std::endl;
            //
            // std::cout<<"intensities before implicit part: "<<intensities(rayposidx, next_freq_idx)<<std::endl;

            intensities(rayposidx, next_freq_idx)=(intensities(rayposidx, next_freq_idx)
            +impl_term*S_next(rayposidx, next_freq_idx)//source term
            +impl_term*impl_freq_der*deltanu/dtau)/(1.0-dIdnu_coef1_next(rayposidx, next_freq_idx)*impl_term*deltanu/dtau); //freq derivative term
            // std::cout<<"implicit part intensities: "<<intensities(rayposidx, next_freq_idx)<<std::endl;
        }
    }

    // std::cout<<"incrementing J"<<std::endl;

    //Finally increment J if and real point check if closest ray
    // if (real_pt[rayposidx]&&closest_ray(rr, nextpointidx)==rayidx)
    if (closest_ray(rr, nextpointidx)==rayidx)
    {
        // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        // std::cout<<"adding to J at rayposidx: "<<rayposidx<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< quadrature number corresponding to frequency
            // TODO: replace with matrix variant when replacing these Vectors
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpointidx, freqid);
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            // std::cout<<"freqid: "<<freqid<<std::endl;
            // std::cout<<"adding: "<<lspec.quadrature.weights[z] * wt * intensities(rayposidx, freqid)<<std::endl;
            // std::cout<<"intensity: "<<intensities(rayposidx, freqid)<<std::endl;
            lspec.J(nextpointidx,k) += lspec.quadrature.weights[z] * wt * intensities(rayposidx, freqid);// Su_()[centre];

            //Computing the ALI element

            //TODO: check whether this is correct (print out Jeff, Jlin in optically thick model)
            const Size curr_point_on_ray_index=start_indices_()(rayposidx, freqid)[0];
            const Size curr_point_idx=nr[curr_point_on_ray_index];
            const Size curr_freq_idx=start_indices_()(rayposidx, freqid)[1];
            const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            const Real deltanu=model.radiation.frequencies.sorted_nu(nextpointidx, freqid)*next_shift-model.radiation.frequencies.sorted_nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real dtau=delta_tau(rayposidx, freqid);
            const Real constant = lspec.linedata.A[k] * lspec.quadrature.weights[z] * model.geometry.rays.weight[rr];//for integrating over both frequency (weighted with profile function) and angle. We also need the einstein A coefficient

            const Real source_term=(dtau+expm1(-dtau));// '/dtau' has already been applied to simplify the fraction
            const Real lambdaterm=constant*source_term/(dtau-dIdnu_coef1_next(rayposidx, freqid)*source_term*deltanu);

            // frq = freqs.nu(nr[n], f) * shift[n];
            // phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
            // L   = constante * frq * phi * L_lower(m,n) * inverse_chi[n];

            lspec.lambda.add_element(nextpointidx, k, nextpointidx, lambdaterm);
        }
    }
}



template<ApproximationType approx>
inline void Solver :: solve_comoving_order_2_sparse (
      Model& model,
      const Size o,//ray origin point
      const Size r,//ray direction index
      const Size rayidx,//ray trace index
      const double dshift_max)
{
    const Size rr=r;
    const Size ar=model.geometry.rays.antipod[rr];

    // std::cout<<"solving for point: "<<o<<std::endl;
    // std::cout<<"solving for rayindex: "<<rayidx<<std::endl;

    // Vector<Real>& eta_c = eta_c_();//current emissivity for all freqs
    // Vector<Real>& eta_n = eta_n_();//next emissivity    for all freqs
    //
    // Vector<Real>& chi_c = chi_c_();//current opacity    for all freqs
    // Vector<Real>& chi_n = chi_n_();//next opacity       for all freqs

    // Vector<Real>& curr_intensity = curr_intensity_();//current intensity for all freqs
    // Vector<Real>& next_intensity = next_intensity_();//storage for next intensity for all freqs

    Matrix<Real>& intensities=intensities_();//will contain all computed intensities

    // Vector<unsigned char>& real_pt= real_pt_();//for denoting whether the given point is a real point
    Vector<Size>& nr=nr_();

    Vector<double>& shift = shift_();// contains the doppler shifts for all points on the ray
    //for the forward ray, do 2.0-shift to get the shift direction I want for my frequencies
    //This is due to the default computed shift mapping from static frame to comoving frame (I want to do the reverse)
    //for the backward ray, this should be correct (as we need to do the reverse)

    // const Real     wt = model.geometry.rays.weight   [rr];//ray direction weight

    double Z = 0.0; //distance along ray
    double dZ= 0.0; //last distance increment

    // std::cout<<"before tracing rays"<<std::endl;

    // //trace the ray, getting all points on the ray and indicating whether they are real
    // first_() = trace_ray_indicate_point(model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    // last_ () = trace_ray_indicate_point(model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;


    Size first_interesting_rayposidx = centre;
    Size last_interesting_rayposidx = centre;

    //trace the ray, getting all points on the ray and indicating whether they are real
    //no more interpolations, so no more checking which points are real
    //Note: technically, it does not matter which frame I use, as long as it is a fixed frame for the entire ray
    first_() = trace_ray_comoving<Rest>(model.geometry, o, rr, rr, rayidx, dshift_max, -1, centre-1, centre-1, first_interesting_rayposidx) + 1;
    last_ () = trace_ray_comoving<Rest>(model.geometry, o, ar, rr, rayidx, dshift_max, +1, centre+1, centre, last_interesting_rayposidx) - 1;

    //old implementation
    //also set ray start as real point; accidentally forgot
    // real_pt_()[centre]=true;

    nr_   ()[centre] = o;
    shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0);
    // TODO: compute shift versus zero velocity!
    n_tot_() = (last_()+1) - first_();

    // std::cout<<"before forward setup"<<std::endl;

    //Now is the perfect time to setup the boundary conditions and data for the forward ray
    comoving_ray_bdy_setup_forward<approx>(model, last_interesting_rayposidx);

    // std::cout<<"dtau first+1: "<<delta_tau_()(first_()+1, 0)<<std::endl;
    //hmm, seem to be correct for now

    // std::cout<<"after forward setup"<<std::endl;

    //check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the if-clause
    if (closest_ray(rr, nr[first_()])==rayidx)
    {
        // std::cout<<"setting bdy intensities at first"<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nr[first_()], freqid);
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            lspec.J(nr[first_()],k) += lspec.quadrature.weights[z] * wt * intensities(first_(), freqid);// Su_()[centre];
            //TODO: compute lambda term
        }
    }
    // else
    // {
    //     std::cout<<"not setting bdy intensities at first"<<std::endl;
    // }

    // std::cout<<"after forward boundary thing"<<std::endl;

    Size rayposidx=first_()+1;//ray position index -> point index through nr[rayposidx]
    //from first_+1 to last_ index, trace the ray in
    // while (rayposidx<=last_())
    while (rayposidx<=last_interesting_rayposidx)
    {
        const Real shift_next=2.0-shift_()[rayposidx];
        const Real shift_curr=2.0-shift_()[rayposidx-1];
        const bool is_upward_disc=(shift_next>=shift_curr);
        // std::cout<<"solving single step"<<std::endl;
        // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        // std::cout<<"dtau rayposidx: "<<delta_tau_()(rayposidx, 0)<<std::endl;
        solve_comoving_single_step (model, rayposidx, rayidx, r, is_upward_disc, true);
        rayposidx++;
    }

    // std::cout<<"after forward iterations"<<std::endl;

    //Same procedure for the backward ray
    // std::cout<<"shift_()[last]: "<<shift_()[last_()]<<std::endl;

    //Now is the perfect time to setup the boundary conditions and data for the backward ray
    comoving_ray_bdy_setup_backward<approx>(model, first_interesting_rayposidx);

    // std::cout<<"after backwards setup"<<std::endl;

    //check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the if-clause
    if (closest_ray(rr, nr[last_()])==rayidx)
    {
        // std::cout<<"setting bdy intensities at last"<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nr[first_()], freqid);
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            lspec.J(nr[last_()],k) += lspec.quadrature.weights[z] * wt * intensities(last_(), freqid);// Su_()[centre];
            //TODO: compute lambda term
        }
    }


    // std::cout<<"after backwards boundary thing"<<std::endl;

    rayposidx=last_()-1;//ray position index -> point index through nr[rayposidx]
    //from last_()-1 to first_ index, trace the ray backwards; Warning: rayposidx is an unsigned int, so +1 necessary to both sides
    // while (rayposidx+1>=first_()+1)
    while (rayposidx+1>=first_interesting_rayposidx+1)
    {
        const Real shift_next=shift_()[rayposidx];
        const Real shift_curr=shift_()[rayposidx+1];
        // std::cout<<"curr ray index: "<<rayposidx+1<<std::endl;
        // std::cout<<"next ray index: "<<rayposidx<<std::endl;
        const bool is_upward_disc=(shift_next>=shift_curr);
        solve_comoving_single_step (model, rayposidx, rayidx, r, is_upward_disc, false);
        rayposidx--;
    }

    // std::cout<<"first ray index: "<<first_()<<std::endl;
    // std::cout<<"last ray index: "<<last_()<<std::endl;
    // std::cout<<"total number points: "<<n_tot_()<<std::endl;

    // std::cout<<"after backwards iteration"<<std::endl;
}



// Note: single line approximation might be useful for determining the splits. IMPLEMENT FANCY VERSION
//Note: first order accurate in space, second order accurate in frequency
//Approximates radiative transfer to be local in frequency space (so no long-range line overlap)
template<ApproximationType approx>
inline void Solver :: solve_comoving_local_approx_order_2_sparse (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        accelerated_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J(o,k) = 0.0;
            }
        })
    }

    std::cout<<"hnrays: "<<model.parameters->hnrays()<<std::endl;
    std::cout<<"number threads: "<<paracabs::multi_threading::n_threads ()<<std::endl;
    // For each ray, solve the radiative transfer equation in a comoving manner
    accelerated_for (rr, model.parameters->hnrays(),
    {
        const Size     ar = model.geometry.rays.antipod  [rr];

        std::cout << "--- rr = " << rr << std::endl;
        //for every ray to trace
        const Size n_rays_to_trace=points_to_trace_ray_through[rr].size();
        for (Size rayidx=0; rayidx<n_rays_to_trace; rayidx++)
        {
        // accelerated_for (rayidx, n_rays_to_trace,
        // {
            //get corresponding origins of the rays
            const Size o = points_to_trace_ray_through[rr][rayidx];
            const double dshift_max = get_dshift_max (model, o);
            solve_comoving_local_approx_order_2_sparse<approx>(model, o, rr, rayidx, dshift_max);//solves both for forward and backward ray
            //complicated to decouple ray tracing from solving, so more logic is moved inside this functions
        // })
        }
    })
}


template<ApproximationType approx>
accel inline void Solver :: solve_comoving_local_approx_order_2_sparse (
      Model& model,
      const Size o,//ray origin point
      const Size r,//ray direction index
      const Size rayidx,//ray trace index
      const double dshift_max)
{
    // std::cout<<"init on ray"<<std::endl;
    const Size rr=r;
    const Size ar=model.geometry.rays.antipod[rr];

    // std::cout<<"solving for point: "<<o<<std::endl;
    // std::cout<<"solving for rayindex: "<<rayidx<<std::endl;

    // Vector<Real>& eta_c = eta_c_();//current emissivity for all freqs
    // Vector<Real>& eta_n = eta_n_();//next emissivity    for all freqs
    //
    // Vector<Real>& chi_c = chi_c_();//current opacity    for all freqs
    // Vector<Real>& chi_n = chi_n_();//next opacity       for all freqs

    // Vector<Real>& curr_intensity = curr_intensity_();//current intensity for all freqs
    // Vector<Real>& next_intensity = next_intensity_();//storage for next intensity for all freqs

    Vector<Real>& cma_intensities=cma_computed_intensities_();//will computed intensities for a single step
    // Vector<Real>& cma_curr_frequencies=cma_curr_frequencies_();//will store the frequencies of the current point for computing the explicit frequency derivative

    // Vector<unsigned char>& real_pt= real_pt_();//for denoting whether the given point is a real point
    Vector<Size>& nr=nr_();

    Vector<double>& shift = shift_();// contains the doppler shifts for all points on the ray
    //for the forward ray, do 2.0-shift to get the shift direction I want for my frequencies
    //This is due to the default computed shift mapping from static frame to comoving frame (I want to do the reverse)
    //for the backward ray, this should be correct (as we need to do the reverse)

    // const Real     wt = model.geometry.rays.weight   [rr];//ray direction weight

    double Z = 0.0; //distance along ray
    double dZ= 0.0; //last distance increment

    // std::cout<<"before tracing rays"<<std::endl;

    // //trace the ray, getting all points on the ray and indicating whether they are real
    // first_() = trace_ray_indicate_point(model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    // last_ () = trace_ray_indicate_point(model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;

    //trace the ray, getting all points on the ray and indicating whether they are real
    //no more interpolations, so no more checking which points are real
    //Note: technically, it does not matter which frame I use, as long as it is a fixed frame for the entire ray
    // first_() = trace_ray<Rest>(model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    // last_ () = trace_ray<Rest>(model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
    //We should also check which points are actually necessary to compute
    Size first_interesting_rayposidx = centre;
    Size last_interesting_rayposidx = centre;
    first_() = trace_ray_comoving<Rest>(model.geometry, o, rr, rr, rayidx, dshift_max, -1, centre-1, centre-1, first_interesting_rayposidx) + 1;
    last_ () = trace_ray_comoving<Rest>(model.geometry, o, ar, rr, rayidx, dshift_max, +1, centre+1, centre, last_interesting_rayposidx) - 1;


    // std::cout<<"traced ray"<<std::endl;

    //old implementation
    //also set ray start as real point; accidentally forgot
    // real_pt_()[centre]=true;

    nr_   ()[centre] = o;
    shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0);
    n_tot_() = (last_()+1) - first_();

    //Doppler shifts are computed in the opposite direction as usual, so correcting for this
    for (Size temp_rayidx=first_(); temp_rayidx<=last_(); temp_rayidx++)
    {
        shift_()[temp_rayidx] = 2.0-shift_()[temp_rayidx];
    }

    // std::cout<<"did the shift"<<std::endl;
    //and also initialize helper variable for solver
    Vector<char>& cma_compute_curr_opacity = cma_compute_curr_opacity_();
    ///DEBUG ONLY HALF OF RAY
    // compute_curr_opacity.fill(compute_curr_opacity.begin(), compute_curr_opacity.end(), 1); grr, Vector is not a vector
    //manually fill instead
    for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
    {
        cma_compute_curr_opacity[freqidx]=true;
        // cma_curr_frequencies[freqidx] = model.radiation.frequencies.sorted_nu(nr_()[first_()], freqidx);//in comoving frame
        //also compute starting intensity ;; err, why do I bother with a shift? It should be in the frame of the point itself!!
        // cma_intensities[freqidx]=boundary_intensity(model, nr_()[first_()], model.radiation.frequencies.sorted_nu(nr_()[first_()], freqidx)*shift_()[first_()]);
        cma_intensities[freqidx]=boundary_intensity(model, nr_()[first_()], model.radiation.frequencies.sorted_nu(nr_()[first_()], freqidx));
        // std::cout<<"cma_I: "<<cma_intensities[freqidx]<<std::endl;
    }

    // std::cout<<"boundary intensities computed"<<std::endl;

    //check if closest ray of the starting boundary point // maybe todo: replace with some weights 0/1 for eliminating the if-clause
    if (closest_ray(rr, nr[first_()])==rayidx)
    {
        // std::cout<<"setting bdy intensities at first"<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nr[first_()], freqid);
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            lspec.J(nr[first_()],k) += lspec.quadrature.weights[z] * wt * cma_intensities[freqid];// Su_()[centre];
            //TODO: compute lambda term; err, for bdy? technically zero as bdy condition is enforced, so ignore
        }
    }
    // else
    // {
    //     std::cout<<"not setting bdy intensities at first"<<std::endl;
    // }

    // std::cout<<"after forward boundary thing"<<std::endl;

    Size rayposidx=first_()+1;//ray position index -> point index through nr[rayposidx]
    //from first_+1 to last_ index, trace the ray in

    // std::cout<<"iterating over all point on ray (forward)"<<std::endl;


    // while (rayposidx<=last_())
    //in order to reduce the amount of redundant computations
    while (rayposidx<=last_interesting_rayposidx)
    {
        // const Real shift_next=2.0-shift_()[rayposidx];
        // const Real shift_curr=2.0-shift_()[rayposidx-1];
        // std::cout<<"rayposidx: "<<rayposidx<<std::endl;

        const Real shift_next=shift_()[rayposidx];
        const Real shift_curr=shift_()[rayposidx-1];
        const bool is_upward_disc=(shift_next>=shift_curr);
        // std::cout<<"solving single step"<<std::endl;
        // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        // std::cout<<"dtau rayposidx: "<<delta_tau_()(rayposidx, 0)<<std::endl;
        // comoving_local_approx_map_data(model, rayposidx, rayidx, r, is_upward_disc, true);
        // std::cout<<"mapping data"<<std::endl;
        comoving_local_approx_map_data<approx>(model, nr_()[rayposidx-1], nr_()[rayposidx], shift_curr, shift_next, is_upward_disc, dZ_()[rayposidx-1], nr_()[first_()]);
        // std::cout<<"solving single step"<<std::endl;
        solve_comoving_local_approx_single_step (model, nr_()[rayposidx], rayidx, r, is_upward_disc);
        rayposidx++;
        // std::cout<<"next iteration"<<std::endl;
    }

    // std::cout<<"after forward iterations"<<std::endl;

    //Same procedure for the backward ray
    // std::cout<<"shift_()[last]: "<<shift_()[last_()]<<std::endl;

    //Now is the perfect time to setup the boundary conditions and data for the backward ray
    // comoving_ray_bdy_setup_backward<approx>(model);

    // std::cout<<"after backwards setup"<<std::endl;

    //Doppler shifts in opposite direction evidently need to have the opposite shift
    for (Size temp_rayidx=first_(); temp_rayidx<=last_(); temp_rayidx++)
    {
        shift_()[temp_rayidx] = 2.0-shift_()[temp_rayidx];
    }

    //manually fill compute_curr_opacity, as we have not yet computed the current opacities on this ray direction (dopper shift changing sign means we cannot reuse the last opacity on the opposite direction)
    for (Size freqidx=0; freqidx<model.parameters->nfreqs(); freqidx++)
    {
        cma_compute_curr_opacity[freqidx]=true;
        // cma_curr_frequencies[freqidx] = model.radiation.frequencies.sorted_nu(nr_()[last_()], freqidx);//in comoving frame
        //also compute starting intensity // shift not necessary, as we are interested in the local reference frame intensity
        // cma_intensities[freqidx]=boundary_intensity(model, nr_()[last_()], model.radiation.frequencies.sorted_nu(nr_()[last_()], freqidx)*shift_()[last_()]);
        cma_intensities[freqidx]=boundary_intensity(model, nr_()[last_()], model.radiation.frequencies.sorted_nu(nr_()[last_()], freqidx));
        // std::cout<<"cma_I: "<<cma_intensities[freqidx]<<std::endl;

    }


    //check if closest ray // maybe todo: replace with some weights 0/1 for eliminating the if-clause
    if (closest_ray(rr, nr[last_()])==rayidx)
    {
        // std::cout<<"setting bdy intensities at last"<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< number of line number corresponding to frequency
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nr[last_()], freqid);
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            lspec.J(nr[last_()],k) += lspec.quadrature.weights[z] * wt * cma_intensities[freqid];// Su_()[centre];
            //TODO: compute lambda term; err, for bdy? technically zero as bdy condition is enforced, so ignore
        }
    }


    // std::cout<<"after backwards boundary thing"<<std::endl;
    // ALSO PUT SHIFT HERE IN CORRECT DIRECTION
    rayposidx=last_()-1;//ray position index -> point index through nr[rayposidx]
    //from last_()-1 to first_ index, trace the ray backwards; Warning: rayposidx is an unsigned int, so +1 necessary to both sides
    // while (rayposidx+1>=first_()+1)
    //in order to reduce the amount of redundant computations
    while (rayposidx+1>=first_interesting_rayposidx+1)
    {
        const Real shift_next=shift_()[rayposidx];
        const Real shift_curr=shift_()[rayposidx+1];
        // std::cout<<"curr ray index: "<<rayposidx+1<<std::endl;
        // std::cout<<"next ray index: "<<rayposidx<<std::endl;
        const bool is_upward_disc=(shift_next>=shift_curr);
        // solve_comoving_single_step (model, rayposidx, rayidx, r, is_upward_disc, false);
        comoving_local_approx_map_data<approx>(model, nr_()[rayposidx+1], nr_()[rayposidx], shift_curr, shift_next, is_upward_disc, dZ_()[rayposidx], nr_()[last_()]);
        solve_comoving_local_approx_single_step (model, nr_()[rayposidx], rayidx, r, is_upward_disc);
        rayposidx--;
    }

    // std::cout<<"first ray index: "<<first_()<<std::endl;
    // std::cout<<"last ray index: "<<last_()<<std::endl;
    // std::cout<<"total number points: "<<n_tot_()<<std::endl;

    // std::cout<<"after backwards iteration"<<std::endl;
}



// Matches the sorted frequency indices such that the frequency difference is minimal
// ALSO relies on the current frequency split
// ERR, should just iterate over all frequencies, ignoring those different splits. Boundary conditions
// is_upward_disc denotes whether we want an upward discretization
// nextpointonrayindex should denote how far the next point lies on the ray
// currpointonrayindex should do the same, but for the previous point//, const Size currpointonrayindex, const Size nextpointonrayindex
template<ApproximationType approx>
accel inline void Solver :: comoving_local_approx_map_data(Model& model, const Size curr_point, const Size next_point, const Real curr_shift, const Real next_shift, const bool is_upward_disc, const Real dZ, const Size bdy_point)
{//, const Size splitindex
    //if upward discretization, we start from the uppermost part
    bool is_in_bounds = false;
    // Size last_used_bdy_index = 0;//just for convenience sake when determining the explicit frequency derivative
    Size in_bounds_count = 0;//to decide whether to artificially add some boundary to the implicit part (if < 2) then the implicit freq derivative cannot be computed up to second order
    Size curr_freq_count = 0;//count for explicit part to determine when to say that something is outside of our bounds
    std::set<Size> encountered_lines =std::set<Size>();
    // const Real max_curr_nu_dist = TODO;
    if (is_upward_disc)
    {
        //Starting from the highest frequency
        //+1 to all indices due to using unsigned ints when looping down (overflow otherwise)
        Size curr_freq_idx=model.parameters->nfreqs()-1+1;
        Size next_freq_idx=model.parameters->nfreqs()-1+1;
        // const Real curr_sorted_nu = model.radiation.frequencies.sorted_nu(currpoint, curr_freq_idx-1);
        Size unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(curr_point, curr_freq_idx-1);//arbitrary frame, as within a single point, one does not need to care about doppler shifts
        Size curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
        encountered_lines.insert(curr_line_idx);

        const Size min_freq_idx=0+1;
        //for determining whether the distance is large enough, just check whether the corresponding line changes ERR, this does not work, unless we assume the lines to lie far enough from eachother (utilize in single line approx)
        //assumes at least a single frequency will not be a boundary condition (is reasonable if we limit the doppler shift)
        //NOTE: oob freqs should be dealt with later on
        while ((curr_freq_idx>=(min_freq_idx+1)) && next_freq_idx>=min_freq_idx)
        {
            //check if the frequency at the previous point index is higher than this frequency (in static frame)
            if (model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx-1)*curr_shift>
                model.radiation.frequencies.sorted_nu(next_point, next_freq_idx-1)*next_shift)
            {
                // //An implementation without break statement exists, but then we need to add some extra helper boolean to signal that the next thing is oob.
                // if (curr_freq_idx==min_freq_idx)
                // {
                //     break;
                // }
                //more lazy implementation: put is_in_bounds_computation in last part
                curr_freq_idx--;
                // std::cout<<"curr_freq_idx"<<curr_freq_idx<<" min_freq_idx: "<<min_freq_idx<<std::endl;

                unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(curr_point, curr_freq_idx-1);//arbitrary frame, as within a single point, one does not need to care about doppler shifts
                curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
                encountered_lines.insert(curr_line_idx);
                curr_freq_count++;

                // //branchless version; ignoring resetting curr_freq_count
                // const bool branchless_if = (curr_freq_count == (encountered_lines.size()-1)*model.parameters->nquads());
                // is_in_bounds = !branchless_if;
                // in_bounds_count *= !branchless_if;//reset, as this will be used for determining boundary during mapping

                //whenever we encounter a new line, the size of encountered_lines will increase
                //This condition is only satisfied whenever we encounter a new line after all quads of all previous lines are encountered
                //Thus we define this to be true when encountering a gap in the current (explicit) frequency quadrature
                if (curr_freq_count == (encountered_lines.size()-1)*model.parameters->nquads())
                {
                    is_in_bounds = false;
                    in_bounds_count = 0;
                    //Also reset count for encountered lines
                    // encountered_lines = std::set<Size>{curr_line_idx};
                    // curr_freq_count = 0;
                }
                else
                {
                    is_in_bounds = true;
                }
            }
            else//freq matching index must be determined if it is just higher than the freq at the previous point
            {
                const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx-1)*curr_shift;//index - 1 due to size
                const Real next_freq = model.radiation.frequencies.sorted_nu(next_point, next_freq_idx-1)*next_shift;//index - 1 due to size
                comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ, is_in_bounds, in_bounds_count, next_freq_idx-1, curr_freq_idx-1, next_freq, curr_freq, next_shift, curr_shift, curr_line_idx, is_upward_disc, bdy_point);

                // TODO GO MAP STUFF (and compute S_dtau, ...)
                // freq_matching_()[next_freq_idx-1]=curr_freq_idx-1;
                // start_indices_()(nextpointonrayindex,next_freq_idx-1)[0]=currpointonrayindex;//stores point on ray index
                // start_indices_()(nextpointonrayindex,next_freq_idx-1)[1]=curr_freq_idx-1;//stores freq
                //for the sake of boundary conditions, these start_indices_ might be overwritten later
                next_freq_idx--;
                in_bounds_count++;
                // curr_freq_count++;
            }
        }
        // is_in_bounds=false;

        //match remaining points as best as possible (with the outermost point)
        while (next_freq_idx>=min_freq_idx)
        {
            //TODO REPLACE CURR_FREQ_IDX WITH MIN_FREQ_IDX
            // const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx-1)*curr_shift;//index - 1 due to size
            const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, min_freq_idx-1)*curr_shift;//index - 1 due to size
            const Real next_freq = model.radiation.frequencies.sorted_nu(next_point, next_freq_idx-1)*next_shift;//index - 1 due to size
            is_in_bounds = (next_freq>=curr_freq);
            comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ, is_in_bounds, in_bounds_count, next_freq_idx-1, curr_freq_idx-1, next_freq, curr_freq, next_shift, curr_shift, curr_line_idx, is_upward_disc, bdy_point);

            // start_indices_()(nextpointonrayindex,next_freq_idx-1)[0]=currpointonrayindex;//stores point
            // start_indices_()(nextpointonrayindex,next_freq_idx-1)[1]=min_freq_idx-1;//stores freq
            next_freq_idx--;
            in_bounds_count++;
        }


    }
    else
    {
        //Starting from the lowest frequency
        Size curr_freq_idx = 0;
        Size next_freq_idx = 0;
        const Size max_freq_idx = model.parameters->nfreqs()-1;

        Size unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(curr_point, curr_freq_idx);//arbitrary frame, as within a single point, one does not need to care about doppler shifts
        Size curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
        encountered_lines.insert(curr_line_idx);

        //assumes at least a single frequency will not be a boundary condition (is reasonable if we limit the doppler shift)
        while ((curr_freq_idx<=(max_freq_idx-1)) && next_freq_idx<=max_freq_idx)
        {
            //check if a previous frequency index exists which is higher than this (in static frame)
            if (model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx)*curr_shift<
                model.radiation.frequencies.sorted_nu(next_point, next_freq_idx)*next_shift)
            {
                //An implementation without break statement exists, but then we need to add some extra helper boolean to signal that the next thing is oob.
                // if (curr_freq_idx==max_freq_idx)
                // {
                //     break;
                // }
                // More lazy implemntation
                curr_freq_idx++;
                // std::cout<<"curr_freq_idx: "<<curr_freq_idx<<" next_freq_idx: "<<next_freq_idx<<"max_freq_idx: "<<max_freq_idx<<std::endl;

                unsorted_freqidx = model.radiation.frequencies.corresponding_nu_index(curr_point, curr_freq_idx);//arbitrary frame, as within a single point, one does not need to care about doppler shifts
                curr_line_idx = model.radiation.frequencies.corresponding_line[unsorted_freqidx];
                curr_freq_count++;
                encountered_lines.insert(curr_line_idx);

                // //branchless version; ignoring resetting curr_freq_count
                // const bool branchless_if = (curr_freq_count == (encountered_lines.size()-1)*model.parameters->nquads());
                // is_in_bounds = !branchless_if;
                // in_bounds_count *= !branchless_if;//reset, as this will be used for determining boundary during mapping

                //whenever we encounter a new line, the size of encountered_lines will increase
                //This condition is only satisfied whenever we encounter a new line after all quads of all previous lines are encountered
                //Thus we define this to be true when encountering a gap in the current (explicit) frequency quadrature
                if (curr_freq_count == (encountered_lines.size()-1)*model.parameters->nquads())
                {
                    is_in_bounds = false;
                    in_bounds_count = 0;
                    //Also reset count for encountered lines
                    // encountered_lines = std::set<Size>{curr_line_idx};
                    // curr_freq_count = 0;
                }
                else
                {
                    is_in_bounds = true;
                }
            }
            else//freq matching index must be determined if it is just higher than the freq at the previous point
            {
                const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx)*curr_shift;
                const Real next_freq = model.radiation.frequencies.sorted_nu(next_point, next_freq_idx)*next_shift;
                // std::cout<<"curr_freq: "<<curr_freq<<" next_freq: "<<next_freq<<std::endl;
                comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ, is_in_bounds, in_bounds_count, next_freq_idx, curr_freq_idx, next_freq, curr_freq, next_shift, curr_shift, curr_line_idx, is_upward_disc, bdy_point);

                // freq_matching_()[next_freq_idx]=curr_freq_idx;
                // start_indices_()(nextpointonrayindex,next_freq_idx)[0]=currpointonrayindex;//stores point
                // start_indices_()(nextpointonrayindex,next_freq_idx)[1]=curr_freq_idx;//stores freq
                //for the sake of boundary conditions, these start_indices_ might be overwritten later
                next_freq_idx++;
                in_bounds_count++;
                // curr_freq_count++;
            }
        }
        // is_in_bounds=false;
        //match remaining points as best as possible (with the outermost point)
        // std::cout<<"matching rest of these"<<std::endl;
        while (next_freq_idx<=max_freq_idx)
        {
            // std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
            //TODO REPLACE CURR_FREQ_IDX WITH MAX_FREQ_IDX
            // const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx)*curr_shift;
            const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, max_freq_idx)*curr_shift;
            const Real next_freq = model.radiation.frequencies.sorted_nu(next_point, next_freq_idx)*next_shift;
            is_in_bounds = (next_freq<=curr_freq);
            // std::cout<<"is_in_bounds: "<<is_in_bounds<<std::endl;
            comoving_approx_map_single_data<approx>(model, curr_point, next_point, dZ, is_in_bounds, in_bounds_count, next_freq_idx, curr_freq_idx, next_freq, curr_freq, next_shift, curr_shift, curr_line_idx, is_upward_disc, bdy_point);

            // start_indices_()(nextpointonrayindex,next_freq_idx)[0]=currpointonrayindex;//stores point
            // start_indices_()(nextpointonrayindex,next_freq_idx)[1]=max_freq_idx;//stores freq
            next_freq_idx++;
            in_bounds_count++;
        }
    }
    // std::cout<<"map explicit freq der"<<std::endl;
    //also map the explicit frequency derivative now
    if (is_upward_disc)//determines the direction of the derivative too choose
    {
        //OR iterate over next_freq_idx and compute corresponding curr_freq -> then decide on
        //iterating over curr_freq_idx such that we later on can match using cma_start_frequency_index
        for (Size next_freq_idx = 0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            const Size curr_freq_idx = cma_start_frequency_index_()[next_freq_idx];
            if (curr_freq_idx>=model.parameters->nfreqs()-2)//making sure that we can compute the explicit second order forward frequency derivative (as we need two extra points in that direction)
            {
                cma_dIdnu_expl_()[next_freq_idx] = 0.0;//not enough data, so assume the derivative to be zero (in a slightly better implementation, we could use a first order derivative for the second to last frequency index)
                continue;//next iteration might still contain some interesting stuff to compute
            }
            const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx)*curr_shift;
            const Real curr_freqp1 = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx+1)*curr_shift;
            const Real curr_freqp2 = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx+2)*curr_shift;
            const Real curr_I = cma_computed_intensities_()[curr_freq_idx];
            const Real curr_Ip1 = cma_computed_intensities_()[curr_freq_idx+1];
            const Real curr_Ip2 = cma_computed_intensities_()[curr_freq_idx+2];

            // Real dfreqsmall = cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+1];
            const Real dfreqsmall = curr_freqp2 - curr_freqp1; //cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+1];
            // (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idxp1))
            //           *shift_curr;
            const Real dfreqlarge = curr_freqp2 - curr_freq; //cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+0];
            // (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idx))
            //           *shift_curr;
            const Real dIdnu_coef3_curr = -dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
            const Real dIdnu_coef2_curr = dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
            const Real dIdnu_coef1_curr = -dIdnu_coef3_curr-dIdnu_coef2_curr;//curr point itself

            const Real dIdnu_expl = curr_I*dIdnu_coef1_curr + curr_Ip1 * dIdnu_coef2_curr + curr_Ip2 * dIdnu_coef3_curr;
            cma_dIdnu_expl_()[next_freq_idx] = dIdnu_expl;
        }
    }
    else
    {
        //OR iterate over next_freq_idx and compute corresponding curr_freq -> then decide on
        //iterating over curr_freq_idx such that we later on can match using cma_start_frequency_index
        for (Size next_freq_idx = 0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            const Size curr_freq_idx = cma_start_frequency_index_()[next_freq_idx];
            if (curr_freq_idx<=1)//making sure that we can compute the explicit second order forward frequency derivative (as we need two extra points in that direction)
            {
                cma_dIdnu_expl_()[next_freq_idx] = 0.0;//not enough data, so assume the derivative to be zero (in a slightly better implementation, we could use a first order derivative for the second frequency index)
                continue;//next iteration might still contain some interesting stuff to compute
            }
            const Real curr_freq = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx)*curr_shift;
            const Real curr_freqm1 = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx-1)*curr_shift;
            const Real curr_freqm2 = model.radiation.frequencies.sorted_nu(curr_point, curr_freq_idx-2)*curr_shift;
            const Real curr_I = cma_computed_intensities_()[curr_freq_idx];
            const Real curr_Im1 = cma_computed_intensities_()[curr_freq_idx-1];
            const Real curr_Im2 = cma_computed_intensities_()[curr_freq_idx-2];

            // Real dfreqsmall = cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+1];
            const Real dfreqsmall = curr_freqm2 - curr_freqm1; //cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+1];
            // (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idxp1))
            //           *shift_curr;
            const Real dfreqlarge = curr_freqm2 - curr_freq; //cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+0];
            // (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idx))
            //           *shift_curr;
            const Real dIdnu_coef3_curr = -dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
            const Real dIdnu_coef2_curr = dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
            const Real dIdnu_coef1_curr = -dIdnu_coef3_curr-dIdnu_coef2_curr;//curr point itself

            const Real dIdnu_expl = curr_I*dIdnu_coef1_curr + curr_Im1 * dIdnu_coef2_curr + curr_Im2 * dIdnu_coef3_curr;
            cma_dIdnu_expl_()[next_freq_idx] = dIdnu_expl;
        }
    }
    // std::cout<<"done mapping exp freq der"<<std::endl;
//? should also return the number/indices of boundary frequencies
}

//next and curr freq includes the doppler shift!
// TODO: ALSO ADD CURR AND NEXT POINT
//Err, curr_line_idx is just some value that might or might not be adequate, as we may be connecting multiple lines together. Ignoring for now, as it only impact the single-line-approx (which fails in exactly the same situation of overlapping lines/very high doppler shifts)
template<ApproximationType approx>
accel inline void Solver :: comoving_approx_map_single_data(Model& model, const Size curr_point, const Size next_point, const Real dZ, const bool is_in_bounds, const Size curr_freq_count, const Size next_freq_idx, const Size curr_freq_idx, const Real next_freq, const Real curr_freq, const Real next_shift, const Real curr_shift, const Size curr_line_idx, const bool is_upward_disc, const Size bdy_point)
{
    //next_freq_idx will be used for indexing
    //If bdy, then map onto freq of next index
    if (!is_in_bounds||curr_freq_count<2)
    {
        // std::cout<<"map boundary idx: "<<next_freq_idx<<std::endl;

        // TODO RENAME EVERYTHING
        //overwrite starting intensities? Might be a bad idea, as this contains the computed intensities
        cma_start_intensities_()[next_freq_idx] = boundary_intensity(model, bdy_point, next_freq/curr_shift);
        // std::cout<<"bdy intensity: "<<cma_start_intensities_()[next_freq_idx]<<std::endl;
        // std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
        cma_start_frequencies_()[next_freq_idx] = next_freq;
        cma_start_frequency_index_()[next_freq_idx] = is_upward_disc ? (model.parameters->nfreqs()-1) : 0;//dummy behavior, signaling that this is a bdy point for the different directions, a different frequency is chosen for simplicity (otherwise I nedd another 'is_bdy' field)
        cma_end_frequencies_()[next_freq_idx] = next_freq;

        // compute_source_dtau
        //l, line index (not necessary), or just map freq to line (or use curr_line_idx)
        //curr_shift, next_shift for computing whether the doppler shift is too large
        // so just compute using the line frequency? No, we might change from one line to the next... TODO CHECK OTHER IMPLEMENTATION
        //🌴 other implementation did not in the slightest care about this, so ?
        //technically, it represents a maximal doppler shift, as any frequency matching will try to reduce the shift. In that sense, we might just keep the (slightly inefficient lazy implementation)
        //Otherwise, we might just replace the delta with 2.0*(next_freq-curr_freq)/(next_freq+curr_freq) (technically, the line frequency should be used for the 2/Denominator part)
        // More correct would be: shift_c=1.0, shift_n = 1.0 + (next_freq-curr_freq)/line_freq
        const Real shift_c=curr_shift;
        const Real shift_n=next_shift;
        // const Real shift_c = 1.0;
        // const Real shift_n = 1.0;
        //dz position increment, we should be able to supply this
        //compute_curr_opacity; err, this might become cumbersome to store: simple algorithm: first init to true, then save every result in the next freq list; finally get results from freqmatched curr freq
        //THE LISTS SHOULD BE SEPERATE to avoid overwriting; (note that we probably may ignore the large velocity gradient branch for now, but we should program this properly)
        //and finally stuff to store (use next_freq_idx as index)
        Real Snext = 0.0;
        Real Scurr = 0.0;
        Real dtau = model.parameters->comoving_min_dtau;
        Real chi_next = model.parameters->min_opacity;
        // Real chi_curr = SAVED_CHI[]
        Real chi_curr = 0.0;//As we are out of bounds either way
        bool compute_curr_opacity = true;//we have no information at the bdy freq index, so we formally need to compute the opacity; mainly important if we ever get around to implementing continuum stuff
        // TODO GET CHI_curr (and save CHI_next) not necessary for this

        //as usual, I define the shift backwards
        compute_source_dtau<approx>(model, curr_point, next_point, curr_line_idx, next_freq/shift_c, next_freq/shift_n, shift_c, shift_n, dZ, compute_curr_opacity, dtau, chi_curr, chi_next, Scurr, Snext);
        //cheat on dIdnu coeffs, as these are all zero on the boundary
        // cma_dIdnu_coef1_curr_()[next_freq_idx] = 0;
        // cma_dIdnu_coef2_curr_()[next_freq_idx] = 0;
        // cma_dIdnu_coef3_curr_()[next_freq_idx] = 0;
        cma_dIdnu_coef1_next_()[next_freq_idx] = 0;
        cma_dIdnu_coef2_next_()[next_freq_idx] = 0;
        cma_dIdnu_coef3_next_()[next_freq_idx] = 0;

        cma_chi_next_()[next_freq_idx] = chi_next;
        cma_compute_next_opacity_()[next_freq_idx] = compute_curr_opacity;//for the next thing
        // cma_chi_curr_()[curr_freq_idx] = chi_curr;??is not part of default frequency discretization
        // cma_compute_curr_opacity_()[curr_freq_idx] = compute_curr_opacity;//if we somehow encounter a current frequency index twice, we set this


        cma_S_next_()[next_freq_idx]=Snext;
        cma_S_curr_()[next_freq_idx]=Scurr;
        // const Real dtau = std::max(trap(chi_curr, chi_next, dZ), COMOVING_MIN_DTAU);
        // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
        dtau = std::max(dtau, model.parameters->comoving_min_dtau);
        cma_delta_tau_()[next_freq_idx]=dtau;
        // chi_temp_copy_()(next_freq_idx)=chi_next;
    }
    else
    {
        // std::cout<<"map nonboundary"<<std::endl;
        //in this case, we need to compute a bit more
        //overwrite starting intensities? Might be a bad idea, as this contains the computed intensities
        cma_start_intensities_()[next_freq_idx] = cma_computed_intensities_()[curr_freq_idx];
        cma_start_frequencies_()[next_freq_idx] = curr_freq;
        cma_start_frequency_index_()[next_freq_idx] = curr_freq_idx;
        cma_end_frequencies_()[next_freq_idx] = next_freq;
        // compute_source_dtau
        //l, line index (not necessary), or just map freq to line (or use curr_line_idx)
        //curr_shift, next_shift for computing whether the doppler shift is too large
        // so just compute using the line frequency? No, we might change from one line to the next... TODO CHECK OTHER IMPLEMENTATION
        //🌴 other implementation did not in the slightest care about this, so ?
        //technically, it represents a maximal doppler shift, as any frequency matching will try to reduce the shift. In that sense, we might just keep the (slightly inefficient lazy implementation)
        //Otherwise, we might just replace the delta with 2.0*(next_freq-curr_freq)/(next_freq+curr_freq) (technically, the line frequency should be used for the 2/Denominator part)
        // More correct would be: shift_c=1.0, shift_n = 1.0 + (next_freq-curr_freq)/line_freq
        //As we may use any frame of reference, we might as well put it a the current point in order to decide whether the doppler shift is large
        const Real shift_c = curr_shift;
        const Real shift_n = next_shift;
        // const Real shift_c=1.0;
        // const Real line_freq = model.lines.line[curr_line_idx];
        // const Real shift_n=1.0 + (next_freq-curr_freq)/line_freq;

        // std::cout<<"computed shift"<<std::endl;
        //dz position increment, we should be able to supply this
        //compute_curr_opacity; err, this might become cumbersome to store: simple algorithm: first init to true, then save every result in the next freq list; finally get results from freqmatched curr freq
        //THE LISTS SHOULD BE SEPERATE to avoid overwriting; (note that we probably may ignore the large velocity gradient branch for now, but we should program this properly)
        //and finally stuff to store (use next_freq_idx as index)
        Real Snext = 0.0;
        Real dtau = model.parameters->comoving_min_dtau;
        Real chi_next = model.parameters->min_opacity;
        Real Scurr = cma_S_curr_()[curr_freq_idx];//if we have already computed the opacity in the low doppler shift case, the current source will not be computed...
        // Real chi_curr = SAVED_CHI[]
        Real chi_curr = cma_chi_curr_()[curr_freq_idx];//retrieving the correct frequency
        bool compute_curr_opacity = cma_compute_curr_opacity_()[curr_freq_idx];
        // std::cout<<"shift_c: "<<shift_c<<" shift_n: "<<shift_n<<std::endl;

        // std::cout<<"retrieved curr opacity"<<std::endl;

        compute_source_dtau<approx>(model, curr_point, next_point, curr_line_idx, curr_freq/shift_c, next_freq/shift_n, shift_c, shift_n, dZ, compute_curr_opacity, dtau, chi_curr, chi_next, Scurr, Snext);

        // std::cout<<"computed sdtau"<<std::endl;

        cma_S_next_()[next_freq_idx]=Snext;
        cma_S_curr_()[next_freq_idx]=Scurr;
        // std::cout<<"snext: "<<cma_S_next_()[next_freq_idx]<<" scurr: "<<cma_S_curr_()[next_freq_idx]<<std::endl;

        cma_chi_next_()[next_freq_idx] = chi_next;
        cma_chi_curr_()[curr_freq_idx] = chi_curr;
        cma_compute_curr_opacity_()[curr_freq_idx] = compute_curr_opacity;//if we somehow encounter a current frequency index twice, we set this
        cma_compute_next_opacity_()[next_freq_idx] = compute_curr_opacity;//for the next thing

        // std::cout<<"set opacities"<<std::endl;

        // const Real dtau = std::max(trap(chi_curr, chi_next, dZ), COMOVING_MIN_DTAU);
        // Floor dtau by COMOVING_MIN_DTAU, due to division by dtau^2
        dtau = std::max(dtau, model.parameters->comoving_min_dtau);
        cma_delta_tau_()[next_freq_idx]=dtau;
        // std::cout<<"deltatau: "<<dtau<<std::endl;
        // chi_temp_copy_()(next_freq_idx)=chi_next;
        // MAP CHINEXT TO USEFUL ARRAY
        // AT END OF FULL MAP, COPY THIS ARRAY *DONE* at end of single position increment

        //⚠ Assumes the frequencies to be set in a fixed order!
        // if is_upward_disc==TRUE, then the frequency indices should be traversed from nfreqs-1 to 0
        // if is_upward_disc==FALSE, then the frequency indices should be traversed from 0 to nfreqs-1
        if (is_upward_disc)
        {
            // // //in order to circumvent annoying frequency derivative stuff bdy stuff, we explicitly compute the freq derivative now
            // // const real curr_freq_intensity = cma_computed_intensities_()[curr_freq];
            // // const Size curr_freq_indexp1 = curr_freq_idx + 1;
            // // const Real curr_freqp1 = TODOcurr_freq_idx + 1;
            // // const Real curr_intensityp1 = cma_computed_intensities_()[curr_freq_idx + 1];
            // // const Size curr_freqp2 = (curr_freq_idx + 1 == model.parameters->nfreqs() - 1) ? cma_start_frequencies_()[last_used_bdy_index]: curr_freq_idx+2;
            // // const Real curr_intensityp2 = cma_start_intensities_()[last_used_bdy_index]: cma_computed_intensities_()[curr_freq_idx + 2];
            //
            // // Real dfreqsmall = cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+1];
            // Real dfreqsmall = cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+1];
            // // (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idxp1))
            // //           *shift_curr;
            // Real dfreqlarge = cma_start_frequencies_()[next_freq_idx+2] - cma_start_frequencies_()[next_freq_idx+0];
            // // (model.radiation.frequencies.nu(currpoint, curr_freq_idxp2)-model.radiation.frequencies.nu(currpoint, curr_freq_idx))
            // //           *shift_curr;
            // cma_dIdnu_coef3_curr_()[next_freq_idx] = -dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
            // cma_dIdnu_coef2_curr_()[next_freq_idx] = dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
            // cma_dIdnu_coef1_curr_()[next_freq_idx] = -cma_dIdnu_coef3_curr_()[next_freq_idx]-cma_dIdnu_coef2_curr_()[next_freq_idx];//curr point itself
            // //index stuff is actually not necessary, as for all points except the boundary, we just need to increment the index.
            // // for the boundary, we do not at all use these values, so this should be fine
            // // dIdnu_index3_curr_()(next_freq_idx)=curr_freq_idxp2;
            // // dIdnu_index2_curr_()(next_freq_idx)=curr_freq_idxp1;
            // // dIdnu_index1_curr_()(next_freq_idx)=curr_freq_idx;

            //And now do the same (but simpler; less index management) for the implicit part
            Real dfreqsmall = cma_end_frequencies_()[next_freq_idx+2] - cma_end_frequencies_()[next_freq_idx+1];
            // (model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx+1))
            //             *shift_next;
            Real dfreqlarge = cma_end_frequencies_()[next_freq_idx+2] - cma_end_frequencies_()[next_freq_idx+0];
            // (model.radiation.frequencies.nu(nextpoint, next_freq_idx+2)-model.radiation.frequencies.nu(nextpoint, next_freq_idx))
            //             *shift_next;

            cma_dIdnu_coef3_next_()[next_freq_idx] = -dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
            cma_dIdnu_coef2_next_()[next_freq_idx] = dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
            cma_dIdnu_coef1_next_()[next_freq_idx] = -cma_dIdnu_coef3_next_()[next_freq_idx]-cma_dIdnu_coef2_next_()[next_freq_idx];//curr point itself;

        }
        else
        {
            // std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
            // //Note: minus signs in front, as we are now computing the first derivative using points on the other side (switched order of term, making sign clearer; just compare to upward disc version)
            // Real dfreqsmall = -(cma_start_frequencies_()[next_freq_idx-1] - cma_start_frequencies_()[next_freq_idx-2]);
            // std::cout<<"dfreqsmall: "<<dfreqsmall<<std::endl;
            // // -(model.radiation.frequencies.nu(currpoint, curr_freq_idxm1)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
            // // *shift_curr;
            // Real dfreqlarge = -(cma_start_frequencies_()[next_freq_idx-0] - cma_start_frequencies_()[next_freq_idx-2]);
            // std::cout<<"dfreqlarge: "<<dfreqlarge<<std::endl;
            // // -(model.radiation.frequencies.nu(currpoint, curr_freq_idx)-model.radiation.frequencies.nu(currpoint, curr_freq_idxm2))
            // // *shift_curr;
            // cma_dIdnu_coef3_curr_()[next_freq_idx] = -dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
            // cma_dIdnu_coef2_curr_()[next_freq_idx] = dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
            // cma_dIdnu_coef1_curr_()[next_freq_idx] = -cma_dIdnu_coef3_curr_()[next_freq_idx]-cma_dIdnu_coef2_curr_()[next_freq_idx];//curr point itself
            // // dIdnu_index3_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm2;
            // // dIdnu_index2_curr_()(rayposidx, next_freq_idx)=curr_freq_idxm1;
            // // dIdnu_index1_curr_()(rayposidx, next_freq_idx)=curr_freq_idx;

            //And now do the same (but simpler; less index management) for the implicit part
            //Note: minus signs in front, as we are now computing the first derivative using points on the other side
            Real dfreqsmall = -(cma_end_frequencies_()[next_freq_idx-1] - cma_end_frequencies_()[next_freq_idx-2]);
            // -(model.radiation.frequencies.nu(nextpoint, next_freq_idx-1)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
            // *shift_next;
            Real dfreqlarge = -(cma_end_frequencies_()[next_freq_idx-0] - cma_end_frequencies_()[next_freq_idx-2]);
            // -(model.radiation.frequencies.nu(nextpoint, next_freq_idx)-model.radiation.frequencies.nu(nextpoint, next_freq_idx-2))
            // *shift_next;
            cma_dIdnu_coef3_next_()[next_freq_idx] = -dfreqsmall/(std::pow(dfreqlarge, 2.0)-dfreqlarge*dfreqsmall);//farthest
            cma_dIdnu_coef2_next_()[next_freq_idx] = dfreqlarge/(-std::pow(dfreqsmall, 2.0)+dfreqlarge*dfreqsmall);//nearer
            cma_dIdnu_coef1_next_()[next_freq_idx] = -cma_dIdnu_coef3_next_()[next_freq_idx]-cma_dIdnu_coef2_next_()[next_freq_idx];//curr point itself;
            // dIdnu_index3_next_()(rayposidx, next_freq_idx)=next_freq_idx-2;
            // dIdnu_index2_next_()(rayposidx, next_freq_idx)=next_freq_idx-1;
            // dIdnu_index1_next_()(rayposidx, next_freq_idx)=next_freq_idx;

        }
        //end of non-boundary data mapping
    }
}


//As stepping a single step should be exactly the same in both directions on the ray, we might as well refactor it out
// //rayposidx==rayposidx_curr_point+-1
//Ergo rayposidx stand for the ray position index of the next point
//rr direction index necessary for determining whether to add intensity to J//TODO? replace with bool denoting whether to add it (precompute if clause somewhere else)
//TODO: figure out whether rayposidx_currpoint is needed!! As it should be replaced with start_indices_...
// inline void Solver :: solve_comoving_single_step (Model& model, const Size raypos-idx, const Size rayidx, const Size rr, const bool is_upward_disc, const bool forward_ray)
accel inline void Solver :: solve_comoving_local_approx_single_step (Model& model, const Size next_point, const Size rayidx, const Size rr, const bool is_upward_disc)
{
    // std::cout<<"is_upward_d isc: "<<is_upward_disc<<std::endl;
    // IF BORED: refactor this function such that it acts the explicit and implicit part on the same time (not first forward for all freqs, than backward for all freqs)
    Vector<Size>& nr=nr_();//stores the exact point indices
    // Vector<unsigned char>& real_pt=real_pt_();//stores whether each point is real point (not added extra for interpolation purposes)
    Vector<double>& shift=shift_();//stores the shift; should already be in the correct direction
    Vector<Real>& cma_delta_tau=cma_delta_tau_();//stores the optical depths
    Vector<Real>& cma_S_curr=cma_S_curr_();
    Vector<Real>& cma_S_next=cma_S_next_();
    Vector<Real>& cma_computed_intensities=cma_computed_intensities_();
    Vector<Real>& cma_start_intensities=cma_start_intensities_();
    Vector<Real>& cma_start_frequencies=cma_start_frequencies_();
    Vector<Real>& cma_end_frequencies=cma_end_frequencies_();
    // TODO: at end of this, do set the start equal to the computed intensities; THIS IS REQUIRED FOR THE NEXT MAPPING STEP
    //*DONE*

    Vector<Real>& cma_dIdnu_expl = cma_dIdnu_expl_();
    // Vector<Real>& cma_dIdnu_coef1_curr=cma_dIdnu_coef1_curr_();//coefficient of the freq point itself
    // Vector<Real>& cma_dIdnu_coef2_curr=cma_dIdnu_coef2_curr_();
    // Vector<Real>& cma_dIdnu_coef3_curr=cma_dIdnu_coef3_curr_();

    // Vector<Size>& cma_start_frequency_index=cma_start_frequency_index_();//direction will tell me anyway which values to use (but do use modulo? to prevent oob access)

    Vector<Real>& cma_dIdnu_coef1_next=cma_dIdnu_coef1_next_();
    Vector<Real>& cma_dIdnu_coef2_next=cma_dIdnu_coef2_next_();
    Vector<Real>& cma_dIdnu_coef3_next=cma_dIdnu_coef3_next_();

    // Matrix<Size>& dIdnu_index1_next=dIdnu_index1_next_();
    // Matrix<Size>& dIdnu_index2_next=dIdnu_index2_next_();
    // Matrix<Size>& dIdnu_index3_next=dIdnu_index3_next_();

    // std::cout<<"after renaming all vars"<<std::endl;
    // std::cout<<"rayposidx=next point on ray idx= "<<rayposidx<<std::endl;

    //rayposidx represents the next position on the ray
    const Size nextpointidx=next_point;//nr[rayposidx];
    // const Real next_shift=(forward_ray) ? 2.0-shift[rayposidx] : shift[rayposidx];//absolute shifts, versus static frame!

    // std::cout<<"after getting shift"<<std::endl;

    // for (Size next_freq_idx = 0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
    // {
    //   std::cout<<"start intensities: "<<cma_start_intensities[next_freq_idx]<<std::endl;
    // }

    //Do the explicit part
    //discretization determines which points are most definitely boundary points (for which the computation is simpler)
    if (is_upward_disc)
    {
        // DO TWO BOUNDARY CONDITIONS at model.parameters->nfreqs()-1,-2 (err, 0, 1 instead here)
        //outermost boundary
        const Real bdy_dtau1 = cma_delta_tau[model.parameters->nfreqs()-1];
        const Real bdy_expl_term1 = (-expm1(-bdy_dtau1)-bdy_dtau1*exp(-bdy_dtau1))/bdy_dtau1;
        cma_computed_intensities[model.parameters->nfreqs()-1]=cma_start_intensities[model.parameters->nfreqs()-1]*exp(-bdy_dtau1)
        +bdy_expl_term1*cma_S_curr[model.parameters->nfreqs()-1]; //source term //no frequency derivative term, as the frequency shift of the boundary is by definition 0
        //second boundary
        const Real bdy_dtau2 = cma_delta_tau[model.parameters->nfreqs()-2];
        const Real bdy_expl_term2 = (-expm1(-bdy_dtau2)-bdy_dtau2*exp(-bdy_dtau2))/bdy_dtau2;
        cma_computed_intensities[model.parameters->nfreqs()-2]=cma_start_intensities[model.parameters->nfreqs()-2]*exp(-bdy_dtau2)
        +bdy_expl_term2*cma_S_curr[model.parameters->nfreqs()-2]; //source term //no frequency derivative term, as the frequency shift of the boundary is by definition 0
        // std::cout<<"bdy start I: "<<cma_start_intensities[model.parameters->nfreqs()-1]<<cma_start_intensities[model.parameters->nfreqs()-2]<<std::endl;

        // for (Size next_freq_idx=2; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs()-2; next_freq_idx++)
        {
            // std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
            // std::cout<<"start I: "<<cma_start_intensities[next_freq_idx]<<std::endl;
            // const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx)[0];
            // std::cout<<"curr_point_on_ray_index: "<<curr_point_on_ray_index<<std::endl;
            // const Size curr_point_idx=nr[curr_point_on_ray_index];
            // const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            // const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame

            // std::cout<<"after computing shift"<<std::endl;

            // const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real deltanu=cma_end_frequencies[next_freq_idx]-cma_start_frequencies[next_freq_idx];
            // std::cout<<"deltanu: "<<deltanu<<std::endl;

            // const Real dtau=delta_tau(rayposidx, next_freq_idx);
            const Real dtau=cma_delta_tau[next_freq_idx];
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"rayposidx: "<<rayposidx<<std::endl;

            // std::cout<<"after computing dtau"<<std::endl;

            const Real expl_term=(-expm1(-dtau)-dtau*exp(-dtau))/dtau;
            // std::cout<<"expl_term: "<<expl_term<<std::endl;
            // std::cout<<"expl term: "<<expl_term<<std::endl;
            // std::cout<<"curr source: "<<S_curr(rayposidx, next_freq_idx)<<std::endl;
            // const Size curr_freq_idx = cma_start_frequency_index_()[next_freq_idx];

            const Real expl_freq_der=cma_dIdnu_expl[next_freq_idx];
            // const Real expl_freq_der=cma_dIdnu_coef1_curr[next_freq_idx]*cma_start_intensities[next_freq_idx]
            //                         +cma_dIdnu_coef2_curr[next_freq_idx]*cma_start_intensities[next_freq_idx+1]
            //                         +cma_dIdnu_coef3_curr[next_freq_idx]*cma_start_intensities[next_freq_idx+2];
            // std::cout<<"expl_freq_der: "<<expl_freq_der<<std::endl;
            // std::cout<<"curr freq der coef 1: "<<dIdnu_coef1_curr(rayposidx, next_freq_idx)<<std::endl;
            // std::cout<<"curr freq der coef 2: "<<dIdnu_coef2_curr(rayposidx, next_freq_idx)<<std::endl;
            // std::cout<<"curr freq der coef 3: "<<dIdnu_coef3_curr(rayposidx, next_freq_idx)<<std::endl;
            // std::cout<<"curr_point_on_ray_index: "<<curr_point_on_ray_index<<std::endl;
            // std::cout<<"respective freq index: "<<dIdnu_index1_curr(rayposidx, next_freq_idx)<<std::endl;
            // std::cout<<"intensity coef 1: "<<intensities(curr_point_on_ray_index, dIdnu_index1_curr(rayposidx, next_freq_idx))<<std::endl;
            // std::cout<<"intensity coef 2: "<<intensities(curr_point_on_ray_index, dIdnu_index2_curr(rayposidx, next_freq_idx))<<std::endl;
            // std::cout<<"intensity coef 3: "<<intensities(curr_point_on_ray_index, dIdnu_index3_curr(rayposidx, next_freq_idx))<<std::endl;
            // std::cout<<"expl_freq_der: "<<expl_freq_der<<std::endl;
            // std::cout<<"full freq der term: "<<expl_term*expl_freq_der*deltanu/dtau<<std::endl;
            //TODO: actually do the explicit part
            // std::cout<<"intensities before explicit part: "<<intensities(curr_point_on_ray_index, curr_freq_idx)<<std::endl;
            // std::cout<<"intensities after optical depth: "<<intensities(curr_point_on_ray_index, curr_freq_idx)*exp(-dtau)<<std::endl;
            cma_computed_intensities[next_freq_idx]=cma_start_intensities[next_freq_idx]*exp(-dtau)
            +expl_term*cma_S_curr[next_freq_idx]//source term
            +expl_term*expl_freq_der*deltanu/dtau;//frequency derivative term
            // std::cout<<"explicit part intensities: "<<cma_computed_intensities[next_freq_idx]<<std::endl;
            //Get the point indices
        }
    }
    else
    {//!is_upward_disc
        // DO TWO BOUNDARY CONDITIONS at 0, 1
        //outermost boundary
        const Real bdy_dtau1 = cma_delta_tau[0];
        const Real bdy_expl_term1 = (-expm1(-bdy_dtau1)-bdy_dtau1*exp(-bdy_dtau1))/bdy_dtau1;
        cma_computed_intensities[0]=cma_start_intensities[0]*exp(-bdy_dtau1)
        +bdy_expl_term1*cma_S_curr[0]; //source term //no frequency derivative term, as the frequency shift of the boundary is by definition 0
        //second boundary
        const Real bdy_dtau2 = cma_delta_tau[1];
        const Real bdy_expl_term2 = (-expm1(-bdy_dtau2)-bdy_dtau2*exp(-bdy_dtau2))/bdy_dtau2;
        cma_computed_intensities[1]=cma_start_intensities[1]*exp(-bdy_dtau2)
        +bdy_expl_term2*cma_S_curr[1]; //source term //no frequency derivative term, as the frequency shift of the boundary is by definition 0

        // std::cout<<"computed explicit boundary"<<std::endl;
        // std::cout<<"bdy start I: "<<cma_start_intensities[0]<<cma_start_intensities[1]<<std::endl;
        // for (Size next_freq_idx=2; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        for (Size next_freq_idx=2; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            // std::cout<<"start of loop"<<std::endl;
            // std::cout<<"next_freq_idx: "<<next_freq_idx<<std::endl;
            // std::cout<<"start I: "<<cma_start_intensities[next_freq_idx]<<std::endl;
            // const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx)[0];
            // std::cout<<"curr_point_on_ray_index: "<<curr_point_on_ray_index<<std::endl;
            // const Size curr_point_idx=nr[curr_point_on_ray_index];
            // const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            // const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame

            // std::cout<<"after computing shift"<<std::endl;

            // const Real deltanu=model.radiation.frequencies.nu(nextpointidx, next_freq_idx)*next_shift-model.radiation.frequencies.nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real deltanu=cma_end_frequencies[next_freq_idx]-cma_start_frequencies[next_freq_idx];

            // const Real dtau=delta_tau(rayposidx, next_freq_idx);
            const Real dtau=cma_delta_tau[next_freq_idx];
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"rayposidx: "<<rayposidx<<std::endl;

            // std::cout<<"after computing dtau"<<std::endl;

            const Real expl_term=(-expm1(-dtau)-dtau*exp(-dtau))/dtau;
            // std::cout<<"expl term: "<<expl_term<<std::endl;
            // std::cout<<"curr source: "<<S_curr(rayposidx, next_freq_idx)<<std::endl;
            // const Size curr_freq_idx = cma_start_frequency_index_()[next_freq_idx];
            // std::cout<<"curr_freq_idx: "<<curr_freq_idx<<std::endl;

            // const Real expl_freq_der=cma_dIdnu_coef1_curr[next_freq_idx]*cma_start_intensities[next_freq_idx]
            //                         +cma_dIdnu_coef2_curr[next_freq_idx]*cma_start_intensities[next_freq_idx-1]
            //                         +cma_dIdnu_coef3_curr[next_freq_idx]*cma_start_intensities[next_freq_idx-2];
            const Real expl_freq_der=cma_dIdnu_expl[next_freq_idx];
            // std::cout<<"curr freq der coef 1: "<<cma_dIdnu_coef1_curr[next_freq_idx]<<std::endl;
            // std::cout<<"curr freq der coef 2: "<<cma_dIdnu_coef2_curr[next_freq_idx]<<std::endl;
            // std::cout<<"curr freq der coef 3: "<<cma_dIdnu_coef3_curr[next_freq_idx]<<std::endl;
            // std::cout<<"curr_point_on_ray_index: "<<curr_point_on_ray_index<<std::endl;
            // std::cout<<"respective freq index: "<<dIdnu_index1_curr(rayposidx, next_freq_idx)<<std::endl;
            // std::cout<<"intensity coef 1: "<<intensities(curr_point_on_ray_index, dIdnu_index1_curr(rayposidx, next_freq_idx))<<std::endl;
            // std::cout<<"intensity coef 2: "<<intensities(curr_point_on_ray_index, dIdnu_index2_curr(rayposidx, next_freq_idx))<<std::endl;
            // std::cout<<"intensity coef 3: "<<intensities(curr_point_on_ray_index, dIdnu_index3_curr(rayposidx, next_freq_idx))<<std::endl;
            // std::cout<<"expl_freq_der: "<<expl_freq_der<<std::endl;
            // std::cout<<"full freq der term: "<<expl_term*expl_freq_der*deltanu/dtau<<std::endl;
            //TODO: actually do the explicit part
            // std::cout<<"intensities before explicit part: "<<intensities(curr_point_on_ray_index, curr_freq_idx)<<std::endl;
            // std::cout<<"intensities after optical depth: "<<intensities(curr_point_on_ray_index, curr_freq_idx)*exp(-dtau)<<std::endl;
            cma_computed_intensities[next_freq_idx]=cma_start_intensities[next_freq_idx]*exp(-dtau)
            +expl_term*cma_S_curr[next_freq_idx]//source term
            +expl_term*expl_freq_der*deltanu/dtau;//frequency derivative term
            // std::cout<<"explicit part intensities: "<<cma_computed_intensities[next_freq_idx]<<std::endl;
            //Get the point indices
        }
    }

    // std::cout<<"after explicit part"<<std::endl;
    //Implicit part ordering depends on the discretization direction
    if (is_upward_disc)
    {
        //TODO DO BOUNDARY STUFF! AT nfreqs-1,-2
        const Real bdy_dtau1=cma_delta_tau[model.parameters->nfreqs()-1];
        const Real impl_term1=(bdy_dtau1+expm1(-bdy_dtau1))/bdy_dtau1;
        cma_computed_intensities[model.parameters->nfreqs()-1]=cma_computed_intensities[model.parameters->nfreqs()-1]
        +impl_term1*cma_S_next[model.parameters->nfreqs()-1];//source term // no freq derivative term
        // +impl_term*impl_freq_der*deltanu/dtau)/(1.0-dIdnu_coef1_next(rayposidx, next_freq_idx-1)*impl_term*deltanu/dtau); //freq derivative term
        const Real bdy_dtau2=cma_delta_tau[model.parameters->nfreqs()-2];
        const Real impl_term2=(bdy_dtau2+expm1(-bdy_dtau2))/bdy_dtau2;
        cma_computed_intensities[model.parameters->nfreqs()-2]=cma_computed_intensities[model.parameters->nfreqs()-2]
        +impl_term2*cma_S_next[model.parameters->nfreqs()-2];//source term

        // std::cout<<"is_upward_disc"<<std::endl;
        //in case of positive doppler shifts, we need to start from the largest frequencies (where we have extra boundary conditions on the end)
        //reverse loop over unsigned int, so indices +1
        // for (Size next_freq_idx=model.parameters->nfreqs(); next_freq_idx>0; next_freq_idx--)
        for (Size next_freq_idx=model.parameters->nfreqs()-2; next_freq_idx>0; next_freq_idx--)
        {
            //TODO: actually do the implicit part
            // const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx-1)[0];
            // const Size curr_point_idx=nr[curr_point_on_ray_index];
            // const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx-1)[1];
            // const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            const Real deltanu=cma_end_frequencies[next_freq_idx-1]-cma_start_frequencies[next_freq_idx-1];

            // const Real dtau=delta_tau(rayposidx, next_freq_idx);
            const Real dtau=cma_delta_tau[next_freq_idx-1];            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"next freq idx: "<<next_freq_idx-1<<std::endl;
            // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
            const Real impl_term=(dtau+expm1(-dtau))/dtau;
            // std::cout<<"impl term: "<<impl_term<<std::endl;
            // std::cout<<"next source: "<<S_next(rayposidx, next_freq_idx-1)<<std::endl;
            //only the parts of the other points (to subtract/add)
            const Real impl_freq_der=cma_dIdnu_coef2_next[next_freq_idx-1]*cma_computed_intensities[next_freq_idx+1-1]
                                    +cma_dIdnu_coef3_next[next_freq_idx-1]*cma_computed_intensities[next_freq_idx+2-1];
            // std::cout<<"impl freq der: "<<impl_freq_der<<std::endl;

            // std::cout<<"intensities before implicit part: "<<intensities(rayposidx, next_freq_idx-1)<<std::endl;

            cma_computed_intensities[next_freq_idx-1]=(cma_computed_intensities[next_freq_idx-1]
            +impl_term*cma_S_next[next_freq_idx-1]//source term
            +impl_term*impl_freq_der*deltanu/dtau)/(1.0-cma_dIdnu_coef1_next[next_freq_idx-1]*impl_term*deltanu/dtau); //freq derivative term

            // std::cout<<"implicit part intensities: "<<cma_computed_intensities[next_freq_idx-1]<<std::endl;
        }
    }
    else
    {
        // std::cout<<"implicit part"<<std::endl;
        //TODO DO BOUNDARY STUFF! AT 0,1
        const Real bdy_dtau1=cma_delta_tau[0];
        const Real impl_term1=(bdy_dtau1+expm1(-bdy_dtau1))/bdy_dtau1;
        cma_computed_intensities[0]=cma_computed_intensities[0]
        +impl_term1*cma_S_next[0];//source term // no freq derivative term

        const Real bdy_dtau2=cma_delta_tau[1];
        const Real impl_term2=(bdy_dtau2+expm1(-bdy_dtau2))/bdy_dtau2;
        cma_computed_intensities[1]=cma_computed_intensities[1]
        +impl_term2*cma_S_next[1];//source term // no freq derivative term

        // std::cout<<"set boundary"<<std::endl;

        for (Size next_freq_idx=2; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
        {
            //TODO: actually do the implicit part
            // const Size curr_point_on_ray_index=start_indices_()(rayposidx, next_freq_idx)[0];
            // const Size curr_point_idx=nr[curr_point_on_ray_index];
            // const Size curr_freq_idx=start_indices_()(rayposidx, next_freq_idx)[1];
            // const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            const Real deltanu=cma_end_frequencies[next_freq_idx]-cma_start_frequencies[next_freq_idx];

            // const Real dtau=delta_tau(rayposidx, next_freq_idx);
            const Real dtau=cma_delta_tau[next_freq_idx];
            // std::cout<<"dtau: "<<dtau<<std::endl;
            // std::cout<<"next freq idx: "<<next_freq_idx<<std::endl;
            const Real impl_term=(dtau+expm1(-dtau))/dtau;
            // std::cout<<"impl_term: "<<impl_term<<std::endl;
            //only the parts of the other points (to subtract/add)
            const Real impl_freq_der=cma_dIdnu_coef2_next[next_freq_idx]*cma_computed_intensities[next_freq_idx-1]
                                    +cma_dIdnu_coef3_next[next_freq_idx]*cma_computed_intensities[next_freq_idx-2];
            // std::cout<<"impl freq der: "<<impl_freq_der<<std::endl;
            // std::cout<<"full freq der term: "<<impl_term*impl_freq_der*deltanu/dtau<<std::endl;
            // std::cout<<"1.0-dIdnu coef times...: "<<1.0-dIdnu_coef1_next(rayposidx, next_freq_idx)*deltanu/dtau<<std::endl;
            //
            // std::cout<<"intensities before implicit part: "<<intensities(rayposidx, next_freq_idx)<<std::endl;

            cma_computed_intensities[next_freq_idx]=(cma_computed_intensities[next_freq_idx]
            +impl_term*cma_S_next[next_freq_idx]//source term
            +impl_term*impl_freq_der*deltanu/dtau)/(1.0-cma_dIdnu_coef1_next[next_freq_idx]*impl_term*deltanu/dtau); //freq derivative term

            // std::cout<<"implicit part intensities: "<<cma_computed_intensities[next_freq_idx]<<std::endl;
        }
    }
    // for (Size next_freq_idx = 0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
    // {
    //   std::cout<<"implicit part intensities: "<<cma_computed_intensities[next_freq_idx]<<std::endl;
    // }


    // std::cout<<"incrementing J"<<std::endl;

    //Finally increment J if and real point check if closest ray
    // if (real_pt[rayposidx]&&closest_ray(rr, nextpointidx)==rayidx)
    if (closest_ray(rr, nextpointidx)==rayidx)
    {
        // std::cout<<"rayposidx: "<<rayposidx<<std::endl;
        // std::cout<<"adding to J at rayposidx: "<<rayposidx<<std::endl;
        const Real     wt = model.geometry.rays.weight   [rr];
        //then obviously add (weighted) to J
        for (Size freqid=0; freqid<model.parameters->nfreqs(); freqid++)
        {
            // Size1 corresponding_l_for_spec;           ///< number of line species corresponding to frequency
            // Size1 corresponding_k_for_tran;           ///< number of transition corresponding to frequency
            // Size1 corresponding_z_for_line;           ///< quadrature number corresponding to frequency
            // TODO: replace with matrix variant when replacing these Vectors
            const Size unsorted_freqidx=model.radiation.frequencies.corresponding_nu_index(nextpointidx, freqid);
            const Size l=model.radiation.frequencies.corresponding_l_for_spec[unsorted_freqidx];
            const Size k=model.radiation.frequencies.corresponding_k_for_tran[unsorted_freqidx];
            const Size z=model.radiation.frequencies.corresponding_z_for_line[unsorted_freqidx];
            LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
            // std::cout<<"freqid: "<<freqid<<std::endl;
            // std::cout<<"adding: "<<lspec.quadrature.weights[z] * wt * intensities(rayposidx, freqid)<<std::endl;
            // std::cout<<"intensity: "<<intensities(rayposidx, freqid)<<std::endl;
            lspec.J(nextpointidx,k) += lspec.quadrature.weights[z] * wt * cma_computed_intensities[freqid];// Su_()[centre];

            //Computing the ALI element

            //TODO: check whether this is correct (print out Jeff, Jlin in optically thick model)
            // const Size curr_point_on_ray_index=start_indices_()(rayposidx, freqid)[0];
            // const Size curr_point_idx=nr[curr_point_on_ray_index];
            // const Size curr_freq_idx=start_indices_()(rayposidx, freqid)[1];
            // const Real curr_shift=(forward_ray) ? 2.0-shift[curr_point_on_ray_index] : shift[curr_point_on_ray_index];//absolute shift, versus static frame
            // const Real deltanu=model.radiation.frequencies.sorted_nu(nextpointidx, freqid)*next_shift-model.radiation.frequencies.sorted_nu(curr_point_idx, curr_freq_idx)*curr_shift;
            const Real deltanu=cma_end_frequencies[freqid]-cma_start_frequencies[freqid];
            const Real dtau=cma_delta_tau[freqid];
            const Real constant = lspec.linedata.A[k] * lspec.quadrature.weights[z] * model.geometry.rays.weight[rr];//for integrating over both frequency (weighted with profile function) and angle. We also need the einstein A coefficient

            const Real source_term=(dtau+expm1(-dtau));// '/dtau' has already been applied to simplify the fraction
            const Real lambdaterm=constant*source_term/(dtau-cma_dIdnu_coef1_next[freqid]*source_term*deltanu);

            // frq = freqs.nu(nr[n], f) * shift[n];
            // phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
            // L   = constante * frq * phi * L_lower(m,n) * inverse_chi[n];

            lspec.lambda.add_element(nextpointidx, k, nextpointidx, lambdaterm);
        }
    }

    //Finally overwrite some variables with their value at the next point
    // in this way, we prepare for the next step
    for (Size next_freq_idx=0; next_freq_idx<model.parameters->nfreqs(); next_freq_idx++)
    {
        //note: for now, I might overwrite a bit too much stuff (which will be recomputed either way)
        //but better be safe than sorry for this complicated solver; I can always check later which parts may safely be not overwritten
        ///these three things are not indexed using curr_freq_idx and are purely some helper vars for storing temp results; so these thing may not be overwritten
        // cma_start_intensities[next_freq_idx]=cma_computed_intensities[next_freq_idx];//is always corrected when mapping the data itself
        // cma_start_frequencies[next_freq_idx]=cma_end_frequencies[next_freq_idx];
        cma_S_curr[next_freq_idx]=cma_S_next[next_freq_idx];
        cma_chi_curr_()[next_freq_idx] = cma_chi_next_()[next_freq_idx];
        cma_compute_curr_opacity_()[next_freq_idx] = cma_compute_next_opacity_()[next_freq_idx];
    }
}




template<ApproximationType approx>
inline void Solver :: solve_shortchar_order_0 (Model& model)
{
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities)
    {
        model.radiation.I.resize (model.parameters->nrays(),  model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize (                           model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    // Initialise mean intensity
    model.radiation.initialize_J();


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            //Approach which just accumulates the intensity contributions as the ray is traced

            solve_shortchar_order_0<approx> (model, o, rr);
            solve_shortchar_order_0<approx> (model, o, ar);

            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
                model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
            }

            //Approach which first traces a ray

            // const Real dshift_max = get_dshift_max (model, o);
            //
            // nr_   ()[centre] = o;
            // shift_()[centre] = 1.0;
            //
            // first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            // last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            // n_tot_() = (last_()+1) - first_();
            //
            // if (n_tot_() > 1)
            // {
            //     solve_shortchar_order_0_ray_forward (model, o, rr);
            //     solve_shortchar_order_0_ray_backward (model, o, ar);
            //
            //     for (Size f = 0; f < model.parameters->nfreqs(); f++)
            //     {
            //         model.radiation.u(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) + model.radiation.I(ar,o,f));
            //         model.radiation.v(rr,o,f) = 0.5 * (model.radiation.I(rr,o,f) - model.radiation.I(ar,o,f));
            //     }
            // }
            // else
            // {
            //     for (Size f = 0; f < model.parameters->nfreqs(); f++)
            //     {
            //         model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
            //         model.radiation.v(rr,o,f)  = 0.0;
            //     }
            // }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.I.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();
}


template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2_uv (Model& model)
{
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities)
    {
        model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.v.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
    }


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            const Real dshift_max = get_dshift_max (model, o);

            nr_   ()[centre] = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

            if (n_tot_() > 1)
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    solve_feautrier_order_2_uv <approx> (model, o, f);

                    model.radiation.u(rr,o,f)  = Su_()[centre];
                    model.radiation.v(rr,o,f)  = Sv_()[centre];
                }
            }
            else
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    model.radiation.u(rr,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.v(rr,o,f)  = 0.0;
                }
            }
        })

        pc::accelerator::synchronize();
    }

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.v.copy_ptr_to_vec();
}


template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2_sparse (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J(o,k) = 0.0;
            }
        })
    }


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        const Real     wt = model.geometry.rays.weight   [rr] * two;
        const Vector3D nn = model.geometry.rays.direction[rr];

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
        {
            threaded_for (o, model.parameters->npoints(),
            {
                const Real dshift_max = get_dshift_max (model, o);

                nr_   ()[centre] = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
                last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
                n_tot_() = (last_()+1) - first_();

                if (n_tot_() > 1)
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            solve_feautrier_order_2 <approx> (model, o,  lspec.nr_line[o][k][z]);

                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * Su_()[centre];

                            update_Lambda<approx> (model, rr, lspec.nr_line[o][k][z]);
                        }
                    }
                }
                else
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                        }
                    }
                }
            })
        }
    }
}


template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2_anis (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.J      .resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_0   .resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_1_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Re.resize(model.parameters->npoints(), lspec.linedata.nrad);
        lspec.J2_2_Im.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J       (o,k) = 0.0;
                lspec.J2_0    (o,k) = 0.0;
                lspec.J2_1_Re (o,k) = 0.0;
                lspec.J2_1_Im (o,k) = 0.0;
                lspec.J2_2_Re (o,k) = 0.0;
                lspec.J2_2_Im (o,k) = 0.0;
            }
        })
    }


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        const Real     wt = model.geometry.rays.weight   [rr];
        const Vector3D nn = model.geometry.rays.direction[rr];

        const Real wt_0    =     inv_sqrt2 * (three * nn.z() * nn.z() - one);
        const Real wt_1_Re =        -sqrt3 *          nn.x() * nn.z();
        const Real wt_1_Im =        -sqrt3 *          nn.y() * nn.z();
        const Real wt_2_Re =  half * sqrt3 * (nn.x() * nn.x() - nn.y() * nn.y());
        const Real wt_2_Im =         sqrt3 *          nn.x() * nn.y();

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
        {
            threaded_for (o, model.parameters->npoints(),
            {
                const Real dshift_max = get_dshift_max (model, o);

                nr_   ()[centre] = o;
                shift_()[centre] = 1.0;

                first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
                last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
                n_tot_() = (last_()+1) - first_();

                if (n_tot_() > 1)
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            solve_feautrier_order_2 <approx> (model, o, lspec.nr_line[o][k][z]);

                            const Real du = lspec.quadrature.weights[z] * wt * Su_()[centre];

                            lspec.J       (o,k) += two     * du;
                            lspec.J2_0    (o,k) += wt_0    * du;
                            lspec.J2_1_Re (o,k) += wt_1_Re * du;
                            lspec.J2_1_Im (o,k) += wt_1_Im * du;
                            lspec.J2_2_Re (o,k) += wt_2_Re * du;
                            lspec.J2_2_Im (o,k) += wt_2_Im * du;
                        }
                    }
                }
                else
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            const Real du = lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));

                            lspec.J       (o,k) += two     * du;
                            lspec.J2_0    (o,k) += wt_0    * du;
                            lspec.J2_1_Re (o,k) += wt_1_Re * du;
                            lspec.J2_1_Im (o,k) += wt_1_Im * du;
                            lspec.J2_2_Re (o,k) += wt_2_Re * du;
                            lspec.J2_2_Im (o,k) += wt_2_Im * du;
                        }
                    }
                }
            })
        }
    }
}


//sparse Feautrier solver, but now does not put point without close lines on the ray
template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2_sparse_pruned_rays (Model& model)
{
    // Initialise variables
    for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
    {
        lspec.lambda.clear();

        lspec.J.resize(model.parameters->npoints(), lspec.linedata.nrad);

        threaded_for (o, model.parameters->npoints(),
        {
            for (Size k = 0; k < lspec.linedata.nrad; k++)
            {
                lspec.J(o,k) = 0.0;
            }
        })
    }


    // For each ray, solve transfer equation
    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size     ar = model.geometry.rays.antipod  [rr];
        const Real     wt = model.geometry.rays.weight   [rr] * two;
        const Vector3D nn = model.geometry.rays.direction[rr];

        cout << "--- rr = " << rr << endl;

        for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
        {
            threaded_for (o, model.parameters->npoints(),
            {
                const Real dshift_max = get_dshift_max (model, o);


                // first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
                // last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
                // n_tot_() = (last_()+1) - first_();

                //Basic ray tracing to check if ray length is at least 2 (so no direct boundary condition required)
                //By definition of the comoving frame, o defines the reference frame. To correctly compute the intensity,
                // we will at least include the origin of the ray and the neighbors in the pruning process.
                const Size tot_ray_length = 1 + model.geometry.get_ray_length <CoMoving>(o, rr, dshift_max) + model.geometry.get_ray_length <CoMoving>(o, ar, dshift_max);

                // if (n_tot_() > 1)
                if (tot_ray_length > 1)
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        const Real line_frequency = lspec.linedata.frequency[k];
                        //Per line, trace which part we actually need to solve
                        //Is identical for all pruned rays
                        nr_   ()[centre] = o;
                        shift_()[centre] = 1.0;

                        first_() = trace_ray_pruned <CoMoving> (model, o, rr, dshift_max, -1, centre-1, centre-1, line_frequency) + 1;
                        last_ () = trace_ray_pruned <CoMoving> (model, o, ar, dshift_max, +1, centre+1, centre  , line_frequency) - 1;
                        n_tot_() = (last_()+1) - first_();
                        // std::cout<<"ntot: "<<n_tot_()<<std::endl;

                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            solve_feautrier_order_2 <approx> (model, o,  lspec.nr_line[o][k][z]);

                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * Su_()[centre];

                            update_Lambda<approx> (model, rr, lspec.nr_line[o][k][z]);
                        }
                    }
                }
                else
                {
                    for (Size k = 0; k < lspec.linedata.nrad; k++)
                    {
                        // Integrate over the line
                        for (Size z = 0; z < model.parameters->nquads(); z++)
                        {
                            lspec.J(o,k) += lspec.quadrature.weights[z] * wt * boundary_intensity(model, o, model.radiation.frequencies.nu(o, lspec.nr_line[o][k][z]));
                        }
                    }
                }
            })
        }
    }
}





template<ApproximationType approx>
inline void Solver :: solve_feautrier_order_2 (Model& model)
{
    // Allocate memory if not pre-allocated
    if (!model.parameters->store_intensities)
    {
        model.radiation.u.resize (model.parameters->hnrays(), model.parameters->npoints(), model.parameters->nfreqs());
        model.radiation.J.resize (                           model.parameters->npoints(), model.parameters->nfreqs());
    }

    // Initialise Lambda operator
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    // Initialise mean intensity
    model.radiation.initialize_J();

    // For each ray, solve transfer equation
    distributed_for (rr, rr_loc, model.parameters->hnrays(),
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            const Real dshift_max = get_dshift_max (model, o);

            nr_   ()[centre] = o;
            shift_()[centre] = 1.0;

            first_() = trace_ray <CoMoving> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
            last_ () = trace_ray <CoMoving> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
            n_tot_() = (last_()+1) - first_();

            if (n_tot_() > 1)
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    solve_feautrier_order_2 <approx> (model, o, f);

                    model.radiation.u(rr_loc,o,f)  = Su_()[centre];
                    model.radiation.J(       o,f) += Su_()[centre] * two * model.geometry.rays.weight[rr];

                    update_Lambda<approx> (model, rr, f);
                }
            }
            else
            {
                for (Size f = 0; f < model.parameters->nfreqs(); f++)
                {
                    model.radiation.u(rr_loc,o,f)  = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
                    model.radiation.J(       o,f) += two * model.geometry.rays.weight[rr] * model.radiation.u(rr_loc,o,f);
                }
            }
        })

        pc::accelerator::synchronize();
    })

    model.radiation.u.copy_ptr_to_vec();
    model.radiation.J.copy_ptr_to_vec();

    cout << "MPI gathering J..." << endl;
    model.radiation.MPI_reduce_J();
    cout << "Done MPI gathering J." << endl;

    cout << "MPI gathering Lambda..." << endl;
    // Gather contributions to the ALO
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.MPI_gather();}
    cout << "Done MPI gathering Lambda." << endl;
}

template<ApproximationType approx>
inline void Solver :: image_feautrier_order_2 (Model& model, const Size rr)
{
    Image image = Image(model.geometry, Intensity, rr);

    const Size ar = model.geometry.rays.antipod[rr];

    accelerated_for (o, model.parameters->npoints(),
    {
        const Real dshift_max = get_dshift_max (model, o);

        nr_   ()[centre] = o;
        shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

        first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
        last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
        n_tot_() = (last_()+1) - first_();

        if (n_tot_() > 1)
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image_feautrier_order_2<approx> (model, o, f);

                image.I(o,f) = two*Su_()[last_()] - boundary_intensity(model, nr_()[last_()], model.radiation.frequencies.nu(o, f));
            }
        }
        else
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image.I(o,f) = boundary_intensity(model, o, model.radiation.frequencies.nu(o, f));
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back (image);
}

template<ApproximationType approx>
inline void Solver :: image_feautrier_order_2_for_point (Model& model, const Size rr, const Size p)
{
    // Redefine p to keep everything as similar to image_feautrier_order_2 as possible
    const Size o = p;


    const Size ar = model.geometry.rays.antipod[rr];

    const Real dshift_max = get_dshift_max (model, o);

    nr_   ()[centre] = o;
    shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

    first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
    last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
    n_tot_() = (last_()+1) - first_();

    model.   S_ray.resize (n_tot_(), model.parameters->nfreqs());
    model.dtau_ray.resize (n_tot_(), model.parameters->nfreqs());
    model.   u_ray.resize (n_tot_(), model.parameters->nfreqs());

    if (n_tot_() > 1)
    {
        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            image_feautrier_order_2_for_point_loc<approx> (model, o, f);
        }
    }

}


template<ApproximationType approx>
inline void Solver :: image_optical_depth (Model& model, const Size rr)
{
    Image image = Image(model.geometry, OpticalDepth, rr);

    const Size ar = model.geometry.rays.antipod[rr];

    accelerated_for (o, model.parameters->npoints(),
    {
        const Real dshift_max = get_dshift_max (model, o);

        nr_   ()[centre] = o;
        shift_()[centre] = model.geometry.get_shift <Rest> (o, rr, o, 0.0);;

        first_() = trace_ray <Rest> (model.geometry, o, rr, dshift_max, -1, centre-1, centre-1) + 1;
        last_ () = trace_ray <Rest> (model.geometry, o, ar, dshift_max, +1, centre+1, centre  ) - 1;
        n_tot_() = (last_()+1) - first_();

        if (n_tot_() > 1)
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image_optical_depth<approx> (model, o, f);

                image.I(o,f) = optical_depth_();
            }
        }
        else
        {
            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                image.I(o,f) = 0.0;
            }
        }
    })

    pc::accelerator::synchronize();

    model.images.push_back (image);
}

//Because of the new method for computing the optical depth, adding extra frequency points for counteracting the large doppler shift is no longer necessary
template <Frame frame>
accel inline Size Solver :: trace_ray (
    const Geometry& geometry,
    const Size      o,
    const Size      r,
    const double    dshift_max,
    const int       increment,
          Size      id1,
          Size      id2 )
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = geometry.get_next (o, r, o, Z, dZ);

    if (geometry.valid_point(nxt))
    {
        Size         crt = o;
        double shift_crt = geometry.get_shift <frame> (o, r, crt, 0.0);
        double shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z  );

        set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = geometry.get_next          (o, r, nxt, Z, dZ);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z    );

            set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}

//Because of the new method for computing the optical depth, adding extra frequency points for counteracting the large doppler shift is no longer necessary
template <Frame frame>
accel inline Size Solver :: trace_ray_comoving (
    const Geometry& geometry,
    const Size      o,
    const Size      r,
    const Size      rr,
    const Size      rayidx,
    const double    dshift_max,
    const int       increment,
          Size      id1,
          Size      id2,
          Size&     outermost_interesting_point_rayidx)
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = geometry.get_next (o, r, o, Z, dZ);

    if (geometry.valid_point(nxt))
    {
        Size         crt = o;
        double shift_crt = geometry.get_shift <frame> (o, r, crt, 0.0);
        double shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z  );

        if (closest_ray(rr, nxt)==rayidx)
        {
            outermost_interesting_point_rayidx = id1;//as the data is set there, this should? be fine
        }

        set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);

        while (geometry.not_on_boundary(nxt))
        {
                  crt =       nxt;
            shift_crt = shift_nxt;

                  nxt = geometry.get_next          (o, r, nxt, Z, dZ);
            shift_nxt = geometry.get_shift <frame> (o, r, nxt, Z    );

            if (closest_ray(rr, nxt)==rayidx)
            {
                outermost_interesting_point_rayidx = id1;//as the data is set there, this should? be fine
            }

            set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
        }
        // //has to be a valid point in order to evaluate this, so inside these brackets
        // if (closest_ray(rr, nxt)==rayidx)
        // {   //subtracting increment required, as id1 always points to the next point after set_data
        //     outermost_interesting_point_rayidx = id1-increment;//as the data is set there, this should? be fine
        // }
    }

    return id1;
}

//Tracing the ray, ignoring all points for which the given CMF frequency lies too far from all line centers
template <Frame frame>
accel inline Size Solver :: trace_ray_pruned (
    const Model&    model,
    const Size      o,
    const Size      r,
    const double    dshift_max,
    const int       increment,
          Size      id1,
          Size      id2,
    const Real      freq )
{
    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z
    double dZcrt = dZ; //previous increment in dZ

    bool is_crt_set = true; //Whether the current point is already set

    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt))
    {
        Size         crt = o;
        // Size         prv = crt;
        double shift_crt = model.geometry.get_shift <frame> (o, r, crt, 0.0);
        // double shift_prv = shift_crt;
        double shift_nxt = model.geometry.get_shift <frame> (o, r, nxt, Z  );

        //As we must make sure that all lines are traced fine
        //check whether position increment has any close lines
        if (check_close_line(shift_crt*freq, shift_nxt*freq, crt, nxt, model))
        // if (check_close_line(shift_prv*freq, shift_crt*freq, shift_nxt*freq, prv, crt, nxt, model))
        {
            set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
            is_crt_set = true;
        }
        else
        {
            is_crt_set=false;
        }

        while (model.geometry.not_on_boundary(nxt))
        {
            //       prv =       crt;
            // shift_prv = shift_crt;
                  crt =       nxt;
            shift_crt = shift_nxt;

            dZcrt = dZ;//current dZ required

                  nxt = model.geometry.get_next          (o, r, nxt, Z, dZ);
            shift_nxt = model.geometry.get_shift <frame> (o, r, nxt, Z    );

            //As we must make sure that all lines are traced fine
            //check whether position increment has any close lines
            if (check_close_line(shift_crt*freq, shift_nxt*freq, crt, nxt, model))
            // if (check_close_line(shift_prv*freq, shift_crt*freq, shift_nxt*freq, prv, crt, nxt, model))
            {
                if (!is_crt_set)
                {
                    //FIXME: FIX FUNCTION DEFINITION SET_DATA
                    set_data (crt, crt, shift_crt, shift_crt, dZcrt, dshift_max, increment, id1, id2);
                }
                set_data (crt, nxt, shift_crt, shift_nxt, dZ, dshift_max, increment, id1, id2);
                is_crt_set = true;
            }
            else
            {
                is_crt_set = false;
            }
        }

        //always set the boundary point (otherwise the boundary condition cannot be evaluated)
        if (!is_crt_set)
        {
            set_data (crt, crt, shift_crt, shift_crt, dZcrt, dshift_max, increment, id1, id2);
        }
    }

    return id1;
}

accel inline bool Solver :: check_close_line (const Real currfreq, const Real nextfreq, const Size currpoint, const Size nextpoint, const Model& model)
// accel inline bool Solver :: check_close_line (const Real prevfreq, const Real currfreq, const Real nextfreq, const Size prevpoint, const Size currpoint, const Size nextpoint, const Model& model)
{
    const Real left_freq = std::min({currfreq, nextfreq});
    const Real right_freq = std::max({currfreq, nextfreq});

    //FIXME?: use wider bounds, as this is just for quantifying whether the line center lies close enough. mhh, the farther we lie from the line center, the less influence the evaluated frequency has on the mean line intensity...
    //using maximum of bounds on the two points to get an upper bound for the line width
    // const Real prev_bound_line_width = model.thermodynamics.profile_width_upper_bound_with_linefreq(prevpoint, right_freq, model.lines.max_inverse_mass);
    const Real curr_bound_line_width = model.thermodynamics.profile_width_upper_bound_with_linefreq(currpoint, right_freq, model.lines.max_inverse_mass);
    const Real next_bound_line_width = model.thermodynamics.profile_width_upper_bound_with_linefreq(nextpoint, right_freq, model.lines.max_inverse_mass);
    const Real upper_bound_line_width = model.parameters->max_distance_opacity_contribution * std::max({curr_bound_line_width, next_bound_line_width});

    const Real left_freq_bound = left_freq - upper_bound_line_width;
    const Real right_freq_bound = right_freq + upper_bound_line_width;

    //apply default search algorithms on the bounds, obtaining iterators
    auto left_line_bound=std::lower_bound(model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound=std::upper_bound(model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);
    //yes, I am comparing pointers; I just want to know whether at least one line lies in the interval // see compute_S_dtau_line_integrated <CloseLines>
    //TODO: refactor such that both trace_ray_pruned and compute_S_dtau_line_integrated <CloseLines> rely on the same code snippet (instead of duplicating this code)
    return (left_line_bound!=right_line_bound);
}


//Because of the new method for computing the optical depth, adding extra frequency points for counteracting the large doppler shift is no longer necessary
accel inline void Solver :: set_data (
    const Size   crt,
    const Size   nxt,
    const double shift_crt,
    const double shift_nxt,
    const double dZ_loc,
    const double dshift_max,
    const int    increment,
          Size&  id1,
          Size&  id2 )
{
    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    const double dshift     = shift_nxt - shift_crt;
    const double dshift_abs = fabs (dshift);

    nr   [id1] = nxt;
    shift[id1] = shift_nxt;
    dZ   [id2] = dZ_loc;

    id1 += increment;
    id2 += increment;
}

///  Gaussian line profile function
///    @param[in] width : profile width
///    @param[in] diff  : frequency difference with line centre
///    @return profile function evaluated with this frequency difference
////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: gaussian (const Real inverse_width, const Real diff) const
{
    const Real sqrt_exp = inverse_width * diff;

    return inverse_width * INVERSE_SQRT_PI * expf (-sqrt_exp*sqrt_exp);
}


///  Planck function
///    @param[in] temp : temperature of the corresponding black body
///    @param[in] freq : frequency at which to evaluate the function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: planck (const Real temp, const Real freq) const
{
    return TWO_HH_OVER_CC_SQUARED * (freq*freq*freq) / expm1f (HH_OVER_KB*freq/temp);
}


///  Getter for the boundary conditions
///    @param[in] model  : reference to model object
///    @param[in] p      : point index of the boundary point
///    @param[in] freq   : frequency at which to evaluate boundary condition
///    @returns incoming radiation intensity at the boundary
////////////////////////////////////////////////////////////////////////////
accel inline Real Solver :: boundary_intensity (const Model& model, const Size p, const Real freq) const
{
    const Size bdy_id = model.geometry.boundary.point2boundary[p];

    switch (model.geometry.boundary.boundary_condition[bdy_id])
    {
        case Zero    : return 0.0;
        case Thermal : return planck (model.geometry.boundary.boundary_temperature[bdy_id], freq);
        default      : return planck (T_CMB, freq);
    }
}


///  Getter for the emissivity (eta) and the opacity (chi)
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
template<>
accel inline void Solver :: get_eta_and_chi <None> (
    const Model& model,
    const Size   p,
    const Size   ll,  // dummy variable
    const Real   freq,
          Real&  eta,
          Real&  chi ) const
{
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    // Set line emissivity and opacity
    for (Size l = 0; l < model.parameters->nlines(); l++)
    {
        const Real diff = freq - model.lines.line[l];
        const Real prof = freq * gaussian (model.lines.inverse_width(p, l), diff);

        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity   (p, l);
    }
}


///  Getter for the emissivity (eta) and the opacity (chi)
///  function uses only the nearby lines to save computation time
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////
template<>
accel inline void Solver :: get_eta_and_chi <CloseLines> (
    const Model& model,
    const Size   p,
    const Size   ll,  // dummy variable
    const Real   freq,
          Real&  eta,
          Real&  chi ) const
{
    // Initialize
    eta = 0.0;
    chi = model.parameters->min_opacity;

    const Real upper_bound_line_width = model.parameters->max_distance_opacity_contribution * model.thermodynamics.profile_width_upper_bound_with_linefreq(p, freq, model.lines.max_inverse_mass);
    const Real left_freq_bound = freq - upper_bound_line_width;
    const Real right_freq_bound = freq + upper_bound_line_width;

    //Just using default search algorithms, obtaining iterators
    auto left_line_bound=std::lower_bound(model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound=std::upper_bound(model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);


    for (auto freq_sort_l = left_line_bound; freq_sort_l != right_line_bound; freq_sort_l++)
    {
        const Size sort_l = freq_sort_l-model.lines.sorted_line.begin();
        // mapping sorted line index to original line index
        const Size l = model.lines.sorted_line_map[sort_l];
        // const Real diff = freq - model.lines.line[l];
        const Real diff = freq - *freq_sort_l; //should be equal to the previous line of code
        const Real inv_width = model.lines.inverse_width(p, l);
        const Real prof = freq * gaussian (model.lines.inverse_width(p, l), diff);
        eta += prof * model.lines.emissivity(p, l);
        chi += prof * model.lines.opacity   (p, l);
    }
}


///  Getter for the emissivity (eta) and the opacity (chi) in the "one line" approixmation
///    @param[in]  model : reference to model object
///    @param[in]  p     : index of the point
///    @param[in]  l     : line index corresponding to the frequency
///    @param[in]  freq  : frequency (in co-moving frame)
///    @param[out] eta   : emissivity
///    @param[out] chi   : opacity
//////////////////////////////////////////////////////////////////////////////////////////
template<>
accel inline void Solver :: get_eta_and_chi <OneLine> (
    const Model& model,
    const Size   p,
    const Size   l,
    const Real   freq,
          Real&  eta,
          Real&  chi ) const
{
    const Real diff = freq - model.lines.line[l];
    const Real prof = freq * gaussian (model.lines.inverse_width(p, l), diff);

    eta = prof * model.lines.emissivity(p, l);
    chi = prof * model.lines.opacity   (p, l) + model.parameters->min_opacity;
}


///  Apply trapezium rule to x_crt and x_nxt
///    @param[in] x_crt : current value of x
///    @param[in] x_nxt : next value of x
///    @param[in] dZ    : distance inscrement along ray
///    @returns integral x over dZ
///////////////////////////////////////////////////////
accel inline Real trap (const Real x_crt, const Real x_nxt, const double dZ)
{
    return half * (x_crt + x_nxt) * dZ;
}




template<ApproximationType approx>
accel inline void Solver :: solve_shortchar_order_0 (
          Model& model,
    const Size   o,
    const Size   r)
{
    Vector<Real>& eta_c = eta_c_();
    Vector<Real>& eta_n = eta_n_();

    Vector<Real>& chi_c = chi_c_();
    Vector<Real>& chi_n = chi_n_();

    Vector<Real>& source_c= source_c_();
    Vector<Real>& source_n= source_n_();

    Vector<Real>& tau = tau_();


    double  Z = 0.0;   // distance along ray
    double dZ = 0.0;   // last distance increment

    Size crt = o;
    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);
    Real term_c, term_n, dtau;
    bool compute_curr_opacity, prev_compute_curr_opacity;

    if (model.geometry.valid_point (nxt))
    {
        double shift_c = 1.0;
        double shift_n = model.geometry.get_shift <CoMoving> (o, r, nxt, Z);

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);
            const Size l    = model.radiation.frequencies.corresponding_line[f];//line index

            compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

            compute_source_dtau<approx>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ, compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
            dtau = std::max(model.parameters->min_dtau, dtau);

            //proper implementation of 2nd order shortchar (not yet times reducing factor of exp(-tau))
            // model.radiation.I(r,o,f) = term_c * (expm1(-dtau)+dtau) / dtau
            //                          + term_n * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau;
            //Rewrite, trying to use less exponentials
            const Real factor = expm1f(-dtau)/dtau;

            model.radiation.I(r,o,f) = factor*(source_c[f]-source_n[f]*(1.0+dtau))
                                     + source_c[f] - source_n[f];
            tau[f] = dtau;

            //Compute local lambda operator
            const Size l_spec = model.radiation.frequencies.corresponding_l_for_spec[f];   // index of species
            const Size k = model.radiation.frequencies.corresponding_k_for_tran[f];   // index of transition
            const Size z = model.radiation.frequencies.corresponding_z_for_line[f];   // index of quadrature point
            const Real w_ang = model.geometry.rays.weight[r];

            LineProducingSpecies &lspec = model.lines.lineProducingSpecies[l_spec];

            const Real freq_line = lspec.linedata.frequency[k];
            const Real invr_mass = lspec.linedata.inverse_mass;
            const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;

            Real eta, chi;//eta is dummy var
            //chi is not necessarily computed, so compute it to be sure
            get_eta_and_chi <approx>(model, o, k, freq_line, eta, chi);
            Real inverse_chi=1.0/chi;
            Real phi = model.thermodynamics.profile(invr_mass, o, freq_line, freq);
            // const Real lambda_factor = (dtau+expm1f(-dtau))/dtau;// If one wants to compute lambda a bit more accurately in case of dtau≃0.
            // Real L   = constante * freq * phi * lambda_factor * inverse_chi;
            Real L   = constante * freq * phi * (factor + 1.0) * inverse_chi;//using factor+1.0, the computed lambda elements can be negative if dtau very small; but then the lambda elements are also negligible
            lspec.lambda.add_element(o, k, o, L);

            //TODO: possible nonlocal lambda part // FIXME: probably incorrect chi used
            // L   = constante * freq * phi * (-factor * (1.0+dtau) - 1.0) * inverse_chi;
            // lspec.lambda.add_element(o, k, nxt, L);
        }

        //For all frequencies, we need to use the same method for computing the optical depth
        // bool prev_compute_curr_opacity=compute_curr_opacity;//technically, we could also keep this bool individually for every frequency
        prev_compute_curr_opacity=compute_curr_opacity;//technically, we could also keep this bool individually for every frequency

        while (model.geometry.not_on_boundary (nxt))
        {
            crt     = nxt;
            shift_c = shift_n;

            model.geometry.get_next (o, r, crt, nxt, Z, dZ, shift_n);

            for (Size f = 0; f < model.parameters->nfreqs(); f++)
            {
                source_c[f]=source_n[f];
                chi_c[f]=chi_n[f];
                const Real freq = model.radiation.frequencies.nu(o, f);
                const Size l    = model.radiation.frequencies.corresponding_line[f];

                compute_curr_opacity=prev_compute_curr_opacity;

                compute_source_dtau<approx>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ, compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
                dtau = std::max(model.parameters->min_dtau, dtau);

                //proper implementation of 2nd order shortchar (not yet times reducing factor of exp(-tau))
                // model.radiation.I(r,o,f) += expf(-tau[f]) *
                //                          ( term_c * (expm1(-dtau)+dtau) / dtau
                //                          + term_n * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau);
                //Rewrite, trying to use less exponentials
                model.radiation.I(r,o,f) += expf(-tau[f]) *
                                           (expm1f(-dtau)/dtau*(source_c[f]-source_n[f]*(1.0+dtau))
                                           + source_c[f] - source_n[f]);
                //TODO: check order of addition, as we might be starting with the largest contributions, before adding the smaller ones...
                tau[f] += dtau;
            }

            //save setting for use for all frequencies for the next interval
            prev_compute_curr_opacity=compute_curr_opacity;
        }

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r,o,f) += boundary_intensity(model, nxt, freq*shift_n) * expf(-tau[f]);
            model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
        }
    }

    else
    {
        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(o, f);

            model.radiation.I(r,o,f)  = boundary_intensity(model, crt, freq);
            model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
        }
    }
}


// accel inline void Solver :: solve_shortchar_order_0_ray_forward (
//           Model& model,
//           const Size   o,
//           const Size   r)
// {
//     Vector<Real>& eta_c = eta_c_();
//     Vector<Real>& eta_n = eta_n_();
//
//     Vector<Real>& chi_c = chi_c_();
//     Vector<Real>& chi_n = chi_n_();
//
//     Vector<Real>& source_c= source_c_();
//     Vector<Real>& source_n= source_n_();
//
//     Vector<Real>& tau = tau_();
//     Vector<double>& shift=shift_();
//     Vector<Size>& nr=nr_();
//     Vector<double>& dZ=dZ_();
//
//     Size crt, nxt;
//     // Size nxt = model.geometry.get_next (o, r, o, Z, dZ);
//     Real term_c, term_n, dtau;
//     bool compute_curr_opacity, prev_compute_curr_opacity;
//     prev_compute_curr_opacity=true;
//
//     // Set boundary condition
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         const Real freq = model.radiation.frequencies.nu(o, f);
//         model.radiation.I(r,o,f)=boundary_intensity(model, nr[first_()], freq*shift[first_()]);
//     }
//     double shift_c, shift_n;
//     //interate until we reach the middle point
//     for (Size idx=first_()+1; idx<=centre; idx++)
//     {
//
//         crt = nr[idx-1];
//         nxt = nr[idx];
//         shift_c = shift[idx-1];
//         shift_n = shift[idx];
//
//         for (Size f = 0; f < model.parameters->nfreqs(); f++)
//         {
//             chi_c[f]=chi_n[f];
//             source_c[f]=source_n[f];
//
//             const Real freq = model.radiation.frequencies.nu(o, f);
//             const Size l    = model.radiation.frequencies.corresponding_line[f];
//
//             compute_curr_opacity = prev_compute_curr_opacity; // for the first point, we need to compute both the curr and next opacity (and source)
//             // compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)
//
//             compute_source_dtau<None>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ[idx-1], compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
//             dtau = std::max(model.parameters->min_dtau, dtau);
//
//             // model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//             //                          + term_c * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau
//             //                          + term_n * (expm1(-dtau)+dtau) / dtau;
//             //slight rewrite, should be more stable
//             model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//                                      + ( source_c[f] * (-expm1(-dtau)-dtau*exp(-dtau))
//                                        + source_n[f] * (expm1(-dtau)+dtau) )/ dtau;
//         }
//         //save setting for use for all frequencies for the next interval
//         prev_compute_curr_opacity=compute_curr_opacity;
//     }
//
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
//     }
// }
//
//
// accel inline void Solver :: solve_shortchar_order_0_ray_backward (
//           Model& model,
//           const Size   o,
//           const Size r)
// {
//     Vector<Real>& eta_c = eta_c_();
//     Vector<Real>& eta_n = eta_n_();
//
//     Vector<Real>& chi_c = chi_c_();
//     Vector<Real>& chi_n = chi_n_();
//
//     Vector<Real>& source_c= source_c_();
//     Vector<Real>& source_n= source_n_();
//
//     Vector<Real>& tau = tau_();
//     Vector<double>& shift=shift_();
//     Vector<Size>& nr=nr_();
//     Vector<double>& dZ=dZ_();
//
//     Size crt, nxt;
//     Real term_c, term_n, dtau;
//     bool compute_curr_opacity, prev_compute_curr_opacity;
//     prev_compute_curr_opacity=true; // for the first point, we need to compute both the curr and next opacity (and source)
//
//     // Set boundary condition
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         const Real freq = model.radiation.frequencies.nu(o, f);
//         model.radiation.I(r,o,f)=boundary_intensity(model, nr[last_()], freq*(shift[last_()]));
//     }
//     double shift_c, shift_n;
//     //interate until we reach the middle point
//     for (Size indexp1=last_(); indexp1>=centre+1; indexp1--)
//     {
//         crt = nr[indexp1];
//         nxt = nr[indexp1-1];
//         shift_c = shift[indexp1];
//         shift_n = shift[indexp1-1];
//
//         for (Size f = 0; f < model.parameters->nfreqs(); f++)
//         {
//             chi_c[f]=chi_n[f];
//             source_c[f]=source_n[f];
//
//             const Real freq = model.radiation.frequencies.nu(o, f);
//             const Size l    = model.radiation.frequencies.corresponding_line[f];
//
//             compute_curr_opacity = prev_compute_curr_opacity; // for the first point, we need to compute both the curr and next opacity (and source)
//             // compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)
//
//             compute_source_dtau<None>(model, crt, nxt, l, freq*shift_c, freq*shift_n, shift_c, shift_n, dZ[indexp1-1], compute_curr_opacity, dtau, chi_c[f], chi_n[f], source_c[f], source_n[f]);
//             dtau = std::max(model.parameters->min_dtau, dtau);
//
//             // model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//             //                          + term_c * (-expm1(-dtau)-dtau*expf(-dtau)) /dtau
//             //                          + term_n * (expm1(-dtau)+dtau) / dtau;
//             //slight rewrite, should be more stable
//             model.radiation.I(r,o,f) = exp(-dtau)*model.radiation.I(r,o,f)
//                                      + ( source_c[f] * (-expm1(-dtau)-dtau*exp(-dtau))
//                                        + source_n[f] * (expm1(-dtau)+dtau) )/ dtau;
//         }
//
//         //save setting for use for all frequencies for the next interval
//         prev_compute_curr_opacity=compute_curr_opacity;
//     }
//
//     for (Size f = 0; f < model.parameters->nfreqs(); f++)
//     {
//         model.radiation.J(  o,f) += model.geometry.rays.weight[r] * model.radiation.I(r,o,f);
//     }
// }


template<ApproximationType approx>
accel inline void Solver :: update_Lambda (Model &model, const Size rr, const Size f)
{
    const Frequencies    &freqs     = model.radiation.frequencies;
    const Thermodynamics &thermodyn = model.thermodynamics;

    if (freqs.appears_in_line_integral[f])
    {
        const Size first = first_();
        const Size last  = last_ ();
        const Size n_tot = n_tot_();

        Vector<Size  >& nr          = nr_         ();
        Vector<double>& shift       = shift_      ();
        Vector<Real  >& L_diag      = L_diag_     ();
        Matrix<Real  >& L_upper     = L_upper_    ();
        Matrix<Real  >& L_lower     = L_lower_    ();
        // Vector<Real  >& inverse_chi = inverse_chi_();

        const Real w_ang = two * model.geometry.rays.weight[rr];

        const Size l = freqs.corresponding_l_for_spec[f];   // index of species
        const Size k = freqs.corresponding_k_for_tran[f];   // index of transition
        const Size z = freqs.corresponding_z_for_line[f];   // index of quadrature point

        LineProducingSpecies &lspec = model.lines.lineProducingSpecies[l];

        const Real freq_line = lspec.linedata.frequency[k];
        const Real invr_mass = lspec.linedata.inverse_mass;
        const Real constante = lspec.linedata.A[k] * lspec.quadrature.weights[z] * w_ang;
        //TODO: approximation if we do not have overlapping lines
        // Real inverse_opacity = 1.0/model.lines.opacity (nr[centre], k); //Inverse line opacity; includes 1/HH_OVER_FOUR_PI
        // Real L   = constante * L_diag[centre] * inverse_opacity;
        Real eta, chi;//eta is dummy var
        get_eta_and_chi <approx>(model, nr[centre], k, freq_line, eta, chi);
        Real inverse_chi=1.0/chi;

        Real frq = freqs.nu(nr[centre], f) * shift[centre];
        Real phi = thermodyn.profile(invr_mass, nr[centre], freq_line, frq);
        Real L   = constante * frq * phi * L_diag[centre] * inverse_chi;

        lspec.lambda.add_element(nr[centre], k, nr[centre], L);

        for (long m = 0; (m < n_off_diag) && (m+1 < n_tot); m++)
        {
            if (centre >= first+m+1) // centre-m-1 >= first
            {
                const long n = centre-m-1;

                //TODO: approximation if we do not have overlapping lines // also check which opacity we need (centre?)
                // inverse_opacity = 1.0/model.lines.opacity (nr[n], k); //Inverse line opacity; no longer 1/HH_OVER_FOUR_PI
                // L   = constante * L_lower(m,n) * inverse_opacity;
                get_eta_and_chi <approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi=1.0/chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_lower(m,n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }

            if (centre+m+1 <= last) // centre+m+1 < last
            {
                const long n = centre+m+1;

                //TODO: approximation if we do not have overlapping lines
                // inverse_opacity = 1.0/model.lines.opacity (nr[n], k); //Inverse line opacity; no includes 1/HH_OVER_FOUR_PI
                // L   = constante * L_upper(m,n) * inverse_opacity;
                get_eta_and_chi <approx>(model, nr[n], k, freq_line, eta, chi);
                Real inverse_chi=1.0/chi;

                frq = freqs.nu(nr[n], f) * shift[n];
                phi = thermodyn.profile (invr_mass, nr[n], freq_line, frq);
                L   = constante * frq * phi * L_upper(m,n) * inverse_chi;

                lspec.lambda.add_element(nr[centre], k, nr[n], L);
            }
        }
    }
}

///  Computer for the optical depth and source function when computing using the formal line integration
///  In case of low velocity differences, almost two times slower
///  For high velocity increments however, this does not need any extra interpolation points (-> way faster) (and some extra optimizations are made in the limit of extremely large doppler shift (as erf goes to its limit value)
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in comoving frame)
///    @param[in] nextfreq : frequency at next point (in comoving frame)
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point to compute
///    @param[out] Snext : source function at next point to compute
/////////////////////////////////////////////////////////////////////
template<>
inline void Solver :: compute_S_dtau_line_integrated <OneLine> (Model& model, Size currpoint, Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext)
{
    dtau=compute_dtau_single_line(model, currpoint, nextpoint, lineidx, currfreq, nextfreq, dZ);
    Scurr=model.lines.emissivity(currpoint, lineidx)/(model.lines.opacity(currpoint, lineidx)+model.parameters->min_opacity);//current source
    Snext=model.lines.emissivity(nextpoint, lineidx)/(model.lines.opacity(nextpoint, lineidx)+model.parameters->min_opacity);//next source
    //note: due to interaction with dtau when computing all sources individually, we do need to recompute Scurr and Snext for all position increments
}

///  Computer for the optical depth and source function when computing using the formal line integration
///  In case of low velocity differences, almost two times slower
///  For high velocity increments however, this does not need any extra interpolation points (-> way faster) (and some extra optimizations are made in the limit of extremely large doppler shift (as erf goes to its limit value)
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in comoving frame)
///    @param[in] nextfreq : frequency at next point (in comoving frame)
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point to compute
///    @param[out] Snext : source function at next point to compute
/////////////////////////////////////////////////////////////////////
template<>
inline void Solver :: compute_S_dtau_line_integrated <None> (Model& model, Size currpoint, Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext)
{
    Real sum_dtau=0.0;
    Real sum_dtau_times_Scurr=0.0;
    Real sum_dtau_times_Snext=0.0;
    for (Size l=0; l<model.parameters->nlines(); l++)
    {
        Real line_dtau=compute_dtau_single_line(model, currpoint, nextpoint, l, currfreq, nextfreq, dZ);
        Real line_Scurr=model.lines.emissivity(currpoint, l)/(model.lines.opacity(currpoint, l)+model.parameters->min_opacity);//current source
        Real line_Snext=model.lines.emissivity(nextpoint, l)/(model.lines.opacity(nextpoint, l)+model.parameters->min_opacity);//next source
        sum_dtau+=line_dtau;
        sum_dtau_times_Scurr+=line_dtau*line_Scurr;
        sum_dtau_times_Snext+=line_dtau*line_Snext;
    }
    dtau=sum_dtau;
    Scurr=sum_dtau_times_Scurr/sum_dtau;
    Snext=sum_dtau_times_Snext/sum_dtau;
}

///  Computer for the optical depth and source function when computing using the formal line integration
///  In case of low velocity differences, almost two times slower
///  For high velocity increments however, this does not need any extra interpolation points (-> way faster) (and some extra optimizations are made in the limit of extremely large doppler shift (as erf goes to its limit value)
///  This function only takes into account the nearby lines, saving some computation time
///    @param[in] curr_point : index of current point
///    @param[in] next_point : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in comoving frame)
///    @param[in] nextfreq : frequency at next point (in comoving frame)
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point to compute
///    @param[out] Snext : source function at next point to compute
/////////////////////////////////////////////////////////////////////
template<>
inline void Solver :: compute_S_dtau_line_integrated <CloseLines> (Model& model, Size currpoint, Size nextpoint, Size lineidx, Real currfreq, Real nextfreq, Real dZ, Real& dtau, Real& Scurr, Real& Snext)
{
    Real sum_dtau=0.0; //division by zero might occur otherwise
    // Real sum_dtau=model.parameters->min_dtau; //division by zero might occur otherwise
    Real sum_dtau_times_Scurr=0.0;
    Real sum_dtau_times_Snext=0.0;

    Real left_freq;
    Real right_freq;

    //err, compiler will probably figure out that I just want these two values ordered
    if (currfreq < nextfreq)
    {
        left_freq = currfreq;
        right_freq = nextfreq;
    }
    else
    {
        right_freq = currfreq;
        left_freq = nextfreq;
    }

    //using maximum of bounds on the two points to get an upper bound for the line width
    const Real curr_bound_line_width = model.parameters->max_distance_opacity_contribution * model.thermodynamics.profile_width_upper_bound_with_linefreq(currpoint, right_freq, model.lines.max_inverse_mass);
    const Real next_bound_line_width = model.parameters->max_distance_opacity_contribution * model.thermodynamics.profile_width_upper_bound_with_linefreq(nextpoint, right_freq, model.lines.max_inverse_mass);
    const Real upper_bound_line_width = std::max(curr_bound_line_width, next_bound_line_width);

    const Real left_freq_bound = left_freq - upper_bound_line_width;
    const Real right_freq_bound = right_freq + upper_bound_line_width;

    //apply default search algorithms on the bounds, obtaining iterators
    auto left_line_bound=std::lower_bound(model.lines.sorted_line.begin(), model.lines.sorted_line.end(), left_freq_bound);
    auto right_line_bound=std::upper_bound(model.lines.sorted_line.begin(), model.lines.sorted_line.end(), right_freq_bound);

    for (auto freq_sort_l = left_line_bound; freq_sort_l != right_line_bound; freq_sort_l++)
    {
        const Size sort_l = freq_sort_l - model.lines.sorted_line.begin();
        // Map sorted line index to original line index
        const Size l = model.lines.sorted_line_map[sort_l];

        Real line_dtau=compute_dtau_single_line(model, currpoint, nextpoint, l, currfreq, nextfreq, dZ);
        Real line_Scurr=model.lines.emissivity(currpoint, l)/(model.lines.opacity(currpoint, l)+model.parameters->min_opacity);//current source
        Real line_Snext=model.lines.emissivity(nextpoint, l)/(model.lines.opacity(nextpoint, l)+model.parameters->min_opacity);//next source
        sum_dtau+=line_dtau;
        sum_dtau_times_Scurr+=line_dtau*line_Scurr;
        sum_dtau_times_Snext+=line_dtau*line_Snext;
    }
    dtau=sum_dtau;
    //needs extra bounding, as nothing may be added in the first place (above for loop may have looped over 0 elements)
    const Real bound_min_dtau = model.parameters->min_opacity * dZ;
    //Correct way of bounding from below; should be able to deal with very minor computation errors around 0.
    if (-bound_min_dtau<dtau)
    {
        dtau=std::max(bound_min_dtau, dtau);
    }
    //Note: 0 source functions can be returned if no lines are nearby; but then the negligible lower bound gets returned
    Scurr=sum_dtau_times_Scurr/dtau;
    Snext=sum_dtau_times_Snext/dtau;

    //note: due to interaction with dtau when computing all sources individually, we do need to recompute Scurr and Snext for all position increments
}

/// Computes the source function and optical depth in a hybrid manner
///    @param[in/out] compute_curr_opacity: for deciding whether we need to compute the current opacity when using the trapezoidal rule
///    @param[in] currpoint : index of current point
///    @param[in] nextpoint : index of next point
///    @param[in] lineidx : index of line to integrate over
///    @param[in] currfreq : frequency at current point (in comoving frame)
///    @param[in] nextfreq : frequency at next point (in comoving frame)
///    @param[in] currshift : shift at curr point
///    @param[in] nextshift : shift at next point
///    @param[in] dZ : position increment
///    @param[out] dtau : optical depth increment to compute
///    @param[out] Scurr : source function at current point to compute
///    @param[out] Snext : source function at next point to compute
/// Warning: depending on how large the doppler shift is, the opacity is NOT computed.
template<ApproximationType approx>
accel inline void Solver :: compute_source_dtau (Model& model, Size currpoint, Size nextpoint, Size line, Real curr_freq, Real next_freq, double curr_shift, double next_shift, Real dZ, bool& compute_curr_opacity, Real& dtaunext, Real& chicurr, Real& chinext, Real& Scurr, Real& Snext)
{
    //deciding which optical depth computation to use, depending on the doppler shift
    const double dshift     = next_shift-curr_shift;//shift[first+1]-shift[first];TODO
    const double dshift_abs = fabs (dshift);
    const double dshift_max = std::min(model.dshift_max[currpoint],model.dshift_max[nextpoint]);
    const bool using_large_shift = (dshift_abs > dshift_max);

    //fancy computation for large doppler shifts
    if (using_large_shift)
    {
        compute_curr_opacity=true;
        compute_S_dtau_line_integrated <approx> (model, currpoint, nextpoint, line, curr_freq, next_freq, dZ, dtaunext, Scurr, Snext);
        //OPACITY IS NOT COMPUTED IN THIS BRANCH!
    }
    else
    {
        //default computation using trapezoidal rule
        if (compute_curr_opacity)//fancy computation does not compute the current opacity, so we might need to recompute it here
        {
            compute_curr_opacity=false;
            Real eta_c=0.0;//current emissivity
            //also get previous opacity (emissivity does not matter)
            get_eta_and_chi <approx> (model, currpoint, line, curr_freq, eta_c, chicurr);
            Scurr=eta_c/chicurr;//might as well compute the source function too
        }

        Real eta_n=0.0;
        // Get new radiative properties
        get_eta_and_chi <approx> (model, nextpoint, line, next_freq, eta_n, chinext);

        Snext = eta_n / chinext;

        dtaunext = half * (chicurr + chinext) * dZ;
    }
}

///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template<ApproximationType approx>
accel inline void Solver :: solve_feautrier_order_2 (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_ ();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

            C[first] = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = 1.0 / C[first];   // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;

    if (n_off_diag == 0)
    {
        if (centre < last)
        {
            /// Write economically: G[last] = (B[last] - A[last]) / A[last];
            GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
            GP[last] = GG[last] / (one + GG[last]);

            for (long n = last-1; n > centre; n--) // use long in reverse loops!
            {
                Su[n] += Su[n+1] * FI[n];

                GG[n] = (C[n] * GP[n+1] + one) * inverse_A[n];
                GP[n] = GG[n] / (one + GG[n]);
            }

            Su    [centre] += Su[centre+1] * FI[centre];
            L_diag[centre]  = inverse_C[centre] / (FF[centre] + GP[centre+1]);
        }
        else
        {
            L_diag[centre] = (one + FF[centre-1]) / (Bl_min_Al + Bl*FF[centre-1]);
        }
    }
    else
    {
        /// Write economically: G[last] = (B[last] - A[last]) / A[last];
        GG[last] = half * Bl_min_Al * dtau_n * dtau_n;
        GI[last] = one / (one + GG[last]);
        GP[last] = GG[last] * GI[last];

        L_diag[last] = (one + FF[last-1]) / (Bl_min_Al + Bl*FF[last-1]);

        for (long n = last-1; n > first; n--) // use long in reverse loops!
        {
            Su[n] += Su[n+1] * FI[n];

            GG[n] = (C[n] * GP[n+1] + one) * inverse_A[n];
            GI[n] = one / (one + GG[n]);
            GP[n] = GG[n] * GI[n];

            L_diag[n] = inverse_C[n] / (FF[n] + GP[n+1]);
        }

        Su    [first] += Su[first+1] * FI[first];
        L_diag[first]  = (one + GG[first+1]) / (Bf_min_Cf + Bf*GG[first+1]);

        for (long n = last-1; n >= first; n--) // use long in reverse loops!
        {
            L_upper(0,n+1) = L_diag[n+1] * FI[n  ];
            L_lower(0,n  ) = L_diag[n  ] * GI[n+1];
        }

        for (Size m = 1; (m < n_off_diag) && (m < n_tot-1); m++)
        {
            for (long n = last-1-m; n >= first; n--) // use long in reverse loops!
            {
                L_upper(m,n+m+1) = L_upper(m-1,n+m+1) * FI[n    ];
                L_lower(m,n    ) = L_lower(m-1,n    ) * GI[n+m+1];
            }
        }
    }
}


///   Computes the optical depth assuming only a single line exists.
///    @param[in] curridx: index of the current point
///    @param[in] nextidx: index of the next point
///    @param[in] lineidx: index of the line for which to compute the optical depth
///    @param[in] curr_freq: current frequency (in comoving frame)
///    @param[in] next_freq: next frequency (in comoving frame)
///    @param[in] dz: distance increment
inline Real Solver :: compute_dtau_single_line(Model& model, Size curridx, Size nextidx, Size lineidx, Real curr_freq, Real next_freq, Real dz)
{
    const Real linefreq=model.lines.line[lineidx];
    const Real average_inverse_line_width=(model.lines.inverse_width(curridx, lineidx)+model.lines.inverse_width(nextidx, lineidx))/2.0;

    //opacity is stored divided by the linefreq, so multiply by it
    const Real curr_line_opacity=linefreq*model.lines.opacity(curridx, lineidx);
    const Real next_line_opacity=linefreq*model.lines.opacity(nextidx, lineidx);

    //if frequencies are equal, division by zero (due to the optical depth formula) happens if we were not to use this branch
    if (curr_freq==next_freq)
    {
        //doing the default computation instead (no shifting)
        const Real diff = curr_freq - model.lines.line[lineidx];//curr_freq==next_freq, so choice is arbitrary
        const Real prof = gaussian(average_inverse_line_width, diff);
        const Real average_opacity = (curr_line_opacity+next_line_opacity)/2.0;

        return dz * (prof * average_opacity + model.parameters->min_opacity);
    }

    //We assume a linear interpolation of these dimensionless frequency positions
    //We will also assume the line width to be somewhat constant, replacing the values with the averages
    const Real next_pos=(linefreq-next_freq)*average_inverse_line_width;
    const Real curr_pos=(linefreq-curr_freq)*average_inverse_line_width;

    //In this way, the diff_pos can be computed quite simple, and we do not have a discrepancy between the interpolation and the bounds
    const Real diff_pos=next_pos-curr_pos;

    /// the more correct approach, taking into account also the line opacity change; however, it does not make too much of a difference in the actual result and is quite a bit slower
    // const Real delta_opacity=(next_line_opacity-curr_line_opacity);
    // const Real deltanu=-next_freq+curr_freq;//differences in curr freqs; +-1 due to shift being defined in the other direction
    //
    // //note: opacity can also be extrapolated; however the correction term (expterm) accounts for that
    // const Real interp_opacity=curr_line_opacity+delta_opacity*(curr_freq-linefreq)/deltanu;
    //
    // //This term is a constant term, giving the usual ... as if the opacity were static
    // const Real erfterm=interp_opacity/diff_pos/2.0*(std::erf(next_pos)-std::erf(curr_pos));
    // //This term corrects for the fact that the opacity between points changes
    // const Real expterm=delta_opacity/2.0*INVERSE_SQRT_PI/diff_pos/diff_pos*(std::exp(-curr_pos*curr_pos)-std::exp(-next_pos*next_pos));
    // return dz*std::max(average_inverse_line_width*(erfterm+expterm), model.parameters->min_opacity);

    //If we instead use an average opacity, the computation is quite a bit faster
    const Real average_opacity=(next_line_opacity+curr_line_opacity)/2.0;
    const Real erfterm=average_opacity/diff_pos/2.0*(std::erff(next_pos)-std::erff(curr_pos));
    //correcting to bound opacity from below to the minimum opacity (assumes positive opacities occuring in the model)
    return dz*(average_inverse_line_width*erfterm+model.parameters->min_opacity);

}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template<ApproximationType approx>
accel inline void Solver :: image_feautrier_order_2 (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_ ();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);


    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;

    // for (long n = last-1; n > first; n--) // use long in reverse loops!
    // {
    //     Su[n] += Su[n+1] * FI[n];
    // }

    // Su[first] += Su[first+1] * FI[first];
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///    @param[in] w : width index
///////////////////////////////////////////////////////////////////////
template<ApproximationType approx>
accel inline void Solver :: image_feautrier_order_2_for_point_loc (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();
    Vector<Real>& GG = GG_();
    Vector<Real>& GI = GI_();
    Vector<Real>& GP = GP_();

    Vector<Real>& L_diag  = L_diag_ ();
    Matrix<Real>& L_upper = L_upper_();
    Matrix<Real>& L_lower = L_lower_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    //err, source function might be slightly different when looking at it from curr and next point
    // this is due to the weighting by the line optical depths; this might not be saved TODO think whether this is correct
    model.   S_ray(0, f) =  term_c;
    model.   S_ray(1, f) =  term_n;

    model.dtau_ray(0, f) = 0.0;
    model.dtau_ray(1, f) = dtau_n;


    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

    C[first] = two * inverse_dtau_f * inverse_dtau_f;

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Su[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);

    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        model. S_ray(n+1-first, f) =  term_n;
        model.dtau_ray(n+1-first, f) = dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;

    model.u_ray(last-first, f) = Su[last];


    for (long n = last-1; n > first; n--) // use long in reverse loops!
    {
        Su[n] += Su[n+1] * FI[n];

        model.u_ray(n-first, f) = Su[n];
    }

    Su[first] += Su[first+1] * FI[first];

    model.u_ray(0, f) = Su[first];
}

/// Image the optical depth
///////////////////////////////////////////////////////////////////////
template<ApproximationType approx>
accel inline void Solver :: image_optical_depth (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Real eta_c, chi_c;
    Real eta_n, chi_n;

    Real tau = 0.0;


    // Get optical properties for first two elements
    get_eta_and_chi <approx> (model, nr[first  ], l, freq*shift[first  ], eta_c, chi_c);
    get_eta_and_chi <approx> (model, nr[first+1], l, freq*shift[first+1], eta_n, chi_n);

    tau += half * (chi_c + chi_n) * dZ[first];


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        eta_c = eta_n;
        chi_c = chi_n;

        // Get new radiative properties
        get_eta_and_chi <approx> (model, nr[n+1], l, freq*shift[n+1], eta_n, chi_n);

        tau += half * (chi_c + chi_n) * dZ[n];
    }

    optical_depth_() = tau;
}


///  Solver for Feautrier equation along ray pairs using the (ordinary)
///  2nd-order solver, without adaptive optical depth increments
///////////////////////////////////////////////////////////////////////
template<ApproximationType approx>
accel inline void Solver :: solve_feautrier_order_2_uv (Model& model, const Size o, const Size f)
{
    const Real freq = model.radiation.frequencies.nu(o, f);
    const Size l    = model.radiation.frequencies.corresponding_line[f];

    Real eta_c, chi_c, dtau_c, term_c;
    Real eta_n, chi_n, dtau_n, term_n;

    const Size first = first_();
    const Size last  = last_ ();
    const Size n_tot = n_tot_();

    Vector<double>& dZ    = dZ_   ();
    Vector<Size  >& nr    = nr_   ();
    Vector<double>& shift = shift_();

    Vector<Real>& inverse_chi = inverse_chi_();

    Vector<Real>& Su = Su_();
    Vector<Real>& Sv = Sv_();

    Vector<Real>& A         = A_        ();
    Vector<Real>& C         = C_        ();
    Vector<Real>& inverse_A = inverse_A_();
    Vector<Real>& inverse_C = inverse_C_();

    Vector<Real>& FF = FF_();
    Vector<Real>& FI = FI_();

    bool compute_curr_opacity = true; // for the first point, we need to compute both the curr and next opacity (and source)

    compute_source_dtau<approx>(model, nr[first], nr[first+1], l, freq*shift[first], freq*shift[first+1], shift[first], shift[first+1], dZ[first], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

    // Set boundary conditions
    const Real inverse_dtau_f = one / dtau_n;

            C[first] = two * inverse_dtau_f * inverse_dtau_f;
    inverse_C[first] = one / C[first];   // Required for Lambda_diag

    const Real Bf_min_Cf = one + two * inverse_dtau_f;
    const Real Bf        = Bf_min_Cf + C[first];
    const Real I_bdy_f   = boundary_intensity (model, nr[first], freq*shift[first]);

    Su[first]  = term_c + two * I_bdy_f * inverse_dtau_f;
    Sv[first]  = two * inverse_dtau_f * (I_bdy_f - term_c);

    Su[first] /= Bf;
    Sv[first] /= Bf;

    /// Write economically: F[first] = (B[first] - C[first]) / C[first];
    FF[first] = half * Bf_min_Cf * dtau_n * dtau_n;
    FI[first] = one / (one + FF[first]);


    /// Set body of Feautrier matrix
    for (Size n = first+1; n < last; n++)
    {
        term_c = term_n;
        dtau_c = dtau_n;
         eta_c =  eta_n;
         chi_c =  chi_n;

        compute_source_dtau<approx>(model, nr[n], nr[n+1], l, freq*shift[n], freq*shift[n+1], shift[n], shift[n+1], dZ[n], compute_curr_opacity, dtau_n, chi_c, chi_n, term_c, term_n);

        const Real dtau_avg = half * (dtau_c + dtau_n);
        inverse_A[n] = dtau_avg * dtau_c;
        inverse_C[n] = dtau_avg * dtau_n;

        A[n] = one / inverse_A[n];
        C[n] = one / inverse_C[n];

        /// Use the previously stored value of the source function
        Su[n] = term_c;

        FF[n] = (A[n] * FF[n-1] * FI[n-1] + one) * inverse_C[n];
        FI[n] = one / (one + FF[n]);
        Su[n] = (A[n] * Su[n-1] + Su[n]) * FI[n] * inverse_C[n];
        Sv[n] = (A[n] * Sv[n-1]        ) * FI[n] * inverse_C[n];
    }


    /// Set boundary conditions
    const Real inverse_dtau_l = one / dtau_n;

    A[last] = two * inverse_dtau_l * inverse_dtau_l;

    const Real Bl_min_Al = one + two * inverse_dtau_l;
    const Real Bl        = Bl_min_Al + A[last];

    const Real denominator = one / (Bl * FF[last-1] + Bl_min_Al);

    const Real I_bdy_l = boundary_intensity (model, nr[last], freq*shift[last]);

    Su[last] = term_n + two * I_bdy_l * inverse_dtau_l;
    Sv[last] = two * inverse_dtau_l * (I_bdy_l - term_n);

    Su[last] = (A[last] * Su[last-1] + Su[last]) * (one + FF[last-1]) * denominator;
    Sv[last] = (A[last] * Sv[last-1]           ) * (one + FF[last-1]) * denominator;

    if (centre < last)
    {
        for (long n = last-1; n >= centre; n--) // use long in reverse loops!
        {
            Su[n] += Su[n+1] * FI[n];
            Sv[n] += Sv[n+1] * FI[n];
        }
    }

}


accel inline void Solver :: set_eta_and_chi (Model& model, const Size rr) const
{
    model.eta.resize (model.parameters->npoints(), model.parameters->nfreqs());
    model.chi.resize (model.parameters->npoints(), model.parameters->nfreqs());

    for (Size p = 0; p < model.parameters->npoints(); p++)
    {
        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            // Extract the Doppler shift
            const double shift = model.geometry.get_shift <Rest> (0, rr, p, 0.0);
            const Real   freq  = model.radiation.frequencies.nu(0, f);
            const Size   l     = model.radiation.frequencies.corresponding_line[f];

            get_eta_and_chi <None> (model, p, l, freq*shift, model.eta(p,f), model.chi(p,f));
        }
    }
}


accel inline void Solver :: set_boundary_condition (Model& model) const
{
    model.boundary_condition.resize (model.parameters->nboundary(), model.parameters->nfreqs());

    for (Size b = 0; b < model.parameters->nboundary(); b++)
    {
        const Size p = model.geometry.boundary.boundary2point[b];

        for (Size f = 0; f < model.parameters->nfreqs(); f++)
        {
            const Real freq = model.radiation.frequencies.nu(0, f);

            model.boundary_condition(b,f) = boundary_intensity (model, p, freq);
        }
    }
}


inline void Solver :: set_column (Model& model) const
{
    model.column.resize (model.parameters->nrays(), model.parameters->npoints());

    for (Size rr = 0; rr < model.parameters->hnrays(); rr++)
    {
        const Size ar = model.geometry.rays.antipod[rr];

        cout << "--- rr = " << rr << endl;

        accelerated_for (o, model.parameters->npoints(),
        {
            model.column(rr, o) = get_column(model, o, rr);
            model.column(ar, o) = get_column(model, o, ar);
        })
    }

}


accel inline Real Solver :: get_column (const Model& model, const Size o, const Size r) const
{
    Real column = 0.0;

    double  Z = 0.0;   // distance from origin (o)
    double dZ = 0.0;   // last increment in Z

    Size nxt = model.geometry.get_next (o, r, o, Z, dZ);

    if (model.geometry.valid_point(nxt))
    {
        Size crt = o;

        column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;

        while (model.geometry.not_on_boundary(nxt))
        {
            crt = nxt;
            nxt = model.geometry.get_next (o, r, nxt, Z, dZ);

            column += 0.5 * (model.density[crt] + model.density[nxt]) * dZ;
        }
    }

    return column;
}
