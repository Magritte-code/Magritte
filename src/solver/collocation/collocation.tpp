#include <cmath>
#include <Eigen/Core>
#include <map>
#include <limits>

//Using a constructed kd tree, directly return the distance towards to farthest of the 'nb_neighbors_to_query' neighbors.
//Warning: the number of neighbors to query must be less or equal than the number of points in the kdtree
inline Real get_kd_tree_basis_radius(Size basis_point, kd_tree& kdtree, Size nb_neighbors_to_query, Geometry& geometry)
{
    std::vector<double> query_pt(3);
    query_pt[0]=geometry.points.position[basis_point].x();
    query_pt[1]=geometry.points.position[basis_point].y();
    query_pt[2]=geometry.points.position[basis_point].z();

    // Note to self, nanoflann requires size_t (unsigned int), so that is why we do not use our Size (long unsinged int) here
    vector<size_t> ret_indexes(nb_neighbors_to_query);//still need to be translated to the actual indices
    vector<double> out_dists_sqr(nb_neighbors_to_query);//note: we are not really using this, but nanoflann requires it

    nanoflann::KNNResultSet<double> resultSet(nb_neighbors_to_query);

    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
    kdtree.index->findNeighbors(resultSet, &query_pt[0],
                                nanoflann::SearchParams(10));

    return (Real) std::sqrt(*max_element(out_dists_sqr.begin(), out_dists_sqr.end()));
}


//Returns all the indices of points within a radius around a given point
inline vector<Size> get_kd_tree_basis_elements_in_radius(Size basis_point, kd_tree& kdtree, Real radius, vector<Size>& index_conversion, Geometry& geometry)
{
    std::vector<double> query_pt(3);
    query_pt[0]=geometry.points.position[basis_point].x();
    query_pt[1]=geometry.points.position[basis_point].y();
    query_pt[2]=geometry.points.position[basis_point].z();

    // Note to self, nanoflann requires size_t (unsigned int), so that is why we do not use our Size (long unsinged int) here
    std::vector<std::pair<long int, double>> ret_matches;//returns index and distance
    nanoflann::SearchParams params;

    // std::cout<<"radius: "<<radius<<std::endl;
    // std::cout<<"radius squared:"<<std::pow(radius,2)<<std::endl;

    //Apparently, nanoflann requires the square of the radius for this function; this is nowhere found in the docs...
    const size_t n_neighbors = kdtree.index->radiusSearch(
          &query_pt[0], std::pow(radius,2), ret_matches, params);

    std::cout<<"ret_matches size: "<<ret_matches.size()<<std::endl;

    // vector<Size> to_return;
    // to_return.resize(n_neighbors);
    // for (size_t i=0; i<n_neighbors; i++)
    // {
    //     to_return[i]=ret_matches[i].first;
    //     // std::cout<<"for point "<<p<<" found neighbor: "<<true_indices[i]<<std::endl;
    // }
    // return to_return

    vector<Size> true_indices;
    true_indices.resize(n_neighbors);
    for (size_t i=0; i<n_neighbors; i++)
    {
        true_indices[i]=index_conversion[ret_matches[i].first];
        // std::cout<<"for point "<<p<<" found neighbor: "<<true_indices[i]<<std::endl;
    }
    return true_indices;

}


/// For now, we will create the kd tree on the fly
/// This is due to default initialization (and anything which remotely looks like it) being impossible for these trees, i guess the data just goes out of scope
//std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 3>, Size1>
inline Eigen::Matrix<double, Eigen::Dynamic, 3> create_mat_for_points_idx(Model& model, vector<Size> &indices_vector, Geometry& geometry)
{
    const Size total_points=indices_vector.size();
    Eigen::Matrix<double, Eigen::Dynamic, 3> temp_mat(total_points, 3);

    Size i=0;//just a temporary index
    for (auto point_to_add:indices_vector)
    {//add each point to the vector of vectors
        // add x, y and z
        temp_mat(i,0)=model.geometry.points.position[point_to_add].x();
        temp_mat(i,1)=model.geometry.points.position[point_to_add].y();
        temp_mat(i,2)=model.geometry.points.position[point_to_add].z();
        i=i+1;
    }

    // For creating the kd_tree, add
    //Set the kd_tree
    // kd_tree mat_index(3, create_mat_for_kd_tree_of_lvl(Size lvl), 10 /* max leaf */);
    // mat_index.index->buildIndex();

    // return std::make_tuple(temp_mat,indices_vector);//index conversion is defined as the second input parameter, so not necessary to return here
    return temp_mat;
}


// Assumes the point indices to be ordered in remaining points; assumes the boundary points to lie first in point indices
// Also assumes the frequency quadrature to be already defined
//ASSUMES: the opacity does not actually depend on the ray direction (or at least not significantly)
inline void Collocation :: prune_positional_bases(Model& model)
{
    // vector<vector<Size>> remaining_points;// temporary
    remaining_points_after_pruning.resize(n_freq_bases);
    radial_basis_functions_having_evaluation_at.resize(n_freq_bases);
    n_point_bases.resize(n_freq_bases);
    cum_n_point_bases.resize(n_freq_bases+1);
    cum_n_point_bases[0]=0;

    //compute kd tree with all points
    // points_in_grid=model.geometry.points.multiscale.get_current_points_in_grid();//MOVED
    // n_points_in_grid=points_in_grid.size();
    auto temp_mat=create_mat_for_points_idx(model, points_in_grid, model.geometry);
    // temp_mat=std::get<0>(tuple);
    // index_conversion=std::get<1>(tuple);
    kd_tree full_mat_index(3, temp_mat, 10 /* max leaf */);
    full_mat_index.index->buildIndex();

    std::map<Size,Size> indexmap;//for simplicity, maps the indices from points_in_grid to consecutive indices
    Size temp_idx=0;
    for(Size point: points_in_grid)
    {
        indexmap.insert(std::pair<Size,Size>(point,temp_idx));//maps from the point in the original grid to the point in the new grid
        temp_idx++;
    }


    //start with assuming all positional bases remain at all frequencies, so we can prune them
    for (Size freqid=0; freqid<n_freq_bases; freqid++)
    {
        // remaining_points[freqid].resize(parameters.npoints());//start with all points
        remaining_points_after_pruning[freqid]=model.geometry.points.multiscale.get_current_points_in_grid();
    }

    rbf_radii.resize(n_freq_bases);
    basis_point_idx.resize(n_freq_bases);//maps to the original grid points
    local_basis_point_idx.resize(n_freq_bases);//maps to the collocation grid points

    //for all frequencies, we can delete basis functions such that a typical optical depth is not negligible
    for (Size freqid=0; freqid<n_freq_bases; freqid++)
    {
        radial_basis_functions_having_evaluation_at[freqid].resize(n_points_in_grid);

        if (PRUNING_POSITIONAL_BASIS_FUNCTIONS)
        {
            // the opacities are computed COMOVING; as pruning is done in order to increase distance between our doppler shifted basis functions
            //when defining static frequency basis functions, this is also valid (as the frequency not longer depends on the location); either way, this just influences some pruning condition
            vector<Real> comoving_opacities;
            comoving_opacities.resize(n_points_in_grid);//should be computed over all points in the current grid
            for (Size pointid=0; pointid<n_points_in_grid; pointid++)
            {
                //ASSUMES: the opacity does not actually depend on the ray direction (or at least not significantly)
                comoving_opacities[pointid]=get_opacity(model, 0, freqid, pointid);
            }
            //for loop get_opacity(model, rayidx, freqidx, pointidx)

            std::set<Size> points_to_prune;//for this specific frequency
            //for all bases not yet deleted, check if neighbors can be deleted
            //for this, we merely iterate over all points, check whether the bases are already deleted, and if not,
            // search around the point in a radius defined by the pruning condition
            for (Size center_point_id=0; center_point_id<n_points_in_grid; center_point_id++)
            {
                const Size center_point=points_in_grid[center_point_id];
                const Real center_opacity=comoving_opacities[center_point_id];
                //if not already pruned, the basis at the point can prune others
                if (points_to_prune.find(center_point)!=points_to_prune.end())
                {
                    //do nothing if already pruned
                    continue;
                }
                const Real search_radius=MIN_OPTICAL_DEPTH_DIFF/center_opacity;
                vector<Size> points_within_radius=get_kd_tree_basis_elements_in_radius(center_point, full_mat_index, search_radius, points_in_grid, model.geometry);//TODO: remove index conversion
                std::cout<<"points_within_radius size: "<<points_within_radius.size()<<std::endl;
                for (Size neighbor: points_within_radius)
                {
                    if (neighbor!=center_point && model.geometry.not_on_boundary(neighbor))//precaution against a basis point deciding to delete itself
                    {// also precaution against deleting boundary points (as the boundary condition equation might also get deleted)
                        //check if distance*opacity<threshold; if so prune the basis at the neighboring location (and remove them from the remaining_points list not yet already done)
                        Vector3D temp_vector=model.geometry.points.position[neighbor]-model.geometry.points.position[center_point];
                        const Real temp_dist=std::sqrt(temp_vector.squaredNorm());
                        const Real other_opacity=get_opacity_bound(model, indexmap.at(neighbor), freqid, center_point_id);
                        // const Real other_opacity=comoving_opacities[indexmap.at(neighbor)];//DO NOT USE THE COMOVING OPACITY FOR THE OTHER, as it results in far too few basis functions for the line wings
                        const Real max_opacity=std::max(center_opacity, other_opacity);//the main disadvantage of this approach, is that the path inbetween is not included...
                        if (max_opacity*temp_dist<MIN_OPTICAL_DEPTH_DIFF)//bound for optical depth
                        {
                            //TODO: delete other basis function; but do check whether it makes sense to readd in least squares sense
                            points_to_prune.insert(neighbor);
                        }
                    }
                }
            }

        // bool deleted_basis=true;//whether we deleted a basis in the previous iteration
        // while (deleted_basis)//if we deleted a basis, we might still need to delete more bases
        // {
        //     deleted_basis=false;
        //
        //     //compute kd tree with current points
        //     auto temp_mat=create_mat_for_points_idx(model, remaining_points_after_pruning[freqid]);
        //     // temp_mat=std::get<0>(tuple);
        //     // index_conversion=std::get<1>(tuple);
        //     kd_tree mat_index(3, temp_mat, 10 /* max leaf */);
        //     mat_index.index->buildIndex();
        //
        //     // vector<Size> points_to_prune;//for the specific frequency
        //     std::set<Size> points_to_prune;//for this specific frequency
        //     // points_to_prune.resize(remaining_points_after_pruning[freqid].size());//upper bound
        //     for (Size pointid: remaining_points_after_pruning[freq_id])
        //     {
        //         //if not already pruned, the basis at the point can prune others
        //         if (points_to_prune.find(pointid)!=points_to_prune.end())
        //         {
        //             //do nothing
        //             continue;
        //         }
        //
        //         //TODO: get neighbors; use kdtree
        //         std::vector<Size> neighbors=model.get_coarser_neighbors_kd_tree(pointid, mat_index, remaining_points_after_pruning[freqid]);//also includes itself
        //         const Real curr_opacity=0;//TODO: add comoving opacity
        //         for (Size neighbor: neighbors)
        //         {
        //             if (neighbor!=pointidx)//precaution against a basis deciding to delete itself
        //             {
        //                 //check if distance*opacity<threshold; if so prune the basis at the neighboring location (and remove them from the remaining_points list not yet already done)
        //                 Vector3D temp_vector=model.geometry.points.position[neighbor]-model.geometry.points.position[pointid];
        //                 Real temp_dist=std::sqrt(temp_vector.squaredNorm());
        //                 Real other_opacity=0;//TODO: precompute all opacities at the 'current' comoving frequency
        //                 Real max_opacity=std::max(curr_opacity, other_opacity);
        //                 if (max_opacity*temp_dist<CONSTANT)//bound for optical depth
        //                 {
        //                     //TODO: delete other basis function; but do check whether it makes sense to readd in least squares sense
        //                     points_to_prune.insert(neighbor);
        //                     deleted_basis=true;
        //                 }
        //             }
        //         }
        //     }

            //and finally prune the basis functions associated to the points
            vector<Size>::iterator it = remaining_points_after_pruning[freqid].begin();

            while(it != remaining_points_after_pruning[freqid].end())
            {
                if(points_to_prune.find(*it)!=points_to_prune.end())
                {
                    it = remaining_points_after_pruning[freqid].erase(it);
                }
                else {it++;}
            }

        }

        //now, we have our remaining center basis points, so we can define the typical radii (using only these remaining center points)
        //compute kd tree with current points
        auto temp_mat=create_mat_for_points_idx(model, remaining_points_after_pruning[freqid], model.geometry);
        kd_tree mat_index(3, temp_mat, 10 /* max leaf */);
        mat_index.index->buildIndex();
        //for each remaining basis point, figure out the typical radius
        Size temp_counter=0;
        const Size nb_positional_basis_funs=remaining_points_after_pruning[freqid].size();//for this frequency, the number of positional basis functions left after pruning
        // const Size nb_neighbors_to_query=N_POINTS_IN_RADIAL_BASIS;//we cannot query the radius if we use more points than there remain in the grid
        const Size nb_neighbors_to_query=std::min(N_POINTS_IN_RADIAL_BASIS,nb_positional_basis_funs);//we cannot query the radius if we use more points than there remain in the grid
        //FIXME: make sure that the rsulting radius is larger than 0
        n_point_bases[freqid]=nb_positional_basis_funs;
        cum_n_point_bases[freqid+1]=cum_n_point_bases[freqid]+nb_positional_basis_funs;
        basis_point_idx[freqid].resize(nb_positional_basis_funs);//maps to the original grid points
        local_basis_point_idx[freqid].resize(nb_positional_basis_funs);//maps to the collocation grid points
        rbf_radii[freqid].resize(nb_positional_basis_funs);
        for (Size basis_center: remaining_points_after_pruning[freqid])
        {
            basis_point_idx[freqid][temp_counter]=basis_center;
            local_basis_point_idx[freqid][temp_counter]=indexmap.at(basis_center);
            rbf_radii[freqid][temp_counter]=get_kd_tree_basis_radius(basis_center, mat_index, nb_neighbors_to_query, model.geometry);
            std::cout<<"rbf radius: "<<rbf_radii[freqid][temp_counter]<<std::endl;
            temp_counter++;
        }
        std::cout<<"nb_positional_basis_funs: "<<nb_positional_basis_funs<<std::endl;
        //now that we have the typical radius of all basis functions, we can figure out what points lie in each basis
        for (Size basis_point_id=0; basis_point_id<nb_positional_basis_funs; basis_point_id++)
        {
            const Size center_point=basis_point_idx[freqid][basis_point_id];
            vector<Size> points_within_radius=get_kd_tree_basis_elements_in_radius(center_point, full_mat_index, rbf_radii[freqid][basis_point_id], points_in_grid, model.geometry);
            //TODO: instead of mapping it from and to the original grid, wouldn't it be better to just use something like iota... as index_conversion parameter?
            std::cout<<"n_point_within_radius: "<<points_within_radius.size()<<std::endl;
            for (Size point_in_range: points_within_radius)
            {
                radial_basis_functions_having_evaluation_at[freqid][indexmap.at(point_in_range)].insert(basis_point_id);
            }
        }
    }
}



// Currently only really sets the interacting frequency bases
inline void Collocation :: set_interacting_bases(Model& model)
{
  // REPLACED by prune_positional_bases
  // //Set the kd_tree
  // Eigen::Matrix<double, Eigen::Dynamic, 3> temp_mat;
  // // Size1 index_conversion;
  // auto tuple=model.create_mat_for_kd_tree_of_lvl(0);
  // temp_mat=std::get<0>(tuple);
  // index_conversion=std::get<1>(tuple);
  // // The reverse map of index_conversion
  // // reverse_index_conversion
  // Size temp_counter=0;
  // for (Size index: index_conversion)
  // {
  //     reverse_index_conversion.insert(std::pair<Size,Size>(index,temp_counter));
  //     temp_counter++;
  // }
  //
  // kd_tree mat_index(3, temp_mat, 10 /* max leaf */);
  // mat_index.index->buildIndex();
  //
  //   for (Size pointid=0; pointid<parameters.npoints(); pointid++)
  //   {
  //       vector<Size> neighbors_coarser_grid=model.get_coarser_neighbors_kd_tree(pointid, mat_index, index_conversion);
  //       Real max_dist=0.0;
  //       Size farthest_point_idx=parameters.npoints();
  //       // Size furthest_point_idx=parameters.npoints();
  //       // TODO: prune furthest point from this
  //       for (Size coarser_neighbor: neighbors_coarser_grid)
  //       {
  //           Vector3D temp_vector=(model.geometry.points.position[coarser_neighbor]-model.geometry.points.position[pointid]);
  //           Real temp_dist=std::sqrt(temp_vector.squaredNorm());
  //           // if (temp_dist>max_dist)
  //           // {
  //           max_dist=temp_dist;
  //           farthest_point_idx=coarser_neighbor;
  //           // }
  //       }
  //       // Vector3D temp_vector=(model.geometry.points.position[neighbors_coarser_grid.back()]-model.geometry.points.position[pointid]);
  //       rbf_radii[pointid]=max_dist;
  //       std::cout<<"radius: "<<rbf_radii[pointid]<<std::endl;
  //       for (Size coarser_neighbor: neighbors_coarser_grid)
  //       {
  //           if (coarser_neighbor!=farthest_point_idx)
  //           {//set for every neighbors that they have a non-zero value in the rbf of pointid
  //             //TODO: needs frequency dependence; and for least squares, it also needs to act on position
  //             //so rename to: radial_basis_function_having_nonzero_evaluation_at[freq][point] ....
  //           other_radial_basis_interacting_with[coarser_neighbor].insert(pointid);//this also includes itself
  //           }
  //       }
  //
  //   }

    // for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    for (Size pointid=0; pointid<n_points_in_grid; pointid++)
    {
        Size point_in_grid=points_in_grid[pointid];
        other_frequency_basis_interacting_with[pointid].resize(n_freq_bases);

        //currently, i am choosing the radius such that the last point lies on the boundary of the domain of the radial basis function
        //Therefore, we do not need to include it in the vector of interacting basis functions
        // vector<Size> subvector={neighbors_coarser_grid.begin(),neighbors_coarser_grid.end()-1};


        // vector<Real> min_doppler_shift(parameters.nrays(),1.0);
        Vector3D velocity_curr_point=model.geometry.points.velocity[point_in_grid];
        // for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
        // {
        // for (Size neighbor: subvector)
        //TODO: possibly change this for the comoving frame...
        // for (Size neighbor: other_radial_basis_interacting_with[pointid])
        // {
        //     max_doppler_shift=std::max(max_doppler_shift,Real(
        // 1+std::sqrt((model.geometry.points.velocity[neighbor]-velocity_curr_point).squaredNorm())));
        //     // min_doppler_shift[rayidx]=std::min(max_doppler_shift,
        // // 1+std::sqrt((model.geometry.points.velocity[neighbor]-velocity_curr_point).dot(model.geometry.rays.direction[rayidx])/CC);
        // }


        // For every frequency, calculate which frequencies lie nearby
        //now using all the quadrature frequencies
        // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
        //TODO: use more efficient estimates?
        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {
            //calculate the max doppler shift possible //TODO: should be compared to the basis functions which have an evaluation in the point
            Real max_doppler_shift=1.0;
            for (Size rbf_center: radial_basis_functions_having_evaluation_at[freqid][pointid])
            {
                max_doppler_shift=std::max(max_doppler_shift,Real(
                  1+std::sqrt((model.geometry.points.velocity[basis_point_idx[freqid][rbf_center]]-velocity_curr_point).squaredNorm())));
                  // min_doppler_shift[rayidx]=std::min(max_doppler_shift,
                  // 1+std::sqrt((model.geometry.points.velocity[neighbor]-velocity_curr_point).dot(model.geometry.rays.direction[rayidx])/CC);
            }
            // const Size corresp_lineid=?; for more efficient estimate, just try whether the line centers lie close enough
            //The line frequencies are not ordered, so just search the whole list...
            // for (Size other_lineid=0; other_lineid<parameters.nlines(); other_lineid++)
            for (Size other_freqid=0; other_freqid<n_freq_bases; other_freqid++)
            {
                // Size l=model.radiation.frequencies.corresponding_l_for_spec[other_freqid];   ///< number of line species corresponding to frequency
                // Size k=model.radiation.frequencies.corresponding_k_for_tran[other_freqid];   ///< number of transition corresponding to frequency
                // Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];   ///< number of line number corresponding to frequency
                // // std::cout<<"l, k, z: "<<l<<", "<<k<<", "<<z<<std::endl;
                // Size corresp_other_lineid=model.lines.line_index(l, k);
                // // Somewhat loose uniform (over direction and neighbors) estimate
                // // Start by bounding the apparent frequency nu' for the point (pointidx, lineid): nu/max_dpplr_shift<nu'<nu*max_dpplr_shift
                // // if ((model.lines.line[lineid]/max_doppler_shift<model.lines.line[other_lineid]+1/model.lines.inverse_width(pointid,other_lineid)*TRUNCATION_SIGMA*max_doppler_shift)
                // // &&  (model.lines.line[lineid]*max_doppler_shift>model.lines.line[other_lineid]-1/model.lines.inverse_width(pointid,other_lineid)*TRUNCATION_SIGMA*max_doppler_shift))
                // //TODO: other possibility is to use non_doppler_shifted_inverse_widths instead of lines.inverse_width
                // //FIXME: use the frequency quadrature
                // // if ((model.radiation.frequencies.nu(pointid,freqid)/max_doppler_shift<model.radiation.frequencies.nu(pointid,other_freqid)+1/model.lines.inverse_width(pointid,corresp_other_lineid)*TRUNCATION_SIGMA*max_doppler_shift)
                // // &&  (model.radiation.frequencies.nu(pointid,freqid)*max_doppler_shift>model.radiation.frequencies.nu(pointid,other_freqid)-1/model.lines.inverse_width(pointid,corresp_other_lineid)*TRUNCATION_SIGMA*max_doppler_shift))
                // if (parameters.nquads()==1)//we have no other relevant inverse width than the inverse line width, so use the old formulation
                // {
                //     //in this case, freqid=lineid, so we can just fill in the freqid in the line id fields
                //     // if ((model.radiation.frequencies.nu(pointid,freqid)/max_doppler_shift<model.radiation.frequencies.nu(pointid,other_freqid)+1/model.lines.inverse_width(pointid,corresp_other_lineid)*TRUNCATION_SIGMA*max_doppler_shift)
                //     // &&  (model.radiation.frequencies.nu(pointid,freqid)*max_doppler_shift>model.radiation.frequencies.nu(pointid,other_freqid)-1/model.lines.inverse_width(pointid,corresp_other_lineid)*TRUNCATION_SIGMA*max_doppler_shift))
                //     if ((model.radiation.frequencies.nu(pointid,freqid)/max_doppler_shift<model.radiation.frequencies.nu(pointid,other_freqid)+1/model.lines.inverse_width(pointid,other_freqid)*TRUNCATION_SIGMA*max_doppler_shift)
                //     &&  (model.radiation.frequencies.nu(pointid,freqid)*max_doppler_shift>model.radiation.frequencies.nu(pointid,other_freqid)-1/model.lines.inverse_width(pointid,other_freqid)*TRUNCATION_SIGMA*max_doppler_shift))
                //     {
                //         other_frequency_basis_interacting_with[pointid][freqid].insert(other_freqid);
                //     }
                // }
                // else
                // {
                if ((frequency_quadrature[freqid]/max_doppler_shift<frequency_quadrature[other_freqid]+1/frequency_inverse_width_quadrature[other_freqid]*TRUNCATION_SIGMA*max_doppler_shift)
                &&  (frequency_quadrature[freqid]*max_doppler_shift>frequency_quadrature[other_freqid]-1/frequency_inverse_width_quadrature[other_freqid]*TRUNCATION_SIGMA*max_doppler_shift))
                {
                    other_frequency_basis_interacting_with[pointid][freqid].insert(other_freqid);
                }
                // }
            }
        }


        // for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        // {
        //     other_frequency_basis_interacting_with[rayid][pointid].resize(subvector.size());
        //     for (Size temp=0; temp<subvector.size(); temp++)
        //     {
        //         other_frequency_basis_interacting_with[rayid][pointid][temp].resize(parameters.nfreqs());
        //     }
        // }
        // also set what interactions there are between the frequency basis functions, which unfortunately also depend on the interacting points
        // Size reduced_index=0;
        // for (Size other_point: subvector)
        // {
        //     for (Size lineid=0; lineid<parameters.nfreqs(); lineid++)
        //     {
        //         for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        //         {
        //             Real curr_freq=frequencies[rayid][pointid][lineid];
        //             vector<Real> searchfreqs=frequencies[rayid][other_point];
        //             for (Size search_freq_idx=0; search_freq_idx<parameters.nfreqs(); search_freq_idx++)
        //             {   // The frequency basis function at the other point and search frequency is only nonzero for the frequency frequencies[pointid][lineid] if the main frequency of that gaussian does not lie to far from it (compared to the width)
        //                 Real abs_freq_diff=std::abs(curr_freq-searchfreqs[search_freq_idx]);
        //                 if (!freq_cutoff_condition(abs_freq_diff,frequencies_widths[rayid][other_point][search_freq_idx]))
        //                 {
        //                     other_frequency_basis_interacting_with[rayid][pointid][reduced_index][lineid].insert(search_freq_idx);
        //                 }
        //             }
        //         }
        //     }
        //     reduced_index++;
        // }
    }

}


/// Dev note: in geomeetry.tpp, similar functions are defined; however, these functions did not state exactly how the ray direction and the velocity difference were specified
///  @param[in] velocity_of_distant_object: the velocity of the distant object minus your local velocity
///  @param[in] ray_direction: the direction of the ray coming towards you
inline Real calculate_doppler_shift_observer_frame(Vector3D velocity_of_distant_object, Vector3D ray_direction)
{
    return 1+(velocity_of_distant_object.dot(ray_direction));//velocity is already divided by the light speed
}



// inline void Collocation :: setup_eval_matrix(Model& model)
// {
//     vector<Triplet<Real, Size>> eval_triplets; //triplets for the evaluation matrix
//     for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
//     {
//         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
//         {
//             for (Size pointidx=0; pointidx<point_locations.size(); pointidx++)
//             {
//                 //Now creating all the triplets, thus we need all nonzero basis triplets
//                 std::set<std::vector<Size>> nonzero_basis_triplets;
//                 get_nonzero_basis_triplets(nonzero_basis_triplets, rayidx, lineidx, pointidx);
//
//                 const Size curr_mat_idx=get_mat_index(rayidx, lineidx, pointidx);
//
//
//             }
//         }
//     }
// }

inline void Collocation :: setup_basis(Model& model)
{
    //Duplicate of all points in the current grid
    points_in_grid=model.geometry.points.multiscale.get_current_points_in_grid();
    n_points_in_grid=points_in_grid.size();
    //then we also set all the radii and locations for the radial basis functions
    //currently just uses all the points // might need to be changed later to all the points in a certain grid // note: currently also adding the boundary points
    point_locations.resize(n_points_in_grid);
    // rbf_radii.resize(parameters.npoints());
    // other_radial_basis_interacting_with.resize(n_points_in_grid);
    //We need all frequencies at ray and each point (can be different at each point due to doppler shift)
    frequencies.resize(parameters.nrays());
    frequencies_inverse_widths.resize(parameters.nrays());
    other_frequency_basis_interacting_with.resize(n_points_in_grid);  //?with a bit more computational effort, this could be reduced a bit in size
    // for (Size i:model.geometry.points.multiscale.get_current_points_in_grid())
    for (Size pointid=0; pointid<n_points_in_grid; pointid++)
    {
        point_locations[pointid]=model.geometry.points.position[points_in_grid[pointid]];
        // other_frequency_basis_interacting_with[pointid].resize(parameters.nlines());
        // other_frequency_basis_interacting_with[pointid].resize(parameters.nfreqs()); TODO PUT SOMEWHERE ELSE/DONE
    }

    // We need the inverse line widths, so we might as well make sure that they are calculated now//TODO: remove this here and just say that we need them before
    model.compute_inverse_line_widths();

    // frequency_quadrature.resize(parameters.nfreqs());
    // frequency_inverse_width_quadrature.resize(parameters.nfreqs());

    if (parameters.nquads()==1)//TODO: figure out whether this still has any place in the code, as we are completely ignoring the default quadratures parameter
    {//in this case, nfreqs=nlines
        frequency_quadrature.resize(parameters.nfreqs());
        frequency_inverse_width_quadrature.resize(parameters.nfreqs());
        for (Size lineid=0; lineid<parameters.nlines(); lineid++)
        {
            center_line_freq_indices.resize(parameters.nlines());
            frequency_quadrature[lineid]=model.lines.line[lineid];
            Real max_inverse_width=0;
            for (Size pointid: points_in_grid)
            {
                const Real curr_inverse_width=model.lines.inverse_width(pointid,lineid);
                max_inverse_width=std::max(max_inverse_width,curr_inverse_width);
            }
            frequency_inverse_width_quadrature[lineid]=max_inverse_width; //(just the maximum of some default for every point)
            center_line_freq_indices[lineid]=lineid;
        }
    }
    else
    {
        // //factors for determining the number of wing frequency quadrature points, independent of the exact line
        // const Real prefactor=(FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR-1.0)/FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR;
        // const Real div_factor=std::log(FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR);
        //
        // vector<Real> number_wing_quadrature_points_needed;//for every line; this might be different per line (depending on line width vs line frequency)
        // number_wing_quadrature_points_needed.resize(parameters.nlines());
        // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
        // {
        //     // So we start from the line center, compute the entire frequency range which needs sampling and then distribute the frequency quadrature (starting from the line center)
        //     // ALSO TODO: find easy formula to compute how many quadrature points we need in order to sample an entire range
        //     Real number_wing_quadrature_points_needed=0;//by definition this is symmetric (same for both sides of the line center)
        //     for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        //     {
        //         for (Size pointid: points_in_grid)
        //         {
        //             const Real line_freq=model.lines.line[lineid];
        //             const Real line_freq_inv_width=model.lines.inverse_width(pointid,lineid);
        //             //TODO: it might be faster to first compute all doppler shifts and then take min and max
        //             const Real temp_doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointid],model.geometry.rays.direction[rayid]);
        //             const Real doppler_width=line_freq*std::abs(temp_doppler_shift-1);//doppler shift applied to center line frequency
        //             std::cout<<"doppler width: "<<doppler_width<<std::endl;
        //
        //             const double local_number_wing_quadratures_bound=std::log(prefactor*(doppler_width/(FREQ_QUADRATURE_CENTER_SEPARATION_DISTANCE))+1.0)/div_factor;
        //             number_wing_quadrature_points_needed_in_line=std::max(number_wing_quadrature_points_needed_in_line,local_number_wing_quadratures_bound);
        //         }
        //     }
        //     //and round up
        //     number_wing_quadrature_points_needed[lineid]=std::ceil(number_wing_quadrature_points_needed_in_line);
        // }
        // Size total_nb_freq_quad_points=0;
        // for (Size lineid=0; lineid<parameters.nlines(); lineid++)//sum for every line how many frequency quadratures we need?
        // {
        //     total_nb_freq_quad_points+=2*number_wing_quadrature_points_needed[lineid]+FREQ_QUADRATURE_CENTER_N_ELEMENTS;
        // }

        // frequency_quadrature.resize(?total_nb_freq_quad_points);
        // frequency_inverse_width_quadrature.resize(?total_nb_freq_quad_points);

        //compute maximal doppler shift; not necessary, just directly compute upper bound for line
        // Real max_doppler_shift=1.0;
        // for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        // {
        //     for (Size pointid: points_in_grid)
        //     {
        //         const Real line_freq_inv_width=model.lines.inverse_width(pointid,lineid);
        //         //TODO: it might be faster to first compute all doppler shifts and then take min and max
        //         const Real temp_doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointid],model.geometry.rays.direction[rayid]);
        //         //compute what the min and max frequency quadratures would be for this point
        //         const Real temp_min_freq_quad=1/temp_doppler_shift*(line_freq-FREQ_QUADRATURE_EDGE/line_freq_inv_width);
        //         const Real temp_max_freq_quad=1/temp_doppler_shift*(line_freq+FREQ_QUADRATURE_EDGE/line_freq_inv_width);
        //
        //         min_freq_quad=std::min(min_freq_quad,temp_min_freq_quad);
        //         max_freq_quad=std::max(max_freq_quad,temp_max_freq_quad);
        //         max_doppler_shift=std::max(max_doppler_shift, std::abs(temp_doppler_shift))
        //     }
        // }

        //now compute the doppler shift
        Real max_abs_doppler_shift=0;
        for (Size rayid=0; rayid<parameters.nrays(); rayid++)
        {
            for (Size pointid: points_in_grid)
            {
                const Real temp_doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointid],model.geometry.rays.direction[rayid]);
                const Real curr_doppler_shift=std::abs(temp_doppler_shift-1);//doppler shift applied to center line frequency
                max_abs_doppler_shift=std::max(max_abs_doppler_shift,curr_doppler_shift);
            }
        }

        // frequency_quadrature[rayid].resize(parameters.nfreqs());
        Size freq_idx_counter=0;
        center_line_freq_indices.resize(parameters.nlines());
        for (Size lineid=0; lineid<parameters.nlines(); lineid++)
        {
            //We start with computing the maximal doppler shift and from there derive how much wing frequency quadrature points we need; err, this is per line
            //TODO: define the new frequency quadrature, centered around the center line frequency
            //TODO: compute how much limb frequencies we need
            const Real line_freq=model.lines.line[lineid];


            //compute min and max quadrature frequencies for this line frequency
            //wait, this assumes non-overlapping lines for numerical stability, or does it? If the lines do interact, the opacity and emissivity are computed exactly anyway,
            // so if we just treat the overlapping lines independently, we do not have any numerical instability; we do have wasted computation effort
            Real min_freq_quad=line_freq;
            Real max_freq_quad=line_freq;

            //For every direction, for every point, compute the quadrature bounds (without considering doppler shifts).
            for (Size rayid=0; rayid<parameters.nrays(); rayid++)
            {
                // for (Size pointid=0; pointid<parameters.npoints(); pointid++)
                for (Size point: points_in_grid)
                {
                    const Real line_freq_inv_width=model.lines.inverse_width(point,lineid);
                    //TODO: it might be faster to first compute all doppler shifts and then take min and max
                    // const Real temp_doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointid],model.geometry.rays.direction[rayid]);
                    //compute what the min and max frequency quadratures would be for this point
                    const Real temp_min_freq_quad=(line_freq-FREQ_QUADRATURE_EDGE/line_freq_inv_width);
                    const Real temp_max_freq_quad=(line_freq+FREQ_QUADRATURE_EDGE/line_freq_inv_width);
                    // const Real temp_min_freq_quad=1/temp_doppler_shift*(line_freq-FREQ_QUADRATURE_EDGE/line_freq_inv_width);
                    // const Real temp_max_freq_quad=1/temp_doppler_shift*(line_freq+FREQ_QUADRATURE_EDGE/line_freq_inv_width);

                    min_freq_quad=std::min(min_freq_quad,temp_min_freq_quad);
                    max_freq_quad=std::max(max_freq_quad,temp_max_freq_quad);
                }
            }
            //frequency difference between two succesive frequency quadrature points
            const Real line_center_freq_diff=(max_freq_quad-min_freq_quad)/(FREQ_QUADRATURE_CENTER_N_ELEMENTS-1);
            //distribute the center quadrature points within the bounds of min and max
            const Real freq_quad_delta=max_freq_quad-min_freq_quad;
            std::cout<<"freq_quad_delta: "<<freq_quad_delta<<std::endl;
            std::cout<<"inverse width: "<<model.lines.inverse_width(0,lineid)<<std::endl;
            std::cout<<"freq_quad_delta/inverse width: "<<freq_quad_delta/model.lines.inverse_width(0,lineid)<<std::endl;
            const Real line_center_inv_width=SQRT_MIN_LN_RATIO*(FREQ_QUADRATURE_CENTER_N_ELEMENTS-1)/freq_quad_delta;
            std::cout<<"line center freq diff:"<<line_center_freq_diff<<std::endl;
            std::cout<<"line center inverse width: "<<line_center_inv_width<<std::endl;

            //factors for determining the number of wing frequency quadrature points, independent of the exact line
            const Real freq_quad_spacing_squared=std::pow(line_center_freq_diff*line_center_inv_width,2);
            std::cout<<"freq quad spacing squared"<<freq_quad_spacing_squared<<std::endl;
            std::cout<<"ln (1/quad ratio)"<<std::log(1.0/MAX_FREQ_QUAD_RATIO)<<std::endl;

            const Real FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR=1.3;//1+std::sqrt(2);
            // const Real FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR=
            //   (freq_quad_spacing_squared
            //   +std::sqrt(std::log(1.0/MAX_FREQ_QUAD_RATIO)*freq_quad_spacing_squared
            //              +std::pow(freq_quad_spacing_squared,2)))
            //   /std::log(1.0/MAX_FREQ_QUAD_RATIO);//In order to sample the line wings sufficiently in the case of doppler shifts, we increase the distance between the frequencies on the line wings.
            std::cout<<"FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR: "<<FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR<<std::endl;

            const Real prefactor=(FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR-1.0)/FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR;
            const Real div_factor=std::log(FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR);

            const Real curr_line_number_wing_quadratures=std::log(prefactor*(line_freq*max_abs_doppler_shift/line_center_freq_diff)+1.0)/div_factor;
            Size number_wing_quadratures=(Size) std::ceil(curr_line_number_wing_quadratures);//and round up// note: compiler will warn against this usage, but we are taking some logarithm, so the resulting value is way smaller than the maximal integer
            const Size total_number_line_quadratures=2*number_wing_quadratures+FREQ_QUADRATURE_CENTER_N_ELEMENTS;

            //for different lines, the number of quadrature points could be different, so this is why I currently just adaptively resize these vectors
            frequency_quadrature.resize(frequency_quadrature.size()+total_number_line_quadratures);
            frequency_inverse_width_quadrature.resize(frequency_inverse_width_quadrature.size()+total_number_line_quadratures);

            //now compute the positions and inverse widths (starting from the edges of the center)
            std::vector<Real> wing_frequency_positions;
            wing_frequency_positions.resize(number_wing_quadratures);
            std::vector<Real> wing_frequency_inverse_widths;//inverse widths go up with the same factor amplification factor for each point further away from the middle
            wing_frequency_inverse_widths.resize(number_wing_quadratures);
            for (Size idx=0; idx<number_wing_quadratures; idx++)
            {
              //explicit non-recursive formula of repeatedly adding factor^(idx+1)
                const Real mult_factor=(std::pow(FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR, idx+2)
                        -FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR)/(FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR-1);
                wing_frequency_positions[idx]=mult_factor*line_center_freq_diff;
                wing_frequency_inverse_widths[idx]=line_center_inv_width/std::pow(FREQ_QUADRATURE_WING_DISTANCE_AMP_FACTOR, idx+1);
                std::cout<<"wing freq positions: "<<wing_frequency_positions[idx]<<std::endl;
                std::cout<<"wing freq inverse width: "<<wing_frequency_inverse_widths[idx]<<std::endl;
                std::cout<<"wing freq positions*center inverse width: "<<wing_frequency_positions[idx]*model.lines.inverse_width(0,lineid)<<std::endl;
            }



            //distribute the left quadrature points
            for (Size quadid=0; quadid<number_wing_quadratures; quadid++)
            {
                Real dist_to_subtract=wing_frequency_positions[number_wing_quadratures-quadid-1];//start from the left end, going towards the middle
                frequency_quadrature[freq_idx_counter]=min_freq_quad-dist_to_subtract;
                frequency_inverse_width_quadrature[freq_idx_counter]=wing_frequency_inverse_widths[number_wing_quadratures-quadid-1];
                freq_idx_counter++;
            }

            //before distributing the rest of the frequency quadratures, write down which frequency corresponds to the central line frequency
            center_line_freq_indices[lineid]=freq_idx_counter+(FREQ_QUADRATURE_CENTER_N_ELEMENTS-1)/2;

            // for (Size quadid=0; quadid<parameters.nquads(); quadid++)
            for (Size quadid=0; quadid<FREQ_QUADRATURE_CENTER_N_ELEMENTS; quadid++)
            {
                // Real frac=((Real)quadid)/(parameters.nquads()-1);//number between 0 and 1 (inclusive)
                Real frac=((Real)quadid)/(FREQ_QUADRATURE_CENTER_N_ELEMENTS-1);//number between 0 and 1 (inclusive)
                frequency_quadrature[freq_idx_counter]=min_freq_quad+freq_quad_delta*frac;
                frequency_inverse_width_quadrature[freq_idx_counter]=line_center_inv_width;
                // frequency_inverse_width_quadrature[freq_idx_counter]=SQRT_MIN_LN_RATIO*(FREQ_QUADRATURE_CENTER_N_ELEMENTS-1)/freq_quad_delta;
                // frequency_inverse_width_quadrature[freq_idx_counter]=SQRT_MIN_LN_RATIO*(parameters.nquads()-1)/freq_quad_delta;

                std::cout<<"frac: "<<frac<<std::endl;

                freq_idx_counter++;
            }

            //and distribute the right quadrature points
            for (Size quadid=0; quadid<number_wing_quadratures; quadid++)
            {
                Real dist_to_add=wing_frequency_positions[quadid];//start from the left end, going towards the middle
                frequency_quadrature[freq_idx_counter]=max_freq_quad+dist_to_add;
                frequency_inverse_width_quadrature[freq_idx_counter]=wing_frequency_inverse_widths[quadid];
                freq_idx_counter++;
            }


            for (Size freqid=0; freqid<total_number_line_quadratures; freqid++)
            {
                std::cout<<"freq quad positions-line freq: "<<frequency_quadrature[freqid]-line_freq<<std::endl;
                std::cout<<"(freq quad positions-line freq)*line inverse width: "<<(frequency_quadrature[freqid]-line_freq)*model.lines.inverse_width(0,lineid)<<std::endl;
                std::cout<<"freq quad inverse widths: "<<frequency_inverse_width_quadrature[freqid]<<std::endl;

            }

        }
    }

    // frequency_quadrature.resize(parameters.nrays());

    //if we only have a single quadrature, we might as well just place it on the line frequency (without doing any fancy calculations)
    // if (parameters.nquads()==1)
    // {
    //     for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    //     {
    //         frequency_quadrature[rayid].resize(parameters.nlines());
    //         for (Size lineid=0; lineid<parameters.nlines(); lineid++)
    //         {
    //             frequency_quadrature[rayid][lineid]=model.lines.line[lineid];
    //         }
    //     }
    // }
    // else
    // {
    //     //For every direction, for every line, compute the quadrature.
    //     for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    //     {
    //         frequency_quadrature[rayid].resize(parameters.nfreqs());
    //         Size freq_idx_counter=0;
    //         for (Size lineid=0; lineid<parameters.nlines(); lineid++)
    //         {
    //             const Real line_freq=model.lines.line[lineid];
    //             //compute min and max doppler shifts for the rays direction
    //
    //             //use point 0 as as starting guess; otherwise the min or max could be never reached might be under/overestimates
    //             const Real pt0_doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[0],model.geometry.rays.direction[rayid]);
    //
    //             Real min_freq_quad=pt0_doppler_shift*line_freq;
    //             Real max_freq_quad=pt0_doppler_shift*line_freq;
    //
    //             // Real min_raydir_doppler_shift=0.0;
    //             // Real max_raydir_doppler_shift=0.0;
    //             for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    //             {
    //                 const Real line_freq_inv_width=model.lines.inverse_width(pointid,lineid);
    //                 //TODO: it might be faster to first compute all doppler shifts and then take min and max
    //                 const Real temp_doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointid],model.geometry.rays.direction[rayid]);
    //                 //compute what the min and max frequency quadratures would be for this point
    //                 const Real temp_min_freq_quad=1/temp_doppler_shift*(line_freq-FREQ_QUADRATURE_EDGE/line_freq_inv_width);
    //                 const Real temp_max_freq_quad=1/temp_doppler_shift*(line_freq+FREQ_QUADRATURE_EDGE/line_freq_inv_width);
    //                 // const Real temp_min_freq_quad=temp_doppler_shift*line_freq-FREQ_QUADRATURE_EDGE/line_freq_inv_width;
    //                 // const Real temp_max_freq_quad=temp_doppler_shift*line_freq+FREQ_QUADRATURE_EDGE/line_freq_inv_width;
    //                 // max_raydir_doppler_shift=std::max(max_raydir_doppler_shift, temp_doppler_shift);
    //                 // min_raydir_doppler_shift=std::min(min_raydir_doppler_shift, temp_doppler_shift);
    //                 min_freq_quad=std::min(min_freq_quad,temp_min_freq_quad);
    //                 max_freq_quad=std::max(max_freq_quad,temp_max_freq_quad);
    //             }
    //
    //             //and now properly distribute the quadrature points within the bounds of min and max
    //             const Real freq_quad_delta=max_freq_quad-min_freq_quad;
    //             std::cout<<"freq_quad_delta: "<<freq_quad_delta<<std::endl;
    //             for (Size quadid=0; quadid<parameters.nquads(); quadid++)
    //             {
    //                 Real frac=((Real)quadid)/(parameters.nquads()-1);//number between 0 and 1 (inclusive)
    //                 frequency_quadrature[rayid][freq_idx_counter]=min_freq_quad+freq_quad_delta*frac;
    //                 std::cout<<"frac: "<<frac<<std::endl;
    //
    //                 freq_idx_counter++;
    //             }
    //
    //         }
    //
    //         // frequency_quadrature[rayid].resize(parameters.npoints());
    //         // for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    //         // {
    //
    //
    //         // Size freq_idx_counter=0;
    //         // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
    //         // {
    //         //     const Real line_freq=model.lines.line[lineid];
    //         //     const Real line_freq_inv_width=model.lines.inverse_width(pointid,lineid);
    //         //
    //         //     const Real min_freq_quad=min_raydir_doppler_shift*line_freq-FREQ_QUADRATURE_EDGE/line_freq_inv_width;
    //         //     const Real max_freq_quad=max_raydir_doppler_shift*line_freq+FREQ_QUADRATURE_EDGE/line_freq_inv_width;
    //         //     const Real freq_quad_delta=max_freq_quad-min_freq_quad;
    //         //     for (Size quadid=0; quadid<paramters.nquads(); quadid++)
    //         //     {
    //         //         Real frac=quadid/(parameters.nquads()-1);//number between 0 and 1 (inclusive)
    //         //         frequency_quadrature[rayid][pointid][freq_idx_counter]=min_freq_quad+freq_quad_delta*frac;
    //         //
    //         //         freq_idx_counter++;
    //         //     }
    //         // }
    //     }
    // }
    n_freq_bases=frequency_quadrature.size();
    std::cout<<"n_freq_bases: "<<n_freq_bases<<std::endl;

    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        frequencies[rayid].resize(n_points_in_grid);
        frequencies_inverse_widths[rayid].resize(n_points_in_grid);

        for (Size pointid=0; pointid<n_points_in_grid; pointid++)
        {
            const Size point_in_grid=points_in_grid[pointid];
            // vector<Real> non_shifted_frequencies=model.lines.line;
            // vector<Real> non_shifted_inverse_widths=model.lines.inverse_width[pointid];
            Real doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[point_in_grid],model.geometry.rays.direction[rayid]);
            // frequencies[rayid][pointid].resize(parameters.nlines());
            // frequencies_inverse_widths[rayid][pointid].resize(parameters.nlines());
            frequencies[rayid][pointid].resize(n_freq_bases);
            frequencies_inverse_widths[rayid][pointid].resize(n_freq_bases);

            // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
            for (Size freqid=0; freqid<n_freq_bases; freqid++)
            {
                // frequencies[rayid][pointid][lineid]=doppler_shift*model.lines.line[lineid];
                frequencies[rayid][pointid][freqid]=doppler_shift*frequency_quadrature[freqid];
                frequencies_inverse_widths[rayid][pointid][freqid]=frequency_inverse_width_quadrature[freqid]/doppler_shift;
                // std::cout<<"frequencies: "<<frequencies[rayid][pointid][freqid]<<std::endl;
                // std::cout<<"frequencies_inverse_widths: "<<frequencies_inverse_widths[rayid][pointid][freqid]<<std::endl;
                // frequencies[rayid][pointid][freqid]=doppler_shift*model.radiation.frequencies.nu[freqid];

                // Size l=model.radiation.frequencies.corresponding_l_for_spec[freqid];   ///< number of line species corresponding to frequency
                // Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];   ///< number of transition corresponding to frequency
                // Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];   ///< number of line number corresponding to frequency
                // // std::cout<<"l, k, z: "<<l<<", "<<k<<", "<<z<<std::endl;
                // Size corresp_lineid=model.lines.line_index(l, k);
                // // model.lines.lineProducingSpecies[l].nr_line[pointid][k][z];
                // // std::cout<<"corresp lineid: "<<corresp_lineid<<std::endl;
                // if (parameters.nquads()==1)
                // {//no better estimate availabe, so just using the default inverse width
                //     frequencies_inverse_widths[rayid][pointid][freqid]=model.lines.inverse_width(pointid,corresp_lineid)/doppler_shift;
                // }
                // else
                // {
                //     // frequencies_inverse_widths[rayid][pointid][freqid]=frequency_inverse_width_quadrature[corresp_lineid]/doppler_shift;
                //     frequencies_inverse_widths[rayid][pointid][freqid]=frequency_inverse_width_quadrature[freqid]/doppler_shift;
                //     std::cout<<"frequencies_inverse_widths: "<<frequencies_inverse_widths[rayid][pointid][freqid]<<std::endl;
                //     std::cout<<"comparison inverse width:"<<model.lines.inverse_width(pointid,corresp_lineid)/doppler_shift<<std::endl;
                // }
                // // frequencies_inverse_widths[rayid][pointid][freqid]=model.lines.inverse_width(pointid,corresp_lineid)/doppler_shift;
            }
        }
    }

    // non_doppler_shifted_frequencies.resize(parameters.nlines());
    // non_doppler_shifted_inverse_widths.resize(parameters.nlines());
    // for (Size lineid=0; lineid<parameters.nlines(); lineid++)
    // {
    //     non_doppler_shifted_frequencies[lineid]=model.lines.line[lineid];
    //     non_doppler_shifted_inverse_widths[lineid].resize(parameters.npoints());
    //     for (Size pointid=0; pointid<parameters.npoints(); pointid++)
    //     {
    //         non_doppler_shifted_inverse_widths[lineid][pointid]=model.lines.inverse_width(pointid,lineid);
    //     }
    // }

    // non_doppler_shifted_frequencies.resize(parameters.nfreqs());
    non_doppler_shifted_frequencies.resize(n_freq_bases);
    for (Size freqid=0; freqid<n_freq_bases; freqid++)
    {
        non_doppler_shifted_frequencies[freqid]=frequency_quadrature[freqid];
        std::cout<<"non doppler shifted frequencies: "<<non_doppler_shifted_frequencies[freqid]<<std::endl;
    }
    // }

    //only for the other frame formulation; FIXME: delete other non-used frame formulation
    non_doppler_shifted_inverse_widths.resize(n_freq_bases);
    for (Size freqid=0; freqid<n_freq_bases; freqid++)
    {
        non_doppler_shifted_inverse_widths[freqid]=frequency_inverse_width_quadrature[freqid];

        // // non_doppler_shifted_frequencies[freqid]=model.radiation.frequencies.nu[freqid];
        // // std::cout<<"non doppler shifted frequencies: "<<non_doppler_shifted_frequencies[freqid]<<std::endl;
        // // non_doppler_shifted_inverse_widths[freqid].resize(parameters.npoints());
        // non_doppler_shifted_inverse_widths[freqid].resize(n_point_bases[freqid]);
        // for (Size pointid=0; pointid<parameters.npoints(); pointid++)
        // {
        //     Size l=model.radiation.frequencies.corresponding_l_for_spec[freqid];   ///< number of line species corresponding to frequency
        //     Size k=model.radiation.frequencies.corresponding_k_for_tran[freqid];   ///< number of transition corresponding to frequency
        //     Size z=model.radiation.frequencies.corresponding_z_for_line[freqid];   ///< number of line number corresponding to frequency
        //     // std::cout<<"l, k, z: "<<l<<", "<<k<<", "<<z<<std::endl;
        //     Size corresp_lineid=model.lines.line_index(l, k);
        //     // std::cout<<"corresp lineid: "<<corresp_lineid<<std::endl;
        //     non_doppler_shifted_inverse_widths[freqid][pointid]=model.lines.inverse_width(pointid,corresp_lineid);
        //     // std::cout<<"inverse_width: "<<non_doppler_shifted_inverse_widths[freqid][pointid]<<std::endl;
        // }
    }


    //Now determine which basis's interact with eachother
    prune_positional_bases(model);
    set_interacting_bases(model);

    // Also do not forget to initialize the basis coefficients // TODO: figure out which exact format (std::vector, Eigen, ...) to use
    // basis_coefficients=VectorXr::Zero(parameters.nrays()*parameters.nlines()*parameters.npoints());
    const Size total_nb_bases=parameters.nrays()*cum_n_point_bases[n_freq_bases];
    basis_coefficients=VectorXr::Zero(total_nb_bases);
    // basis_renorm_factor=VectorXr::Zero(total_nb_bases);//renormalizes the bases; should just be interpreted as a diagonal matrix for most purposes
    basis_renorm_factor.resize(total_nb_bases,total_nb_bases);
    basis_renorm_factor.reserve(Eigen::VectorXi::Constant(total_nb_bases,1));
    // basis_coefficients=VectorXr::Zero(parameters.nrays()*parameters.nfreqs()*parameters.npoints());

    //err, it turns out I first need to determine the typical rbf radius in order to determine whether we need a slope...
    // So that is why this seems out of place.
    //TODO: when finally deciding to remove the observer frame formulation, please try to remove the ray dependence
    rbf_using_slope.resize(parameters.nrays());
    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        rbf_using_slope[rayid].resize(n_freq_bases);
        std::cout<<"n freq bases: "<<n_freq_bases<<std::endl;
        //use basis_point_idx; switch order of for loop
        // for (Size freqid=0; freqid<parameters.nfreqs(); freqid++)
        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {
            std::cout<<"n point bases: "<<n_point_bases[freqid]<<std::endl;
            rbf_using_slope[rayid][freqid].resize(n_point_bases[freqid]);

            // for (Size pointid=0; pointid<parameters.npoints(); pointid++)
            for (Size basis_point_id=0; basis_point_id<n_point_bases[freqid]; basis_point_id++)
            {
                const Real basis_radius=rbf_radii[freqid][basis_point_id];
                const Size basis_center_id=local_basis_point_idx[freqid][basis_point_id];
                const Real point_opacity=get_opacity(model, rayid, freqid, basis_center_id);
                //TODO: maybe do something else than a binary choice between slope and no slope
                if (point_opacity*basis_radius>SLOPE_THRESHOLD)
                {
                    rbf_using_slope[rayid][freqid][basis_point_id]=true;//forcefully enabling slope
                    // rbf_using_slope[rayid][freqid][basis_point_id]=false;
                }
                else
                {
                    rbf_using_slope[rayid][freqid][basis_point_id]=true;
                    // rbf_using_slope[rayid][freqid][basis_point_id]=false; //disabling the slope
                }
                std::cout<<"rid, fid, pid: "<<rayid<<", "<<freqid<<", "<<basis_point_id<<std::endl;
                std::cout<<"using slope: "<<rbf_using_slope[rayid][freqid][basis_point_id]<<std::endl;
            }
        }
    }

    // set_eval_matrix(model);
}


inline bool Collocation :: freq_cutoff_condition(Real abs_freq_diff, Real freq_basis_inverse_width)
{
    return false;//debug
    return (abs_freq_diff*freq_basis_inverse_width>TRUNCATION_SIGMA);
}

/// The basis function associated with each direction
inline Real Collocation :: basis_direction(Size basis_ray_id)
{
    //as we do not need directional derivatives, piecewise constant basis functions are sufficient;
    return 1;
}

/// The basis function associated with each direction integrated over the solid angle
inline Real Collocation :: basis_direction_int(Geometry& geometry, Size basis_ray_id)
{
    return FOUR_PI*geometry.rays.weight[basis_ray_id];
}

///  The basis function associated with each frequency, evaluated at frequency currfreq
inline Real Collocation :: basis_freq(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Real currfreq)
{
    //TODO: figure out what to use exactly
    //also implement some logical data structure (possibly first read out all the sorted frequencies, then construct the basis functions)
    Real abs_freq_diff=0;
    Real freq_inv_width=0;
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];

    // Real background=boundary_intensity(Model& model, Size rayidx, Size freqidx, Size pointidx);
    if (USING_COMOVING_FRAME)
    {
        freq_inv_width=non_doppler_shifted_inverse_widths[basis_freq_id];
        abs_freq_diff=std::abs(currfreq-non_doppler_shifted_frequencies[basis_freq_id]-FREQ_OFFSET/freq_inv_width);
        // std::cout<<"freq inverse width: "<<freq_inv_width<<std::endl;
        // std::cout<<"shift: "<<FREQ_OFFSET/freq_inv_width<<std::endl;
        // std::cout<<"base freq: "<<non_doppler_shifted_frequencies[freqidx]<<std::endl;
        // std::cout<<"abs_freq_diff: "<<abs_freq_diff<<std::endl;
        // std::cout<<"freqidx: "<<freqidx<<std::endl;
    }
    else
    {
        abs_freq_diff=std::abs(currfreq-frequencies[basis_ray_id][center_point_id][basis_freq_id]);
        freq_inv_width=frequencies_inverse_widths[basis_ray_id][center_point_id][basis_freq_id];
    }

    //FIXME: freq offset somewhat makes the cutoff condition a bit nonsymmetrical
    if (freq_cutoff_condition(abs_freq_diff,freq_inv_width))
    {
        return 0;
    }
    else
    {
        // std::cout<<"basis_freq: "<<INVERSE_SQRT_PI*frequencies_inverse_widths[rayidx][pointidx][lineidx]*std::exp(-std::pow(abs_freq_diff*frequencies_inverse_widths[rayidx][pointidx][lineidx],2))<<std::endl;
        return INVERSE_SQRT_PI*freq_inv_width*std::exp(-std::pow(abs_freq_diff*freq_inv_width,2));
    }
}

inline Real Collocation :: basis_freq_lp_int(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Real lp_freq, Real lp_freq_inv_width)
{
    Real point_freq=0;//frequencies[rayidx][pointidx][lineidx];
    Real point_freq_inv_width=0;//frequencies_inverse_widths[rayidx][pointidx][lineidx];
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];

    // Real abs_freq_diff=0;
    // Real freq_inv_width=0;
    if (USING_COMOVING_FRAME)
    {
        point_freq_inv_width=non_doppler_shifted_inverse_widths[basis_freq_id];
        point_freq=non_doppler_shifted_frequencies[basis_freq_id]+FREQ_OFFSET/point_freq_inv_width;
    }
    else
    {
        point_freq=frequencies[basis_ray_id][center_point_id][basis_freq_id];
        point_freq_inv_width=frequencies_inverse_widths[basis_ray_id][center_point_id][basis_freq_id];
    }

    const Real inv_variance=1/(1/std::pow(point_freq_inv_width,2)+1/std::pow(lp_freq_inv_width,2));// note: still includes the factor (1/sqrt(2)) squared (because variance)
    const Real delta_freq=lp_freq-point_freq;

    // std::cout<<"rel_delta_freq: "<<point_freq_inv_width*delta_freq<<std::endl;
    // std::cout<<"freq int: "<<INVERSE_SQRT_PI*std::sqrt(2)*std::sqrt(inv_variance)*std::exp(-(std::pow(delta_freq,2)*inv_variance))<<std::endl;
    //NOTE: factor *std::sqrt(2) is maybe necessary, but i don't know for sure... TODO: recalculate this analytically
    return INVERSE_SQRT_PI*std::sqrt(inv_variance)*std::exp(-(std::pow(delta_freq,2)*inv_variance));

}




///  The derivative of the frequency basis function
inline Real Collocation :: basis_freq_der(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Real currfreq)
{
    Real abs_freq_diff=0;
    Real freq_diff=0;
    Real freq_inv_width=0;
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];

    if (USING_COMOVING_FRAME)
    {
        freq_inv_width=non_doppler_shifted_inverse_widths[basis_freq_id];//TEST:/2 for testing wider function
        freq_diff=currfreq-non_doppler_shifted_frequencies[basis_freq_id]-FREQ_OFFSET/freq_inv_width;
        // abs_freq_diff=std::abs(currfreq-non_doppler_shifted_frequencies[freqidx]-FREQ_OFFSET/freq_inv_width);
    }
    else
    {
        freq_diff=currfreq-frequencies[basis_ray_id][center_point_id][basis_freq_id];
        // abs_freq_diff=std::abs(currfreq-frequencies[rayidx][pointidx][freqidx]);
        freq_inv_width=frequencies_inverse_widths[basis_ray_id][center_point_id][basis_freq_id];
    }
    abs_freq_diff=std::abs(freq_diff);//TODO: maybe replace def of freq_cutoff_condition by accepting the freq diff instead of the abs freq diff

    if (freq_cutoff_condition(abs_freq_diff,freq_inv_width))
    {
        return 0;
    }
    else
    {
        return INVERSE_SQRT_PI*freq_inv_width*-2*(freq_diff*freq_inv_width)*freq_inv_width*std::exp(-std::pow(abs_freq_diff*freq_inv_width,2));
    }
}

//DEPRECATED
// ///  Function to manipulate the (already scaled) distance
// inline Real Collocation :: distance_manip(Size pointid, Real original_distance)
// {
//     //test option to check that everything works correctly7
//     // return original_distance;
//     //option to scale by using the n-th root
//     // return std::pow(original_distance, DISTANCE_EXPONENT);
//     //option to scale by using a truncated logarithm TODO: define eps for each point, (look at closest distance)
//     // if (original_distance<DISTANCE_EPS)
//     // {return 0;}
//     // else
//     // {return 1-std::log(original_distance)/std::log(DISTANCE_EPS);}
//     //truncated logarithm fitted exactly into the domain and range [0,1]
//     return (1-std::log(original_distance+DISTANCE_EPS)/std::log(DISTANCE_EPS))/(1-std::log(1+DISTANCE_EPS)/std::log(DISTANCE_EPS));
// }
//
// inline Real Collocation :: distance_manip_der(Size pointid, Real original_distance)
// {
//     //test option to check that everything works correctly
//     // return 1;
//     //option to scale by using the n-th root
//     // return DISTANCE_EXPONENT*std::pow(original_distance, DISTANCE_EXPONENT-1);//TODO: check if nans occur if using pow in this way (negative exponent)
//     //option to scale by using a truncated logarithm TODO: define eps for each point, (look at closest distance)
//     //if (original_distance<eps)
//     //{return 0;}
//     //else
//     //{return -1/(original_distance*std::log(eps);}
//     //truncated logarithm fitted exactly into the domain and range [0,1]
//     return -1/(original_distance+DISTANCE_EPS)/std::log(DISTANCE_EPS)/(1-std::log(1+DISTANCE_EPS)/std::log(DISTANCE_EPS));
//
// }

///Returns the slope
///   @param[in]: relative distance on ray compared to the rbf radius (should lie between -1 and 1, inclusive)
//TODO: also add the other slope candidate(s)
inline Real Collocation :: get_slope(Real x, Size basis_ray_id, Size basis_point_id, Size basis_freq_id)
{
   if (!rbf_using_slope[basis_ray_id][basis_freq_id][basis_point_id])
   {
      // std::cout<<"not using slope"<<std::endl;
      return 1;
   }
   else
   {
      // A sinus as slope
      // Real slope=1+std::sin(PI*x);
      // The slope 1+x/SLOPE_FACTOR
      // Real slope=1+x/SLOPE_FACTOR;
      // The sigmoid-inspired slope
      Real slope=1.0/2.0+1.0/(1.0+std::exp(-SLOPE_STEEPNESS*x));
      // A slope which has derivative 0 at the ends
      // Real slope=1.0+3.0/2.0*x-std::pow(x,3)/2.0;
      return slope;
   }
}

///Returns the slope derivative
///   @param[in]: relative distance on ray compared to the rbf radius (should lie between -1 and 1, inclusive)
///NOTE: the usual factor /radius should still be multiplied manually
//TODO: also add the other slope candidate(s)
inline Real Collocation :: get_slope_derivative(Real x, Size basis_ray_id, Size basis_point_id, Size basis_freq_id)
{
    if (!rbf_using_slope[basis_ray_id][basis_freq_id][basis_point_id])
    {
        // std::cout<<"not using slope"<<std::endl;
        return 0;
    }
    else
    {
        // A sinus as slope
        // Real slope_derivative=PI*std::cos(PI*x);
        // the slope 1+x/SLOPE_FACTOR
        // Real slope_derivative=1/SLOPE_FACTOR;
        // The sigmoid-inspired slope
        Real exp_factor=std::exp(-SLOPE_STEEPNESS*x);
        Real slope_derivative=SLOPE_STEEPNESS*exp_factor/std::pow(1+exp_factor,2);
        // A slope which has derivative 0 at the ends
        // Real slope_derivative=3.0/2.0-3.0/2.0*std::pow(x,2);
        return slope_derivative;
    }

}

// Multiplies a Vector3D with a Real
inline Vector3D multiply_vector3D_real(Vector3D vector3d, Real real)
{
    Vector3D toreturn=Vector3D(real*vector3d.x(), real*vector3d.y(), real*vector3d.z());
    return toreturn;
}

///  The radial basis function for the position (perpendicular to the ray)
inline Real Collocation :: basis_point_perp(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry)
{
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];
    Vector3D diff_vector=location-point_locations[center_point_id];
    Vector3D raydirection=geometry.rays.direction[basis_ray_id];
    Vector3D ray_component=multiply_vector3D_real(raydirection, raydirection.dot(diff_vector));
    Vector3D diff_vector_minus_ray_component=diff_vector-ray_component;
    Real radius=rbf_radii[basis_freq_id][basis_point_id];
    Real perp_distance=std::sqrt(diff_vector_minus_ray_component.dot(diff_vector_minus_ray_component))/radius;

    //TODO: add basis function, fill in perp_distance (normed to 1)

    return 1;//TODO: add better basis function (lazy as it doesn't matter in 1D), in 3D it definitely does ()
}

///  The radial basis function for the position (symmetric to the ray)
inline Real Collocation :: basis_point_symm(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry)
{
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];
    Vector3D diff_vector=location-point_locations[center_point_id];
    Vector3D raydirection=geometry.rays.direction[basis_ray_id];
    // std::cout<<"diff vector.x: "<<diff_vector.x()<<std::endl;
    // std::cout<<"raydirection.x: "<<raydirection.x()<<std::endl;
    // Vector3D ray_component=multiply_vector3D_real(raydirection, raydirection.dot(diff_vector));
    Real radius=rbf_radii[basis_freq_id][basis_point_id];
    // Real ray_distance=std::sqrt(ray_component.dot(ray_component))/radius;
    Real ray_distance=raydirection.dot(diff_vector)/radius;
    Real abs_ray_distance=std::abs(ray_distance);

    // std::cout<<"ray_distance: "<<ray_distance<<std::endl;
    // std::cout<<"abs_ray_distance: "<<abs_ray_distance<<std::endl;

    bool switch_sign=false;//true if ray_distance negative, false if positive; just for switching the sign of the integral
    if (ray_distance<0)
    {
        switch_sign=true;
    }

    //TODO: add basis function, fill in ray_distance

    //Integral of: Wendland basis function (3,1), but made symmetric around 0
    Real mid_value=1.0/3.0;//actually double, to shift the integral
    Real integral=-(4.0*abs_ray_distance+1.0)/5.0*std::pow(1.0-abs_ray_distance,5)-2.0/15.0*std::pow(1.0-abs_ray_distance,6)+1.0/3.0;//with constant to make zero at x=0
    if (switch_sign)//in this case, maybe using std::copysign makes sense...
    {
        integral=-integral;
    }
    std::cout<<"integral: "<<integral<<std::endl;
    std::cout<<"switch_sign: "<<switch_sign<<std::endl;
    std::cout<<"sum: "<<mid_value+integral<<std::endl;

    return mid_value+integral;

    // //Wendland basis function (3,1), but made symmetric around 0
    // Real basis_wend=std::pow(1.0-abs_ray_distance,4)*(1.0+4.0*abs_ray_distance);
    // //Some polynomial function such that the derivative is nonzero at 0
    // Real basis_pol=-3.0*std::pow(abs_ray_distance,5)+8.0*std::pow(abs_ray_distance,4)-6.0*std::pow(abs_ray_distance,3)+abs_ray_distance;
    // if (ray_distance<0)
    // {
    //     basis_pol=-basis_pol;
    // }
    //
    // return 0.02*basis_wend+basis_pol;

    // Real slope=get_slope(ray_distance, rayindex, centerpoint, freqidx);
    // return slope*basis;
}

///  The derivative of the radial basis function for the position (symmetric to the ray)
inline Real Collocation :: basis_point_symm_der(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry)
{
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];
    Vector3D diff_vector=location-point_locations[center_point_id];
    Vector3D raydirection=geometry.rays.direction[basis_ray_id];
    // Vector3D ray_component=multiply_vector3D_real(raydirection, raydirection.dot(diff_vector));
    Real radius=rbf_radii[basis_freq_id][basis_point_id];
    Real ray_distance=raydirection.dot(diff_vector)/radius;
    Real abs_ray_distance=std::abs(ray_distance);

    //Wendland basis function (3,1), but made symmetric around 0
    Real basis=std::pow(1.0-abs_ray_distance,4)*(1.0+4.0*abs_ray_distance);

    return basis/radius;

    // // Real basis=std::pow(1.0-abs_ray_distance,4)*(1.0+4.0*abs_ray_distance);
    // Real basis_wend_der=-4*std::pow(1.0-abs_ray_distance,3)*(1.0+4.0*abs_ray_distance)+4.0*std::pow(1.0-abs_ray_distance,4);
    // Real basis_pol_der=-15*std::pow(abs_ray_distance,4)+32.0*std::pow(abs_ray_distance,3)-18.0*std::pow(abs_ray_distance,2)+1;
    // if (ray_distance<0)//switch sign due to symmetrizing the basis
    // {
    //    basis_wend_der=-basis_wend_der;
    // }
    // return (0.02*basis_wend_der+basis_pol_der)/radius;

    // Real slope=get_slope(ray_distance, rayindex, centerpoint, freqidx);
    // Real slope_der=get_slope_derivative(ray_distance, rayindex, centerpoint, freqidx);
    // return (slope*basis_der+slope_der*basis)/radius;
}


///  The radial basis function for the position
inline Real Collocation :: basis_point(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry)
{
    // return basis_point_symm(centerpoint, location, rayindex, freqidx, geometry)*basis_point_perp(centerpoint, location, rayindex, freqidx, geometry);
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];
    Vector3D diff_vector=location-point_locations[center_point_id];
    Real distance=std::sqrt(diff_vector.dot(diff_vector));
    Real radius=rbf_radii[basis_freq_id][basis_point_id];
    // std::cout<<"distance: "<<distance<<std::endl;
    // std::cout<<"radius: "<<radius<<std::endl;
    //if the location is too far from the center, the evalutation of the compact radial basis function is 0.
    if (distance>=radius)
    {
        return 0;
    }
    else
    {
        Real rel_dist=distance/radius;
        // Real manip_rel_dist=distance_manip(centerpoint, rel_dist);
        // std::cout<<"rel_dist: "<<rel_dist<<std::endl;
        //This basis function has nice continuous derivatives at r=0 and r=1 (being 0 at both places);
        Real basis=-4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        // Real basis=1-rel_dist;
        // Real basis=std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));
        Vector3D raydirection=geometry.rays.direction[basis_ray_id];

        // Real x=0;
        // if (rel_dist>0)
        // {
        //     Real costheta=raydirection.dot(diff_vector)/std::sqrt(diff_vector.squaredNorm());
        //     Real x=costheta*manip_rel_dist;
        // }
        Real x=raydirection.dot(diff_vector)/radius;
        //the slope 1+x in the direction of the ray;
        // Real slope=1+x/SLOPE_FACTOR;
        //the slope 1/2+1/(1+e^(-kx))
        // Real slope=1.0/2.0+1.0/(1.0+std::exp(-SLOPE_STEEPNESS*x));
        // std::cout<<"slope: "<<slope<<std::endl;
        // Real slope=1;
        Real slope=get_slope(x, basis_ray_id, basis_point_id, basis_freq_id);
        return basis*slope;
        // return basis;
        // return 1-rel_dist;
    }
}

///  The directional derivative of the position basis function
inline Real Collocation :: basis_point_der(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry)
{
    // return basis_point_symm_der(centerpoint, location, rayindex, freqidx, geometry)*basis_point_perp(centerpoint, location, rayindex, freqidx, geometry);
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];
    Vector3D diff_vector=location-point_locations[center_point_id];
    Real distance=std::sqrt(diff_vector.dot(diff_vector));
    // Real distance=distance_manip(centerpoint, std::sqrt(diff_vector.dot(diff_vector)));
    // Real manip_der_factor=distance_manip_der(centerpoint, std::sqrt(diff_vector.dot(diff_vector)));
    Real radius=rbf_radii[basis_freq_id][basis_point_id];
    //if the location is too far from the center, the evalutation of the compact radial basis function is 0.
    if (distance>=radius)
    {
        return 0;
    }
    else if (distance==0)
    {
        // std::cout<<"distance is zero"<<std::endl;
        Real rel_dist=distance/radius;
        // std::cout<<"manip der factor: "<<manip_der_factor<<std::endl;
        // Real manip_rel_dist=distance_manip(centerpoint, rel_dist);
        // Real manip_der_factor=distance_manip_der(centerpoint, 0);

        // return -2/radius;//for 1-2x+x^2 as basis function
        // Real basis=1-rel_dist;//for 1-x as basis function
        // return 0; //derivative at r=0 is defined to be zero

        // the slope 1+x/SLOPE_FACTOR
        // Real slope_derivative=1/SLOPE_FACTOR;

        // the slope 1/2+1/(1+e^(-kx))
        // Real exp_factor=std::exp(-SLOPE_STEEPNESS*x);
        // Real slope_derivative=SLOPE_STEEPNESS*std::exp(-SLOPE_STEEPNESS*x)/std::pow(1+std::exp(-SLOPE_STEEPNESS*x),2);
        //evaluated at x=0
        // Real slope_derivative=SLOPE_STEEPNESS/4;
        Real slope_derivative=get_slope_derivative(0, basis_ray_id, basis_point_id, basis_freq_id);
        // Real slope_derivative=0;
        Real basis=-4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        // Real basis=std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));
        return slope_derivative*basis/radius;
        // return manip_der_factor*slope_derivative*basis/radius;
    }
    else
    {
        Real rel_dist=distance/radius;
        // std::cout<<"rel_distance: "<<rel_dist<<std::endl;
        // Real manip_der_factor=distance_manip_der(centerpoint, rel_dist);
        // std::cout<<"manip der factor: "<<manip_der_factor<<std::endl;
        // Real manip_rel_dist=distance_manip(centerpoint, rel_dist);

        Vector3D raydirection=geometry.rays.direction[basis_ray_id];
        Real costheta=raydirection.dot(diff_vector)/std::sqrt(diff_vector.squaredNorm());
        // Real x=raydirection.dot(diff_vector)/radius;
        Real x=costheta*rel_dist;
        //the slope 1+x in the direction of the ray;
        // Real slope=1+x/SLOPE_FACTOR;
        // Real slope_derivative=1/SLOPE_FACTOR;
        //the slope 1/2+1/(1+e^(-kx)) in direction of the ray
        // Real slope=1.0/2.0+1.0/(1.0+std::exp(-SLOPE_STEEPNESS*x));
        // Real exp_factor=std::exp(-SLOPE_STEEPNESS*x);
        // Real slope_derivative=SLOPE_STEEPNESS*exp_factor/std::pow(1+exp_factor,2);
        Real slope=get_slope(x, basis_ray_id, basis_point_id, basis_freq_id);
        Real slope_derivative=get_slope_derivative(x, basis_ray_id, basis_point_id, basis_freq_id);

        // Real slope=1;
        // Real slope_derivative=0;
        Real basis=-4.0/(1.0+std::pow(rel_dist,3))+6.0/(1.0+std::pow(rel_dist,2))-1;
        // Real basis=1-rel_dist;
        // Real basis=std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));

        // std::cout<<"raydirection=[x,y,z]: "<<raydirection.x()<<", "<<raydirection.y()<<", "<<raydirection.z()<<std::endl;
        //The derivatives in the perpendicular direction (to the radial direction) are 0, so we can calculate the derivative
        // by taking the radial derivative and multiplying it by cos(theta) in which theta is the angle between the radial and the actual direction
        //Below, we see the simple rule for calculating the cosine
        //Note: The ray directions are normalised, so we can omit the normalization of the ray direction
        // Minus sign because we actually need the gradient (not -1 times the gradient); this is due to the definition of the diff_vector, which is natural for RBF's but does actually not give us the angle we want

        //derivative of -4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        Real radial_derivative = 12*(std::pow(rel_dist,2))/(std::pow(1.0+std::pow(rel_dist,3),2))-12*rel_dist/(std::pow(1.0+std::pow(rel_dist,2),2));
        // Real slope_derivative = 1;
        //derivative of 1-rel_dist
        // Real radial_derivative = -1;
        //derivative of 1-2rel_dist+rel_dist^2
        // Real radial_derivative = -2+2*rel_dist;
        //derivative of e^-(x*trunc)^2
        // Real radial_derivative = -std::pow(TRUNCATION_SIGMA,2)*2*rel_dist*std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));
        // Because we switched the perspective of the basis and point, we need the minus sign
        // Normally you take a basis, and compute the derivative at a given point; in this case, we do the opposite:
        // we take a point and calculate the derivative at that position in some basis functions.
        // std::cout<<"der without slope: "<<costheta*radial_derivative/radius<<std::endl;
        // std::cout<<"der with slope: "<<costheta*radial_derivative*slope/radius+slope_derivative*basis/radius<<std::endl;
        //currently, the slope uses the default distance, not the manipulated one
        return costheta*radial_derivative*slope/radius+slope_derivative*basis/radius;
        // return costheta*radial_derivative/radius;
    }
}

///Second order derivative (does not support a slope)
///  The directional derivative of the position basis function
inline Real Collocation :: basis_point_der2(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, Vector3D& location, Geometry& geometry)
{
    // return basis_point_symm_der(centerpoint, location, rayindex, freqidx, geometry)*basis_point_perp(centerpoint, location, rayindex, freqidx, geometry);
    const Size center_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];
    Vector3D diff_vector=location-point_locations[center_point_id];
    Real distance=std::sqrt(diff_vector.dot(diff_vector));
    // Real distance=distance_manip(centerpoint, std::sqrt(diff_vector.dot(diff_vector)));
    // Real manip_der_factor=distance_manip_der(centerpoint, std::sqrt(diff_vector.dot(diff_vector)));
    Real radius=rbf_radii[basis_freq_id][basis_point_id];
    //if the location is too far from the center, the evalutation of the compact radial basis function is 0.
    if (distance>=radius)
    {
        return 0;
    }
    else if (distance==0)//TODO: think about removing this, and instead setting costheta to 1... (same effect, less branching?)
    {
        // std::cout<<"distance is zero"<<std::endl;
        // basis -4.0/(1.0+std::pow(rel_dist,3))+6.0/(1.0+std::pow(rel_dist,2))-1;
        return -12/std::pow(radius,2);

        // Real rel_dist=distance/radius;

        // Real slope_derivative=get_slope_derivative(0, rayindex, centerpoint, freqidx);
        // // Real slope_derivative=0;
        // Real basis=-4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        // // Real basis=std::exp(-std::pow(rel_dist*TRUNCATION_SIGMA,2));
        // return slope_derivative*basis/std::pow(radius,2);
    }
    else
    {
        Real rel_dist=distance/radius;
        // std::cout<<"rel_distance: "<<rel_dist<<std::endl;
        // Real manip_der_factor=distance_manip_der(centerpoint, rel_dist);
        // std::cout<<"manip der factor: "<<manip_der_factor<<std::endl;
        // Real manip_rel_dist=distance_manip(centerpoint, rel_dist);

        Vector3D raydirection=geometry.rays.direction[basis_ray_id];
        Real costheta=raydirection.dot(diff_vector)/std::sqrt(diff_vector.squaredNorm());
        // Real x=raydirection.dot(diff_vector)/radius;
        Real x=costheta*rel_dist;

        // Real slope=get_slope(x, rayindex, centerpoint, freqidx);
        // Real slope_derivative=get_slope_derivative(x, rayindex, centerpoint, freqidx);

        // std::cout<<"raydirection=[x,y,z]: "<<raydirection.x()<<", "<<raydirection.y()<<", "<<raydirection.z()<<std::endl;
        //The derivatives in the perpendicular direction (to the radial direction) are 0, so we can calculate the derivative
        // by taking the radial derivative and multiplying it by cos(theta) in which theta is the angle between the radial and the actual direction

        //2nd derivative of -4/(1+std::pow(rel_dist,3))+6/(1+std::pow(rel_dist,2))-1;
        const Real oneplusr2=1.0+std::pow(rel_dist,2);
        const Real oneplusr3=1.0+std::pow(rel_dist,3);
        Real radial_derivative = 12*(std::pow(rel_dist,2))/(std::pow(oneplusr3,2))-12*rel_dist/(std::pow(oneplusr2,2));
        Real radial_derivative2 = 24.0*rel_dist/std::pow(oneplusr3,2)-72.0*std::pow(rel_dist,4)/std::pow(oneplusr3,3)
                                 -12.0/std::pow(oneplusr2,2)+48.0*std::pow(rel_dist,2)/std::pow(oneplusr2,3);

        return std::pow(costheta,2)*radial_derivative2/std::pow(radius,2)+radial_derivative*((1.0-std::pow(costheta,2))/rel_dist)/std::pow(radius,2);
        // TODO: if you do multiply with some function, recompute the second derivative
        // return std::pow(costheta,2)*radial_derivative*slope/radius+slope_derivative*basis/radius;
        // return costheta*radial_derivative/radius;
    }
}


/// Fills the triplets which correspond to the point (in (raydirection, frequency, position) space) at a given triplet location
///   @param[in] pointid:
///   @param[in/out] basis_triplets_to_fill: An empty set which will contain the basis triplets after execution
inline void Collocation :: get_nonzero_basis_triplets(std::set<std::vector<Size>>& basis_triplets_to_fill, Size rayid, Size freqid, Size pointid)
{
    // std::set<std::vector<Size>> basis_triplets;
    // I currently assume that the different ray direction bases do not interact with eachother
    // std::cout<<"n_interacting_freqs: "<<other_frequency_basis_interacting_with[pointid][freqid].size()<<std::endl;
    for (Size interacting_freq: other_frequency_basis_interacting_with[pointid][freqid])
    // for (Size interacting_points: other_radial_basis_interacting_with[pointidx])
    {
        // std::cout<<"n_interacting_points: "<<radial_basis_functions_having_evaluation_at[interacting_freq][pointid].size()<<std::endl;
        for (Size interacting_position_basis_id: radial_basis_functions_having_evaluation_at[interacting_freq][pointid])
        // for (Size interacting_freqs: other_frequency_basis_interacting_with[pointidx][freqidx])
        {   //TODO: we are constantly allocating new memory; this might be slow; better to just compute this once...
            vector<Size> temp_triplet{rayid, interacting_freq, interacting_position_basis_id};
            basis_triplets_to_fill.insert(temp_triplet);
        }
    }
    // return basis_triplets;
}

// ///  Returns the matrix index corresponding to the general index (rayidx, lineidx, pointidx)
// inline Size Collocation :: get_mat_index(Size rayidx, Size freqidx, Size pointidx)
// {
//     // return rayidx*parameters.nlines()*point_locations.size()+lineidx*point_locations.size()+pointidx;
//     // return rayidx*parameters.nfreqs()*point_locations.size()+freqidx*point_locations.size()+pointidx;
//     return rayidx*parameters.nfreqs()*point_locations.size()+pointidx*parameters.nfreqs()+freqidx;
// }

/// TODO: number of point basis functions changes with frequency, so change order of point and frequency
///  Returns the matrix index corresponding to the general index (rayidx, lineidx, pointidx)//for basis function index
inline Size Collocation :: get_mat_col_index(Size basis_ray_id, Size basis_freq_id, Size basis_point_id)
{
    // return rayidx*parameters.nlines()*point_locations.size()+lineidx*point_locations.size()+pointidx;
    // return rayidx*parameters.nfreqs()*point_locations.size()+freqidx*point_locations.size()+pointidx;
    // return rayidx*parameters.nfreqs()*point_locations.size()+pointidx*parameters.nfreqs()+freqidx;
    return basis_ray_id*cum_n_point_bases[n_freq_bases]+cum_n_point_bases[basis_freq_id]+basis_point_id;
    //TODO: replace nfreqs() with number frequency quadratures
}

///  Returns the matrix index corresponding to the general index (rayidx, lineidx, pointidx)//for evaluation point index
//TOOD: for non-least squares formulation, this might need to change (or just use get_mat_col_index...)
inline Size Collocation :: get_mat_row_index(Size rayid, Size freqid, Size pointid)
{
    // return rayidx*parameters.nlines()*point_locations.size()+lineidx*point_locations.size()+pointidx;
    // return rayidx*parameters.nfreqs()*point_locations.size()+freqidx*point_locations.size()+pointidx;
    // return rayidx*parameters.nfreqs()*parameters.npoints()+eval_pointidx*parameters.nfreqs()+freqidx;
    return rayid*n_freq_bases*n_points_in_grid+freqid*n_points_in_grid+pointid;
    //TODO: replace nfreqs() with number frequency quadratures
}

///  Returns the matrix index corresponding to the general index (rayidx, lineidx, pointidx) in the second order feautrier form
/// Can use the same as above if we make the agreement that one uses the first nrays/2 for the u indices and the last nrays/2 for the v indices
/// However it is more convenient to explicitly ask whether we are using u or v
///  Returns the matrix index corresponding to the general index (rayidx (nrays/2), lineidx, pointidx, u_or_v)
inline Size Collocation :: get_mat_col_index_2nd_feautrier(Size basis_ray_id, Size basis_freq_id, Size basis_point_id, bool using_v)
{
    // return rayidx*parameters.nlines()*point_locations.size()+lineidx*point_locations.size()+pointidx;
    // return rayidx*parameters.nfreqs()*point_locations.size()+freqidx*point_locations.size()+pointidx;
    if (basis_ray_id>parameters.nrays()/2)
    {
      std::cout<<"this index only uses up to nrays/2; you have made a programming error";
    }
    //implicit bool to int conversion
    // return using_v*(parameters.nrays()/2)*parameters.nfreqs()*point_locations.size()+
    //        rayidx*parameters.nfreqs()*point_locations.size()+pointidx*parameters.nfreqs()+freqidx;
    return using_v*(parameters.nrays()/2)*cum_n_point_bases[n_freq_bases]+
           basis_ray_id*cum_n_point_bases[n_freq_bases]+cum_n_point_bases[basis_freq_id]+basis_point_id;
    //TODO: replace nfreqs() with number frequency quadratures
}

///  Returns the matrix index corresponding to the general index (rayidx, lineidx, pointidx) in the second order feautrier form
inline Size Collocation :: get_mat_row_index_2nd_feautrier(Size rayid, Size freqid, Size pointid, bool is_v_eq)
{
    // return rayidx*parameters.nlines()*point_locations.size()+lineidx*point_locations.size()+pointidx;
    // return rayidx*parameters.nfreqs()*point_locations.size()+freqidx*point_locations.size()+pointidx;
    if (rayid>parameters.nrays()/2)
    {
        std::cout<<"this index only uses up to nrays/2; you have made a programming error";
    }
    //implicit bool to int conversion
    // return is_v_eq*(parameters.nrays()/2)*parameters.nfreqs()*point_locations.size()+
    //        rayidx*parameters.nfreqs()*point_locations.size()+pointidx*parameters.nfreqs()+freqidx;
    return is_v_eq*(parameters.nrays()/2)*n_freq_bases*n_points_in_grid+
           rayid*n_freq_bases*n_points_in_grid+freqid*n_points_in_grid+pointid;
    //TODO: replace nfreqs() with number frequency quadratures
}

// inline Real Collocation :: calculate_doppler_shift(Size pointidx, Size rayidx, Geometry& geometry)
// {
//     Vector3D raydirection=model.geometry.rays.direction[rayindex];
//     return 1+std::sqrt(raydirection.dot())
// }

//DEPRECATED
//Computes the balancing factor needed in order to make the boundary condition equation the same magnitude as the other nearby equations in the matrix
//maybe TODO: also take the scattering stuff into account (when implemented)
// inline Real Collocation :: compute_balancing_factor_boundary(Size rayidx, Size freqidx, Size pointidx, Real local_velocity_gradient, Model& model)
// {
//     return 1;
//     // const Real numerator=basis_point_der(pointidx, point_locations[pointidx], rayidx, freqidx, model.geometry)*basis_freq(rayidx, freqidx, pointidx, non_doppler_shifted_frequencies[freqidx])-local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]*basis_point(pointidx, point_locations[pointidx], rayidx, freqidx, model.geometry)*basis_freq_der(rayidx, freqidx, pointidx, non_doppler_shifted_frequencies[freqidx]);
//     // const Real denominator=basis_point(pointidx, point_locations[pointidx], rayidx, freqidx, model.geometry)*get_opacity(model, rayidx, freqidx, pointidx)*basis_freq(rayidx, freqidx, pointidx, non_doppler_shifted_frequencies[freqidx]);
//     // // std::cout<<"numerator/denominator: "<<numerator/denominator<<std::endl;
//     //
//     // return 1 + numerator/denominator;
// }

///  Returns whether the boundary condition of the raditive transfer equation is need for a specific direction on a boundary point
///     @param[in] rayidx: the id of the ray direction on which we are currently looking
///     @param[in] point_on_boundary: id of point in original grid
inline bool Collocation :: is_boundary_condition_needed_for_direction_on_boundary(Model& model, Size rayidx, Size point_on_boundary)
{
    //For all neighboring boundary points, vote whether the direction is appropriate for the boundary condition
    bool has_boundary_neighbor=false;
    Size same_direction_votes=0;//if same_direction_votes*2>total_votes, then it is appropriate
    Size total_votes=0;

    const Vector3D curr_position=model.geometry.points.position[point_on_boundary];
    const Vector3D raydirection=model.geometry.rays.direction[rayidx];

    std::tuple<Size*,Size> temp_tuple=model.geometry.points.multiscale.get_intern_neighbors(point_on_boundary);
    Size* start_neighbors=std::get<0>(temp_tuple);
    Size n_neighbors=std::get<1>(temp_tuple);
    // for (Size n:temp_neighbors)
    for (Size i = 0; i < n_neighbors; i++)
    {
        const Size n=*(start_neighbors+i);
        if (!model.geometry.not_on_boundary(n))
        {
            has_boundary_neighbor=true;
            const Vector3D relative_diff_vector=model.geometry.points.position[n]-curr_position;//this needs not be normalized, as we only check whether inner product is positive/negative
            if (raydirection.dot(relative_diff_vector)>0.0)
            {
                same_direction_votes+=1;
            }
            total_votes+=1;
        }

    }
    std::cout<<"total_votes: "<<total_votes<<std::endl;
    std::cout<<"same_direction_votes: "<<same_direction_votes<<std::endl;
    //If no neighboring boundary points (e.g. in 1D), instead just check whether some point exists behind the boundary point
    if (!has_boundary_neighbor)
    {
        double dist=0;//these two variables are actually not used, but required in the function definition get_next
        double ddist=0;

        Size antipod=model.geometry.rays.antipod[rayidx];

        const Size point_behind=model.geometry.get_next(point_on_boundary, antipod, point_on_boundary, dist, ddist);
        std::cout<<"point_behind: "<<point_behind<<std::endl;
        //if there lies no point behind this point, we need to evaluate the boundary condition here
        return (point_behind==parameters.npoints());
    }
    else
    {
        return ((2*same_direction_votes)>=total_votes);//if at least half of the votes agree, we may put the boundary condition in this direction (in case of a tie, it is ambiguous whether we may use the boundary condition or not)
    }
}

/// Sets up the basis matrix for the default radiative transfer equation
inline void Collocation :: setup_basis_matrix_Eigen(Model& model)
{
    //We are currently not using the lambda operator, so set to zero
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    // Size mat_size=parameters.nrays()*parameters.nfreqs()*point_locations.size();
    const Size evaluation_size=parameters.nrays()*n_freq_bases*n_points_in_grid;//number of points at which we evaluate=is number of rows
    const Size tot_n_basis_functions=parameters.nrays()*cum_n_point_bases[n_freq_bases];//total number of basis functions=number of columns//TODO: define somewhere in general
    vector<Triplet<Real, Size>> eigen_triplets;
    // eigen_triplets.reserve(mat_size*16*nrays());//matrix size times number of interacting points times the number of rays (scattering)
    eigen_triplets.reserve(evaluation_size*N_POINTS_IN_RADIAL_BASIS);//matrix size times number of interacting points (does not include scattering) and frequencies
    //TODO: also multiply by estimate of average number/upper bound of interacting frequencies; currently we will need to reallocate

    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        std::cout<<"rayid: "<<rayid<<std::endl;
        std::cout<<"raydirection: "<<model.geometry.rays.direction[rayid].x()<<std::endl;
        // for (Size freqidx=0; freqidx<parameters.nfreqs(); freqidx++)
        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {

            for (Size pointid=0; pointid<n_points_in_grid; pointid++)
            {
                const Size point_in_grid=points_in_grid[pointid];
                //Now creating all the triplets, thus we need all nonzero basis triplets
                std::set<std::vector<Size>> nonzero_basis_triplets;
                get_nonzero_basis_triplets(nonzero_basis_triplets, rayid, freqid, pointid);

                // const Size curr_mat_idx=get_mat_index(rayidx, freqidx, pointidx);
                const Size mat_row_idx=get_mat_row_index(rayid, freqid, pointid);
                const Real curr_opacity=get_opacity(model, rayid, freqid, pointid);
                Vector3D curr_location=point_locations[pointid];

                //test output
                // const Real radius=rbf_radii[pointidx];
                // std::cout<<"1/(R*chi)"<<1.0/(radius*curr_opacity)<<std::endl;
                // std::cout<<"curr_opacity: "<<curr_opacity<<std::endl;
                // TODO: think about whether the frequency should be evaluated in the frame of the neighbors...//no, as we evaluate it at a single point each time (no interpolation stuff required)

                //TODO: maybe refactor this part
                Real local_velocity_gradient=0;

                double dist=0;
                double ddist=0;
                Size next_point=model.geometry.get_next(point_in_grid, rayid, point_in_grid, dist, ddist);
                if (next_point!=parameters.npoints())//if there actually exists a next point
                {
                    // std::cout<<"dist: "<<dist<<std::endl;
                    local_velocity_gradient=model.geometry.rays.direction[rayid].dot(model.geometry.points.velocity[next_point]-model.geometry.points.velocity[point_in_grid])/dist;
                    // std::cout<<"local velocity grad: "<<local_velocity_gradient<<std::endl;
                }
                else
                {
                    double dist=0;
                    double ddist=0;
                    Size prev_point=model.geometry.get_next(point_in_grid, model.geometry.rays.antipod[rayid], point_in_grid, dist, ddist);
                    if (prev_point!=parameters.npoints())
                    {
                        local_velocity_gradient=model.geometry.rays.direction[rayid].dot(model.geometry.points.velocity[point_in_grid]-model.geometry.points.velocity[prev_point])/dist;
                    }//otherwise we only have a single point on the ray in this direction; so defining a non-zero velocity gradient seems silly
                }
                //TODO: refactor until here


                Real curr_freq=0;
                if (USING_COMOVING_FRAME)
                {
                    curr_freq=non_doppler_shifted_frequencies[freqid];
                }
                else
                {
                    curr_freq=frequencies[rayid][pointid][freqid];
                }

                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(point_in_grid))
                {
                    const bool bdy_condition_req_in_direction=is_boundary_condition_needed_for_direction_on_boundary(model, rayid, point_in_grid);
                    boundary_condition_required=bdy_condition_req_in_direction;
                }

                if (boundary_condition_required)
                {
                    //Merely evaluating the intensity at the boundary
                    for (std::vector<Size> triplet: nonzero_basis_triplets)
                    {
                        const Size mat_col_idx=get_mat_col_index(triplet[0],triplet[1],triplet[2]);
                        // std::cout<<"mat col index: "<<mat_col_idx<<std::endl;
                        // std::cout<<"triplet indices: "<<triplet[0]<<","<<triplet[1]<<","<<triplet[2]<<std::endl;

                        if (USING_COMOVING_FRAME)
                        {
                            // const Real doppler_shifted_frequency=curr_freq*calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx]
                            //                                                                                       , model.geometry.rays.direction[rayidx]);
                            // Real doppler_shifted_frequency=curr_freq;//try it out non-doppler shifted
                            // Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point(triplet[2], curr_location, triplet[0], model.geometry);
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, basis_eval));
                        }
                        else
                        {
                            //In order to make this term of the same size as the others, multiply by the opacity (only in emissivity formulation)
                            // Real basis_eval=curr_opacity*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            Real boundary_condition_mult_factor=boundary_condition_LS_importance_factor;
                            // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                            // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, basis_eval*boundary_condition_mult_factor));
                        }
                        std::cout<<"bdy condition triplets: "<<triplet[0]<<", "<<triplet[1]<<", "<<triplet[2]<<std::endl;
                        std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                        std::cout<<"basis_point: "<<basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)<<std::endl;

                    }
                }
                else
                {
                    for (std::vector<Size> triplet: nonzero_basis_triplets)
                    {
                        const Size mat_col_idx=get_mat_col_index(triplet[0],triplet[1],triplet[2]);
                        // std::cout<<"mat col index: "<<mat_col_idx<<std::endl;
                        // std::cout<<"triplet indices: "<<triplet[0]<<","<<triplet[1]<<","<<triplet[2]<<std::endl;

                        if (USING_COMOVING_FRAME)
                        {   //TODO: give a more solid reasoning why the doppler shift should be calculated using the rayidx, instead of triplet[0]
                            // const Real doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx]
                            //                                                                 , model.geometry.rays.direction[rayidx]);
                            // const Real doppler_shifted_frequency=curr_freq*doppler_shift;
                            // const Real doppler_shifted_opacity=get_opacity(model, rayidx, lineidx, pointidx, doppler_shift);
                            // Real doppler_shifted_frequency=curr_freq;//try it out non-doppler shifted
                            // Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry)/curr_opacity;///doppler_shifted_opacity;///curr_opacity;
                            // Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point(triplet[2], curr_location, triplet[0], model.geometry);
                            Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)/curr_opacity;///doppler_shifted_opacity;///curr_opacity;
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            // std::cout<<"curr_opacity: "<<curr_opacity<<std::endl;
                            // std::cout<<"basis_eval+dir_der: "<<basis_eval+directional_der<<std::endl;
                            // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"dir_der: "<<directional_der<<std::endl;
                            // std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_point: "<<basis_point(triplet[2], curr_location, triplet[0], model.geometry)<<std::endl;
                            Real doppler_shift_term=0;//TODO: add frequency quadrature;; TODO: think about approximating velocity gradient somehow at the point itself?

                            doppler_shift_term=-local_velocity_gradient*curr_freq/curr_opacity*
                                  basis_direction(triplet[0])*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);

                            // std::cout<<"velocity grad: "<<local_velocity_gradient<<std::endl;
                            // std::cout<<"freq: "<<non_doppler_shifted_frequencies[freqidx]<<std::endl;
                            // std::cout<<"1/opacity: "<<1.0/curr_opacity<<std::endl;
                            // std::cout<<"doppler_shift_spat_der_ratio"<<local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)/basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], model.geometry)/basis_point_der(triplet[2], curr_location, triplet[0], model.geometry)<<std::endl;
                            // std::cout<<"inverse width: "<<non_doppler_shifted_inverse_widths[triplet[1]][triplet[2]]<<std::endl;
                            // std::cout<<"doppler_shift_term: "<<doppler_shift_term<<std::endl;
                            // std::cout<<"doppler_shift_basis_ratio: "<<-local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]/curr_opacity*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)/basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // doppler_shift_term=0;

                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, directional_der+basis_eval+doppler_shift_term));

                            // //For every scattering direction! add scattering//also divide by curr_opacity
                            // for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                            // {
                            //     const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                            //     //TODO: add toggle to whether we are actually using scattering
                            //     Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                            //     Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                            //     //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                            //     // eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                            // }
                        }
                        else
                        {
                            //Start with the directional derivative (in source function formulation divided by the opacity)
                            // Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry);
                            Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)/curr_opacity;
                            // std::cout<<"directional_der: "<<directional_der<<std::endl;
                            // directional_der=0;
                            //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                            // Real opacity=curr_opacity*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                            // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_point_eval: "<<basis_point(triplet[2], curr_location)<<std::endl;
                            // basis_eval*=0.01;

                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, directional_der+basis_eval));

                            // //For every scattering direction! add scattering//also divide by curr_opacity
                            // for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                            // {
                            //     const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                            //     //TODO: add toggle to whether we are actually using scattering
                            //     Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                            //     Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                            //     //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                            //     // eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                            // }

                        }

                    }
                }
            }
        }
    }

    // std::cout<<"resizing collocation mat"<<std::endl;

    //after all this, construct the matrix
    collocation_mat.resize (evaluation_size, tot_n_basis_functions);

    // std::cout<<"setting collocation mat"<<std::endl;
    collocation_mat.setFromTriplets (eigen_triplets.begin(), eigen_triplets.end());

    //TODO: add toggle for computing svd and condition number
    if (model.COMPUTING_SVD_AND_CONDITION_NUMBER)
    {
        Eigen::JacobiSVD<MatrixXr> svd(collocation_mat);
        Real cond = svd.singularValues()(0)
        / svd.singularValues()(svd.singularValues().size()-1);
        std::cout<<"condition number: "<<cond<<std::endl;
        // collocation_mat.data().squeeze();
        std::cout << MatrixXr(collocation_mat) << std::endl;
    }
}


/// Sets up the basis matrix for the default radiative transfer equation
/// Using only the evaluation points at all basis locations
inline void Collocation :: setup_basis_square_matrix_Eigen(Model& model)
{
    //We are currently not using the lambda operator, so set to zero
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    // Size mat_size=parameters.nrays()*parameters.nfreqs()*point_locations.size();
    // const Size evaluation_size=parameters.nrays()*n_freq_bases*n_points_in_grid;//number of points at which we evaluate=is number of rows
    const Size tot_n_basis_functions=parameters.nrays()*cum_n_point_bases[n_freq_bases];//total number of basis functions=number of columns//TODO: define somewhere in general
    const Size evaluation_size=tot_n_basis_functions;
    vector<Triplet<Real, Size>> eigen_triplets;
    eigen_triplets.reserve(evaluation_size*N_POINTS_IN_RADIAL_BASIS);//matrix size times number of interacting points (does not include scattering) and frequencies
    //TODO: also multiply by estimate of average number/upper bound of interacting frequencies; currently we will need to reallocate

    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        std::cout<<"rayid: "<<rayid<<std::endl;
        std::cout<<"raydirection: "<<model.geometry.rays.direction[rayid].x()<<std::endl;
        // for (Size freqidx=0; freqidx<parameters.nfreqs(); freqidx++)
        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {

            for (Size point_basis_id=0; point_basis_id<n_point_bases[freqid]; point_basis_id++)
            {
                const Size pointid=local_basis_point_idx[freqid][point_basis_id];
                const Size point_in_grid=basis_point_idx[freqid][point_basis_id];
                //Now creating all the triplets, thus we need all nonzero basis triplets
                std::set<std::vector<Size>> nonzero_basis_triplets;
                get_nonzero_basis_triplets(nonzero_basis_triplets, rayid, freqid, pointid);

                // const Size curr_mat_idx=get_mat_index(rayidx, freqidx, pointidx);
                // const Size mat_row_idx=get_mat_row_index(rayid, freqid, pointid);
                const Size mat_row_idx=get_mat_col_index(rayid, freqid, point_basis_id);//As we want to make a square matrix, we replace the row indices with column indices
                const Real curr_opacity=get_opacity(model, rayid, freqid, pointid);
                Vector3D curr_location=point_locations[pointid];

                //TODO: maybe refactor this part
                Real local_velocity_gradient=0;

                double dist=0;
                double ddist=0;
                Size next_point=model.geometry.get_next(point_in_grid, rayid, point_in_grid, dist, ddist);
                if (next_point!=parameters.npoints())//if there actually exists a next point
                {
                    // std::cout<<"dist: "<<dist<<std::endl;
                    local_velocity_gradient=model.geometry.rays.direction[rayid].dot(model.geometry.points.velocity[next_point]-model.geometry.points.velocity[point_in_grid])/dist;
                    // std::cout<<"local velocity grad: "<<local_velocity_gradient<<std::endl;
                }
                else
                {
                    double dist=0;
                    double ddist=0;
                    Size prev_point=model.geometry.get_next(point_in_grid, model.geometry.rays.antipod[rayid], point_in_grid, dist, ddist);
                    if (prev_point!=parameters.npoints())
                    {
                        local_velocity_gradient=model.geometry.rays.direction[rayid].dot(model.geometry.points.velocity[point_in_grid]-model.geometry.points.velocity[prev_point])/dist;
                    }//otherwise we only have a single point on the ray in this direction; so defining a non-zero velocity gradient seems silly
                }
                //TODO: refactor until here


                Real curr_freq=0;
                if (USING_COMOVING_FRAME)
                {
                    curr_freq=non_doppler_shifted_frequencies[freqid];
                }
                else
                {
                    curr_freq=frequencies[rayid][pointid][freqid];
                }

                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(point_in_grid))
                {
                    const bool bdy_condition_req_in_direction=is_boundary_condition_needed_for_direction_on_boundary(model, rayid, point_in_grid);
                    boundary_condition_required=bdy_condition_req_in_direction;
                }

                if (boundary_condition_required)
                {
                    //Merely evaluating the intensity at the boundary
                    for (std::vector<Size> triplet: nonzero_basis_triplets)
                    {
                        const Size mat_col_idx=get_mat_col_index(triplet[0],triplet[1],triplet[2]);
                        // std::cout<<"mat col index: "<<mat_col_idx<<std::endl;
                        // std::cout<<"triplet indices: "<<triplet[0]<<","<<triplet[1]<<","<<triplet[2]<<std::endl;

                        if (USING_COMOVING_FRAME)
                        {
                            // const Real doppler_shifted_frequency=curr_freq*calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx]
                            //                                                                                       , model.geometry.rays.direction[rayidx]);
                            // Real doppler_shifted_frequency=curr_freq;//try it out non-doppler shifted
                            // Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point(triplet[2], curr_location, triplet[0], model.geometry);
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, basis_eval));
                        }
                        else
                        {
                            //In order to make this term of the same size as the others, multiply by the opacity (only in emissivity formulation)
                            // Real basis_eval=curr_opacity*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                            // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, basis_eval));
                        }
                        std::cout<<"bdy condition triplets: "<<triplet[0]<<", "<<triplet[1]<<", "<<triplet[2]<<std::endl;
                        std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                        std::cout<<"basis_point: "<<basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)<<std::endl;

                    }
                }
                else
                {
                    for (std::vector<Size> triplet: nonzero_basis_triplets)
                    {
                        const Size mat_col_idx=get_mat_col_index(triplet[0],triplet[1],triplet[2]);
                        // std::cout<<"mat col index: "<<mat_col_idx<<std::endl;
                        // std::cout<<"triplet indices: "<<triplet[0]<<","<<triplet[1]<<","<<triplet[2]<<std::endl;

                        if (USING_COMOVING_FRAME)
                        {   //TODO: give a more solid reasoning why the doppler shift should be calculated using the rayidx, instead of triplet[0]
                            // const Real doppler_shift=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[pointidx]
                            //                                                                 , model.geometry.rays.direction[rayidx]);
                            // const Real doppler_shifted_frequency=curr_freq*doppler_shift;
                            // const Real doppler_shifted_opacity=get_opacity(model, rayidx, lineidx, pointidx, doppler_shift);
                            // Real doppler_shifted_frequency=curr_freq;//try it out non-doppler shifted
                            // Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry)/curr_opacity;///doppler_shifted_opacity;///curr_opacity;
                            // Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency)*basis_point(triplet[2], curr_location, triplet[0], model.geometry);
                            Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)/curr_opacity;///doppler_shifted_opacity;///curr_opacity;
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            // std::cout<<"curr_opacity: "<<curr_opacity<<std::endl;
                            // std::cout<<"basis_eval+dir_der: "<<basis_eval+directional_der<<std::endl;
                            // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"dir_der: "<<directional_der<<std::endl;
                            // std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_point: "<<basis_point(triplet[2], curr_location, triplet[0], model.geometry)<<std::endl;
                            Real doppler_shift_term=0;//TODO: add frequency quadrature;; TODO: think about approximating velocity gradient somehow at the point itself?

                            doppler_shift_term=-local_velocity_gradient*curr_freq/curr_opacity*
                                  basis_direction(triplet[0])*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);

                            // std::cout<<"velocity grad: "<<local_velocity_gradient<<std::endl;
                            // std::cout<<"freq: "<<non_doppler_shifted_frequencies[freqidx]<<std::endl;
                            // std::cout<<"1/opacity: "<<1.0/curr_opacity<<std::endl;
                            // std::cout<<"doppler_shift_spat_der_ratio"<<local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)/basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location, triplet[0], model.geometry)/basis_point_der(triplet[2], curr_location, triplet[0], model.geometry)<<std::endl;
                            // std::cout<<"inverse width: "<<non_doppler_shifted_inverse_widths[triplet[1]][triplet[2]]<<std::endl;
                            // std::cout<<"doppler_shift_term: "<<doppler_shift_term<<std::endl;
                            // std::cout<<"doppler_shift_basis_ratio: "<<-local_velocity_gradient*non_doppler_shifted_frequencies[freqidx]/curr_opacity*basis_freq_der(triplet[0], triplet[1], triplet[2], curr_freq)/basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // doppler_shift_term=0;

                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, directional_der+basis_eval+doppler_shift_term));

                            // //For every scattering direction! add scattering//also divide by curr_opacity
                            // for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                            // {
                            //     const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                            //     //TODO: add toggle to whether we are actually using scattering
                            //     Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                            //     Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                            //     //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                            //     // eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                            // }
                        }
                        else
                        {
                            //Start with the directional derivative (in source function formulation divided by the opacity)
                            // Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], model.geometry);
                            Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)/curr_opacity;
                            // std::cout<<"directional_der: "<<directional_der<<std::endl;
                            // directional_der=0;
                            //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                            // Real opacity=curr_opacity*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                            Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                            // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                            // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                            // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                            // std::cout<<"basis_point_eval: "<<basis_point(triplet[2], curr_location)<<std::endl;
                            // basis_eval*=0.01;

                            eigen_triplets.push_back (Triplet<Real, Size> (mat_row_idx, mat_col_idx, directional_der+basis_eval));

                            // //For every scattering direction! add scattering//also divide by curr_opacity
                            // for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                            // {
                            //     const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                            //     //TODO: add toggle to whether we are actually using scattering
                            //     Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                            //     Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                            //     //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                            //     // eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                            // }

                        }

                    }
                }
            }
        }
    }

    // std::cout<<"resizing collocation mat"<<std::endl;

    //after all this, construct the matrix
    collocation_mat.resize (evaluation_size, tot_n_basis_functions);

    // std::cout<<"setting collocation mat"<<std::endl;
    collocation_mat.setFromTriplets (eigen_triplets.begin(), eigen_triplets.end());

    //TODO: add toggle for computing svd and condition number
    if (model.COMPUTING_SVD_AND_CONDITION_NUMBER)
    {
        Eigen::JacobiSVD<MatrixXr> svd(collocation_mat);
        Real cond = svd.singularValues()(0)
        / svd.singularValues()(svd.singularValues().size()-1);
        std::cout<<"condition number: "<<cond<<std::endl;
        // collocation_mat.data().squeeze();
        std::cout << MatrixXr(collocation_mat) << std::endl;
    }
}




/// Sets up the basis matrix for the second order feautrier version of the radiative transfer equation
/// Only defined for the static frame
inline void Collocation :: setup_basis_matrix_Eigen_2nd_feautrier(Model& model)
{

    if (USING_COMOVING_FRAME)//In the comoving frame, the equations will get way too complicated
    {
        throw std::runtime_error("Second order feautrier collocation method only defined for the non-comoving frame");
    }

    //We are currently not using the lambda operator, so set to zero
    for (auto &lspec : model.lines.lineProducingSpecies) {lspec.lambda.clear();}

    const Size evaluation_size=parameters.nrays()*n_freq_bases*n_points_in_grid;//number of points at which we evaluate=is number of rows
    const Size tot_n_basis_functions=parameters.nrays()*cum_n_point_bases[n_freq_bases];//total number of basis functions=number of columns//TODO: define somewhere in general

    vector<Triplet<Real, Size>> eigen_triplets;//will first contain all triplets corresponding to the u equations (and the boundary conditions). At the end, the v_eigen_triplets will be apprended
    vector<Triplet<Real, Size>> v_eigen_triplets;//Temporary container for all triplets corresponding to the v equations.
    // eigen_triplets.reserve(mat_size*16*nrays());//matrix size times number of interacting points times the number of rays (scattering)
    eigen_triplets.reserve(evaluation_size*N_POINTS_IN_RADIAL_BASIS);//matrix size times number of interacting points (does not include scattering) and frequencies
    v_eigen_triplets.reserve(evaluation_size*N_POINTS_IN_RADIAL_BASIS/2);//only stores half, so half the space needs to be reserved

    //TODO: in the Eigen docs is stated that the triplets vector does not necessarilly need to be ordered. So I could use a single vector instead
    // eigen_triplets.reserve(mat_size*16*nrays());//matrix size times number of interacting points times the number of rays (scattering)
    // eigen_triplets.reserve(evaluation_size*16);//matrix size times number of interacting points (does not include scattering)
    // v_eigen_triplets.reserve(evaluation_size*8);//only stores half, so half the space needs to be reserved

    for (Size rayid=0; rayid<parameters.hnrays(); rayid++)
    {
        std::cout<<"rayid: "<<rayid<<std::endl;
        // std::cout<<"raydirection: "<<model.geometry.rays.direction[rayidx].x()<<std::endl;

        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {

            for (Size pointid=0; pointid<n_points_in_grid; pointid++)
            {
                const Size point_in_grid=points_in_grid[pointid];
                //Now creating all the triplets, thus we need all nonzero basis triplets
                std::set<std::vector<Size>> nonzero_basis_triplets;
                get_nonzero_basis_triplets(nonzero_basis_triplets, rayid, freqid, pointid);

                // const Size curr_mat_idx=get_mat_index(rayidx, freqidx, pointidx);
                // const Size mat_row_idx=get_mat_row_index(rayidx, freqidx, pointidx);
                const Size u_eq_mat_row_idx=get_mat_row_index_2nd_feautrier(rayid, freqid, pointid, false);
                const Size v_eq_mat_row_idx=get_mat_row_index_2nd_feautrier(rayid, freqid, pointid, true);
                const Real curr_opacity=get_opacity(model, rayid, freqid, pointid);
                std::cout<<"curr_opacity: "<<curr_opacity<<std::endl;

                const Real rel_opacity_grad=get_opacity_grad(model, rayid, freqid, pointid)/curr_opacity;//relative opacity gradient
                std::cout<<"rel_opacity_grad: "<<rel_opacity_grad<<std::endl;
                const Real curr_freq=frequencies[rayid][pointid][freqid];
                Vector3D curr_location=point_locations[pointid];

                bool boundary_condition_required=false;
                bool bdy_condition_req_in_direction=false;
                if (!model.geometry.not_on_boundary(point_in_grid))
                {
                    bdy_condition_req_in_direction=is_boundary_condition_needed_for_direction_on_boundary(model, rayid, point_in_grid);
                    boundary_condition_required=true;
                    //on all boundary points, we need a boundary condition (far simpler than the default approach)
                    //But we still need to determine in which direction the boundary condition lies
                }

                if (boundary_condition_required)
                {
                    //Merely evaluating the intensity at the boundary in u (+-v) (replacing 2nd order equation for u) and using the 2nd order equation for v
                    for (std::vector<Size> triplet: nonzero_basis_triplets)
                    {
                        // const Size mat_col_idx=get_mat_col_index(triplet[0],triplet[1],triplet[2]);
                        const Size u_mat_col_idx=get_mat_col_index_2nd_feautrier(triplet[0], triplet[1], triplet[2], false);
                        const Size v_mat_col_idx=get_mat_col_index_2nd_feautrier(triplet[0], triplet[1], triplet[2], true);

                        //In order to make this term of the same size as the others, multiply by the opacity (only in emissivity formulation)
                        // Real basis_eval=curr_opacity*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location);
                        Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                        // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                        // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                        // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                        //add triplet (mat_idx(triplet[0],triplet[1],triplet[2]),directional_der)
                        eigen_triplets.push_back (Triplet<Real, Size> (u_eq_mat_row_idx, u_mat_col_idx, basis_eval));

                        Real v_mult_factor=-1.0;
                        std::cout<<"bdy req in direction nhat: "<<bdy_condition_req_in_direction<<std::endl;
                        if (bdy_condition_req_in_direction)//v changes sign on the different boundary conditions in different directions
                        {
                            v_mult_factor=1.0;
                        }

                        eigen_triplets.push_back (Triplet<Real, Size> (u_eq_mat_row_idx, v_mat_col_idx, v_mult_factor*basis_eval));

                        std::cout<<"bdy condition triplets: "<<triplet[0]<<", "<<triplet[1]<<", "<<triplet[2]<<std::endl;
                        std::cout<<"basis_freq: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                        std::cout<<"basis_point: "<<basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)<<std::endl;


                        //And also compute the standard second order feautrier term for the v equations
                        //As we use the same basis functions for u and v, we do not need to make a distinction between them
                        //TODO: check if previous sentence is true

                        // Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[2], curr_location, triplet[0], triplet[1], model.geometry)
                        // *rel_opacity_grad/std::pow(curr_opacity,2);
                        // Real directional_der2=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der2(triplet[2], curr_location, triplet[0], triplet[1], model.geometry)
                        // /std::pow(curr_opacity,2);

                        //only simple first order term necessary on boundary
                        Real simple_directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)
                        /curr_opacity;
                        //APPROXIMATED BOUNDARY CONDITION
                        // simple_directional_der=0;

                        // v_eigen_triplets.push_back (Triplet<Real, Size> (v_eq_mat_row_idx, v_mat_col_idx, basis_eval-directional_der2+directional_der));
                        v_eigen_triplets.push_back (Triplet<Real, Size> (v_eq_mat_row_idx, v_mat_col_idx, basis_eval-v_mult_factor*simple_directional_der));
                    }
                }
                else
                {
                    for (std::vector<Size> triplet: nonzero_basis_triplets)
                    {
                        // const Size mat_col_idx=get_mat_col_index(triplet[0],triplet[1],triplet[2]);
                        const Size u_mat_col_idx=get_mat_col_index_2nd_feautrier(triplet[0], triplet[1], triplet[2], false);
                        const Size v_mat_col_idx=get_mat_col_index_2nd_feautrier(triplet[0], triplet[1], triplet[2], true);

                        Real basis_eval=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point(triplet[0], triplet[1], triplet[2], curr_location, model.geometry);
                        Real directional_der=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)
                        *rel_opacity_grad/std::pow(curr_opacity,2);
                        Real directional_der2=basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)*basis_point_der2(triplet[0], triplet[1], triplet[2], curr_location, model.geometry)
                        /std::pow(curr_opacity,2);
                        std::cout<<"dir_der2: "<<directional_der2<<std::endl;
                        // basis_eval=0;
                        // directional_der=0;
                        // directional_der2*=10;

                        // std::cout<<"basis_eval: "<<basis_eval<<std::endl;
                        // std::cout<<"basis_dir_eval: "<<basis_direction(triplet[0])<<std::endl;
                        // std::cout<<"basis_freq_eval: "<<basis_freq(triplet[0], triplet[1], triplet[2], curr_freq)<<std::endl;
                        // std::cout<<"basis_point_eval: "<<basis_point(triplet[2], curr_location)<<std::endl;

                        // std::cout<<"basis/(grad/opacity): "<<basis_eval/directional_der<<std::endl;

                        eigen_triplets.push_back (Triplet<Real, Size> (u_eq_mat_row_idx, u_mat_col_idx, basis_eval-directional_der2+directional_der));
                        v_eigen_triplets.push_back (Triplet<Real, Size> (v_eq_mat_row_idx, v_mat_col_idx, basis_eval-directional_der2+directional_der));

                        // //For every scattering direction! add scattering//also divide by curr_opacity
                        // for (Size temp_rayidx=0; temp_rayidx<parameters.nrays(); temp_rayidx++)
                        // {
                        //     const Size scatt_mat_idx=get_mat_index(temp_rayidx,triplet[1],triplet[2]);
                        //     //TODO: add toggle to whether we are actually using scattering
                        //     Real integral_scatt_redistr=0;//FIXME: add me //in uniform scattering, this term is proportional to 4*pi/N_rays
                        //     Real scattering=-basis_freq(temp_rayidx, triplet[1], triplet[2], curr_freq)*basis_point(triplet[2], curr_location)*integral_scatt_redistr;
                        //     //add triplet (mat_idx(temp_rayidx,triplet[1],triplet[2]),scattering)
                        //     // eigen_triplets.push_back (Triplet<Real, Size> (curr_mat_idx, scatt_mat_idx, scattering));
                        // }
                    }
                }
            }
        }
    }
    // and finally adding the triplets together
    eigen_triplets.insert(eigen_triplets.end(), v_eigen_triplets.begin(), v_eigen_triplets.end());

    // std::cout<<"resizing collocation mat"<<std::endl;

    //after all this, construct the matrix
    collocation_mat.resize (evaluation_size, tot_n_basis_functions);

    // std::cout<<"setting collocation mat"<<std::endl;
    collocation_mat.setFromTriplets (eigen_triplets.begin(), eigen_triplets.end());

    //TODO: add toggle for computing svd and condition number
    if (model.COMPUTING_SVD_AND_CONDITION_NUMBER)
    {
        Eigen::JacobiSVD<MatrixXr> svd(collocation_mat);
        Real cond = svd.singularValues()(0)
        / svd.singularValues()(svd.singularValues().size()-1);
        std::cout<<"condition number: "<<cond<<std::endl;
        // collocation_mat.data().squeeze();
        std::cout << MatrixXr(collocation_mat) << std::endl;
    }
}







// Returns the opacity (takes into account the gaussian line profile)
inline Real Collocation :: get_opacity(Model& model, Size rayid, Size freqid, Size pointid)
{
    Real toreturn=0;
    const Size point_in_grid=points_in_grid[pointid];
    // if using other line profile/frequency basis function than gaussian, please replace the for(...) with the following line:
    // for (Size other_line_idx=0; other_line_idx<parameters.nlines(); other_line_idx++)
    // for (Size other_line_idx: other_frequency_basis_interacting_with[pointidx][lineidx])
    for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
    {
        Real delta_freq=0;
        // Real inv_width=0;
        if (USING_COMOVING_FRAME)
        {
            // inv_width=model.lines.inverse_width(pointidx,lineidx);
            delta_freq=non_doppler_shifted_frequencies[freqid]-model.lines.line[lineidx];
        }
        else
        {
            // Real doppler_shift_factor=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointidx],model.geometry.rays.direction[rayidx]);
            //local opacity should use the local frequencies, not the global ones...; otherwise it goes very quickly to zero if we have a minor doppler shift during the ray
            // I just don't feel like doppler shifting everything -> only minor correction in the inv_width factor in front is countered by the doppler shifted line frequency
            delta_freq=non_doppler_shifted_frequencies[freqid]-model.lines.line[lineidx];
            // delta_freq=frequencies[rayidx][pointidx][freqidx]-model.lines.line[lineidx]*doppler_shift_factor;
            // inv_width=model.lines.inverse_width(pointidx,lineidx)/doppler_shift_factor;
        }
        const Real inv_width=model.lines.inverse_width(point_in_grid,lineidx);
        toreturn+=INVERSE_SQRT_PI*inv_width*std::exp(-std::pow(delta_freq*inv_width,2))
                  *model.lines.line[lineidx]*model.lines.opacity(point_in_grid, lineidx);
    }
    //with a too low minimum opacity, the intensity gets completely ignored compared to the derivative
    return std::max(toreturn,MIN_OPACITY);
    // return toreturn;//disabling the minimal opacity, as it might wreck havoc with very off-line center frequencies; they also need to be treated correctly
}


// Returns the opacity (takes into account the gaussian line profile)
// The frequency should include the doppler shift
inline Real Collocation :: get_opacity_bound(Model& model, Size pointid, Size freqid, Size other_pointid)
{
    Real toreturn=0;
    const Size point_in_grid=points_in_grid[pointid];
    const Size other_point_in_grid=points_in_grid[other_pointid];
    // if using other line profile/frequency basis function than gaussian, please replace the for(...) with the following line:
    // for (Size other_line_idx=0; other_line_idx<parameters.nlines(); other_line_idx++)
    // for (Size other_line_idx: other_frequency_basis_interacting_with[pointidx][lineidx])
    //compute bound of opacity over all directions
    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        const Real doppler_shift_factor=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[other_point_in_grid],model.geometry.rays.direction[rayid]);
        Real temp_total_opacity=0;
        for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
        {
            Real delta_freq=0;
            Real inv_width=0;

            if (USING_COMOVING_FRAME)
            {
                delta_freq=non_doppler_shifted_frequencies[freqid]-model.lines.line[lineidx];
                inv_width=model.lines.inverse_width(point_in_grid,lineidx);
            }
            else
            {
                //local opacity should use the local frequencies, not the global ones...; otherwise it goes very quickly to zero if we have a minor doppler shift during the ray
                delta_freq=frequencies[rayid][freqid][pointid]-model.lines.line[lineidx]*doppler_shift_factor;
                inv_width=model.lines.inverse_width(other_point_in_grid,lineidx)/doppler_shift_factor;
            }

            temp_total_opacity+=INVERSE_SQRT_PI*inv_width*std::exp(-std::pow(delta_freq*inv_width,2))
                                *model.lines.line[lineidx]*model.lines.opacity(other_point_in_grid, lineidx);
            }
        toreturn=std::max(toreturn, temp_total_opacity);

    }
    //with a too low minimum opacity, the intensity gets completely ignored compared to the derivative
    return std::max(toreturn,MIN_OPACITY);
    // return toreturn;//disabling the minimal opacity, as it might wreck havoc with very off-line center frequencies; they also need to be treated correctly
}



//Returns a vector with two nearby points on the specified ray around the specified point
// In the fail case no two nearby points could be found, can return the oob values parameters.npoints().
inline std::vector<Size> Collocation :: get_two_nearby_points_on_ray(Model& model, Size rayidx, Size pointidx)
{
    std::vector<Size> toreturn;
    toreturn.reserve(2);
    double dist=0;
    double ddist=0;
    Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);

    dist=0;
    ddist=0;
    Size prev_point=model.geometry.get_next(pointidx, model.geometry.rays.antipod[rayidx], pointidx, dist, ddist);

    if (next_point==parameters.npoints()&&prev_point==parameters.npoints())//In the case we have no neighbors, return the default value parameters.npoints() for both points
    {
        std::cout<<"points on both sides do not exist; therefore returing"<<std::endl;
        toreturn=std::vector<Size>{next_point, prev_point};
        return toreturn;
    }

    if (next_point==parameters.npoints())//if next point in direction does not exist, just replace by another point in the other direction
    {
        std::cout<<"next_point==npoints"<<std::endl;
        dist=0;
        ddist=0;
        Size next_point=model.geometry.get_next(prev_point, model.geometry.rays.antipod[rayidx], prev_point, dist, ddist);
        //if this does not exist, it gains the value parameters.npoints()
    }
    else if (prev_point==parameters.npoints())//if next point in direction does not exist, just replace by another point in the other direction
    {
        std::cout<<"prev_point==npoints"<<std::endl;
        double dist=0;
        double ddist=0;
        Size prev_point=model.geometry.get_next(next_point, rayidx, next_point, dist, ddist);
        //if this does not exist, it gains the value parameters.npoints()
    }

    toreturn=std::vector<Size>{next_point, prev_point};
    return toreturn;

}

//Computes the opacity gradient along a ray direction at a given position by using a quadrature rule
// Uses an implicit total derivative (TODO: ref to own derivation)
//Warning: In the case we are on the boundary, it assumes we can find two points in a single direction
//Warning: assumes we use the same grid in both the original code and the collocation code (is valid, unless the original grid has changed after)
inline Real Collocation :: get_opacity_grad(Model& model, Size rayid, Size freqid, Size pointid)
{
    const Size point_in_grid=points_in_grid[pointid];
    // //TODO: refactor this part of finding two near points on the ray.
    //
    // double dist=0;
    // double ddist=0;
    // Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
    //
    // dist=0;
    // ddist=0;
    // Size prev_point=model.geometry.get_next(pointidx, model.geometry.rays.antipod[rayidx], pointidx, dist, ddist);
    //
    // if (next_point==parameters.npoints())//if next point in direction does not exist, just replace by another point in the other direction
    // {
    //     dist=0;
    //     ddist=0;
    //     Size next_point=model.geometry.get_next(pointidx, model.geometry.rays.antipod[rayidx], prev_point, dist, ddist);
    //     //ASSUMED to exist; naming becomes a bit non-sensical, but we just need to other points on the ray
    // }
    // else if (prev_point==parameters.npoints())//if next point in direction does not exist, just replace by another point in the other direction
    // {
    //     double dist=0;
    //     double ddist=0;
    //     Size prev_point=model.geometry.get_next(next_point, rayidx, pointidx, dist, ddist);
    //     //ASSUMED to exist; naming becomes a bit non-sensical, but we just need to other points on the ray
    // }
    vector<Size> nearby_points=get_two_nearby_points_on_ray(model, rayid, point_in_grid);
    Size next_point=nearby_points[0];
    Size prev_point=nearby_points[1];

    if (next_point==parameters.npoints()||prev_point==parameters.npoints())//if at least one of the nearby points do no exist, then no gradient can be computed
    {//I do not expect this to happen at all, but either way,
        std::cout<<"A ray through a point does not go through at least 3 points"<<std::endl;
        std::cout<<"pointid: "<<pointid<<std::endl;
        std::cout<<"rayid: "<<rayid<<std::endl;
        std::cout<<"next point: "<<next_point<<std::endl;
        std::cout<<"prev point: "<<prev_point<<std::endl;
        return 0;
    }

    Vector3D raydirection=model.geometry.rays.direction[rayid];
    Vector3D next_vector=model.geometry.points.position[next_point]-model.geometry.points.position[point_in_grid];
    Vector3D prev_vector=model.geometry.points.position[prev_point]-model.geometry.points.position[point_in_grid];
    const Real next_position=raydirection.dot(next_vector);
    const Real prev_position=raydirection.dot(prev_vector);

    const Real position_frac=next_position/prev_position;
    Real prev_coef=1.0/(prev_position*(1.0-(1.0/position_frac)));
    Real next_coef=1.0/(next_position*(1.0-position_frac));
    Real curr_coef=-prev_coef-next_coef;

    Real prev_opacity=get_opacity(model, rayid, freqid, prev_point);//<-Warning: same grid assumed
    Real curr_opacity=get_opacity(model, rayid, freqid, pointid);
    Real next_opacity=get_opacity(model, rayid, freqid, next_point);//<-Warning: same grid assumed
    //comoving opacities are used

    return prev_coef*prev_opacity
          +curr_coef*curr_opacity
          +next_coef*next_opacity;
}

//Computes the gradient of the source function along a ray direction at a given position by using a quadrature rule
// Uses an implicit total derivative (TODO: ref to own derivation)
//Warning: In the case we are on the boundary, it assumes we can find two points in a single direction
//Warning: assumes we use the same grid in both the original code and the collocation code (is valid, unless the original grid has changed after)
inline Real Collocation :: get_source_grad(Model& model, Size rayid, Size freqid, Size pointid)
{
    const Size point_in_grid=points_in_grid[pointid];
    vector<Size> nearby_points=get_two_nearby_points_on_ray(model, rayid, point_in_grid);
    Size next_point=nearby_points[0];
    Size prev_point=nearby_points[1];

    if (next_point==parameters.npoints()||prev_point==parameters.npoints())//if at least one of the nearby points do no exist, then no gradient can be computed
    {//I do not expect this to happen at all, but either way,
        std::cout<<"A ray through a point does not go through at least 3 points"<<std::endl;
        std::cout<<"pointid: "<<pointid<<std::endl;
        std::cout<<"rayid: "<<rayid<<std::endl;
        std::cout<<"next point: "<<next_point<<std::endl;
        std::cout<<"prev point: "<<prev_point<<std::endl;
        return 0;
    }

    Vector3D raydirection=model.geometry.rays.direction[rayid];
    Vector3D next_vector=model.geometry.points.position[next_point]-model.geometry.points.position[point_in_grid];
    Vector3D prev_vector=model.geometry.points.position[prev_point]-model.geometry.points.position[point_in_grid];
    const Real next_position=raydirection.dot(next_vector);
    const Real prev_position=raydirection.dot(prev_vector);

    const Real position_frac=next_position/prev_position;
    Real prev_coef=1.0/(prev_position*(1.0-(1.0/position_frac)));
    Real next_coef=1.0/(next_position*(1.0-position_frac));
    Real curr_coef=-prev_coef-next_coef;

    Real prev_opacity=get_opacity(model, rayid, freqid, prev_point);//<-Warning: same grid assumed
    Real prev_emissivity=get_emissivity(model, rayid, freqid, prev_point);//<-Warning: same grid assumed
    Real curr_opacity=get_opacity(model, rayid, freqid, pointid);
    Real curr_emissivity=get_emissivity(model, rayid, freqid, pointid);
    Real next_opacity=get_opacity(model, rayid, freqid, next_point);//<-Warning: same grid assumed
    Real next_emissivity=get_emissivity(model, rayid, freqid, next_point);//<-Warning: same grid assumed
    //comoving source functions are used

    return prev_coef*prev_emissivity/prev_opacity
          +curr_coef*curr_emissivity/curr_opacity
          +next_coef*next_emissivity/next_opacity;
}



// Returns the emissivity (takes into account the gaussian line profile)
inline Real Collocation :: get_emissivity(Model& model, Size rayid, Size freqid, Size pointid)
{
    Real toreturn=0;
    const Size point_in_grid=points_in_grid[pointid];
    // if using other line profile/frequency basis function than gaussian, please replace the for(...) with the following line:
    // for (Size other_line_idx=0; other_line_idx<parameters.nlines(); other_line_idx++)
    // for (Size other_line_idx: other_frequency_basis_interacting_with[pointidx][lineidx])
    for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
    {
        Real delta_freq=0;
        // Real inv_width=0;
        if (USING_COMOVING_FRAME)
        {
            // inv_width=model.lines.inverse_width(pointidx,lineidx);
            delta_freq=non_doppler_shifted_frequencies[freqid]-model.lines.line[lineidx];
        }
        else
        {
          // Real doppler_shift_factor=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[pointidx],model.geometry.rays.direction[rayidx]);
          //local opacity should use the local frequencies, not the global ones...; otherwise it goes very quickly to zero if we have a minor doppler shift during the ray
          // I just don't feel like doppler shifting everything -> only minor correction in the inv_width factor in front is countered by the doppler shifted line frequency
          delta_freq=non_doppler_shifted_frequencies[freqid]-model.lines.line[lineidx];
          // delta_freq=frequencies[rayidx][pointidx][freqidx]-model.lines.line[lineidx]*doppler_shift_factor;
          // inv_width=model.lines.inverse_width(pointidx,lineidx)/doppler_shift_factor;
        }
        const Real inv_width=model.lines.inverse_width(point_in_grid,lineidx);
        toreturn+=INVERSE_SQRT_PI*inv_width*std::exp(-std::pow(delta_freq*inv_width,2))
                  *model.lines.line[lineidx]*model.lines.emissivity(point_in_grid, lineidx);
    }
    return toreturn;
}

inline void Collocation :: setup_rhs_Eigen(Model& model)
{
    rhs = VectorXr::Zero (collocation_mat.rows());
    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {
            for (Size pointid=0; pointid<n_points_in_grid; pointid++)
            {
                const Size point_in_grid=points_in_grid[pointid];
                const Size curr_mat_row_idx=get_mat_row_index(rayid, freqid, pointid);
                //FIXME: check if boundary condition needs to apply

                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(point_in_grid))
                {
                    const bool bdy_condition_req_in_direction=is_boundary_condition_needed_for_direction_on_boundary(model, rayid, point_in_grid);
                    boundary_condition_required=bdy_condition_req_in_direction;

                    // double dist=0;//these two variables are actually not used, but required in the function definition get_next
                    // double ddist=0;
                    //
                    // Size antipod=model.geometry.rays.antipod[rayidx];
                    //
                    // const Size point_behind=model.geometry.get_next(pointidx, antipod, pointidx, dist, ddist);
                    // const Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                    // if (point_behind==parameters.npoints())
                    // {//if there lies no point behind this point, we need to evaluate the boundary condition here
                    //     boundary_condition_required=true;
                    //
                    //     // dist=0;//these two variables are actually not used, but required in the function definition get_next
                    //     // ddist=0;
                    //     // //also compute velocity gradient for the balancing factor
                    //     // Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                    //     // std::cout<<"next point: "<<next_point<<std::endl;
                    //     // if (next_point!=parameters.npoints())//if there actually exists a next point
                    //     // {
                    //     //     std::cout<<"computing velocity grad"<<std::endl;
                    //     //     local_velocity_gradient=model.geometry.rays.direction[rayidx].dot(model.geometry.points.velocity[next_point]-model.geometry.points.velocity[pointidx])/dist;
                    //     // }
                    // }
                }
                // boundary_condition_required=false;
                if (boundary_condition_required)
                {
                    Real boundary_condition_mult_factor=boundary_condition_LS_importance_factor;
                    rhs[curr_mat_row_idx]=boundary_intensity(model, rayid, freqid, pointid)*boundary_condition_mult_factor;
                }
                else
                {
                    //TODO: possibly some doppler shift required here? //no, we do not interpolate these values
                    rhs[curr_mat_row_idx]=get_emissivity(model, rayid, freqid, pointid)/get_opacity(model, rayid, freqid, pointid);
                    std::cout<<"emissivity: "<<get_emissivity(model, rayid, freqid, pointid)<<std::endl;
                    std::cout<<"opacity: "<<get_opacity(model, rayid, freqid, pointid)<<std::endl;
                }
            }
        }
    }
    std::cout<<"rhs: "<<rhs<<std::endl;

    rescale_matrix_and_rhs_Eigen(model);

    std::cout<<"rescaled rhs: "<<rhs<<std::endl;
}

// Sets up the right-hand side for the square matrix
inline void Collocation :: setup_rhs_square_Eigen(Model& model)
{
    rhs = VectorXr::Zero (collocation_mat.rows());
    std::cout<<"nb collocation rows: "<<collocation_mat.rows()<<std::endl;
    for (Size rayid=0; rayid<parameters.nrays(); rayid++)
    {
        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {
            for (Size point_basis_id=0; point_basis_id<n_point_bases[freqid]; point_basis_id++)
            {
                const Size pointid=local_basis_point_idx[freqid][point_basis_id];
                const Size point_in_grid=basis_point_idx[freqid][point_basis_id];
                // const Size curr_mat_row_idx=get_mat_row_index(rayid, freqid, pointid);
                const Size curr_mat_row_idx=get_mat_col_index(rayid, freqid, point_basis_id);//As we want to make a square matrix, we replace the row indices with column indices
                //FIXME: check if boundary condition needs to apply

                bool boundary_condition_required=false;
                if (!model.geometry.not_on_boundary(point_in_grid))
                {
                    const bool bdy_condition_req_in_direction=is_boundary_condition_needed_for_direction_on_boundary(model, rayid, point_in_grid);
                    boundary_condition_required=bdy_condition_req_in_direction;

                    // double dist=0;//these two variables are actually not used, but required in the function definition get_next
                    // double ddist=0;
                    //
                    // Size antipod=model.geometry.rays.antipod[rayidx];
                    //
                    // const Size point_behind=model.geometry.get_next(pointidx, antipod, pointidx, dist, ddist);
                    // const Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                    // if (point_behind==parameters.npoints())
                    // {//if there lies no point behind this point, we need to evaluate the boundary condition here
                    //     boundary_condition_required=true;
                    //
                    //     // dist=0;//these two variables are actually not used, but required in the function definition get_next
                    //     // ddist=0;
                    //     // //also compute velocity gradient for the balancing factor
                    //     // Size next_point=model.geometry.get_next(pointidx, rayidx, pointidx, dist, ddist);
                    //     // std::cout<<"next point: "<<next_point<<std::endl;
                    //     // if (next_point!=parameters.npoints())//if there actually exists a next point
                    //     // {
                    //     //     std::cout<<"computing velocity grad"<<std::endl;
                    //     //     local_velocity_gradient=model.geometry.rays.direction[rayidx].dot(model.geometry.points.velocity[next_point]-model.geometry.points.velocity[pointidx])/dist;
                    //     // }
                    // }
                }
                // boundary_condition_required=false;
                if (boundary_condition_required)
                {
                    rhs[curr_mat_row_idx]=boundary_intensity(model, rayid, freqid, pointid);
                }
                else
                {
                    //TODO: possibly some doppler shift required here? //no, we do not interpolate these values
                    rhs[curr_mat_row_idx]=get_emissivity(model, rayid, freqid, pointid)/get_opacity(model, rayid, freqid, pointid);
                    std::cout<<"emissivity: "<<get_emissivity(model, rayid, freqid, pointid)<<std::endl;
                    std::cout<<"opacity: "<<get_opacity(model, rayid, freqid, pointid)<<std::endl;
                }
            }
        }
    }
    std::cout<<"rhs: "<<rhs<<std::endl;

    rescale_matrix_and_rhs_Eigen(model);

    std::cout<<"rescaled rhs: "<<rhs<<std::endl;
}

//Sets up the rhs correctly for the second order feautrier version of the collocation method
inline void Collocation :: setup_rhs_Eigen_2nd_feautrier(Model& model)
{
    rhs = VectorXr::Zero (collocation_mat.rows());
    for (Size rayid=0; rayid<parameters.hnrays(); rayid++)//for u and v at the same time
    {
        for (Size freqid=0; freqid<n_freq_bases; freqid++)
        {
            for (Size pointid=0; pointid<n_points_in_grid; pointid++)
            {
                const Size point_in_grid=points_in_grid[pointid];
                const Size u_mat_row_idx=get_mat_row_index_2nd_feautrier(rayid, freqid, pointid, false);
                const Size v_mat_row_idx=get_mat_row_index_2nd_feautrier(rayid, freqid, pointid, true);
                //FIXME: check if boundary condition needs to apply
                Real local_velocity_gradient=0;
                bool boundary_condition_required=false;
                bool bdy_condition_req_in_direction=false;
                if (!model.geometry.not_on_boundary(point_in_grid))//TODO: refactor this (instead of extra if-clause)
                {
                    bdy_condition_req_in_direction=is_boundary_condition_needed_for_direction_on_boundary(model, rayid, point_in_grid);
                    boundary_condition_required=true;

                    // double dist=0;//these two variables are actually not used, but required in the function definition get_next
                    // double ddist=0;
                    //
                    // Size antipod=model.geometry.rays.antipod[rayidx];
                    //
                    // const Size point_behind=model.geometry.get_next(pointidx, antipod, pointidx, dist, ddist);
                    // if (point_behind==parameters.npoints())
                    // {//if there lies no point behind this point, we need to evaluate the boundary condition here
                    //     boundary_condition_required=true;
                    // }
                }
                // boundary_condition_required=false;
                if (boundary_condition_required)//boundary conditions only applied at the u equations, but also slightly different equation for v
                {
                    // std::cout<<"local velocity grad: "<<local_velocity_gradient<<std::endl;
                    //This term makes sure that the magnitude of the matrix terms corresponding to this boundary condition is similar to the other nearby terms
                    // const Real balancing_factor=compute_balancing_factor_boundary(rayidx, freqidx, pointidx, local_velocity_gradient, model);
                    Real bdy_intensity=boundary_intensity(model, rayid, freqid, pointid);//boundary intensity
                    rhs[u_mat_row_idx]=bdy_intensity;
                    Real source=get_emissivity(model, rayid, freqid, pointid)/get_opacity(model, rayid, freqid, pointid);//source function
                    //TODO: refactor such that no if statement is required
                    if (bdy_condition_req_in_direction)
                    {
                        rhs[v_mat_row_idx]=bdy_intensity-source;
                    }
                    else
                    {
                        rhs[v_mat_row_idx]=source-bdy_intensity;
                    }
                }
                else
                {
                    //TODO: possibly some doppler shift required here? //no, we do not interpolate these values
                    rhs[u_mat_row_idx]=get_emissivity(model, rayid, freqid, pointid)/get_opacity(model, rayid, freqid, pointid);
                    // std::cout<<"emissivity: "<<get_emissivity(model, rayidx, freqidx, pointidx)<<std::endl;
                    // std::cout<<"opacity: "<<get_opacity(model, rayidx, freqidx, pointidx)<<std::endl;
                    rhs[v_mat_row_idx]=-get_source_grad(model, rayid, freqid, pointid)/get_opacity(model, rayid, freqid, pointid);
                    std::cout<<"-source fun grad/opacity= "<<rhs[v_mat_row_idx]<<std::endl;
                }

            }
        }
    }
    std::cout<<"rhs: "<<rhs<<std::endl;

    rescale_matrix_and_rhs_Eigen(model);

    std::cout<<"rescaled rhs: "<<rhs<<std::endl;
}


// Rescale the matrix rows such that there are ones on the diagonal
// Practically, this is the same as applying a the inverse of the diagonal as left preconditioner
// TODO: the explicit multiplication factor for the boundary conditions might be a bit useless now
//Assumes that both the eigen matrix and rhs have already been constructed
inline void Collocation :: rescale_matrix_and_rhs_Eigen(Model& model)
{
    //FIXME: in least squares, we need to normalize the columns instead
    //So collect all corresponding 'diagonal elements' (for each basis the evaluation centered at that basis)
    //Then only normalize the matrix, not the rhs; after computation rescale the basis coefficients

    // std::cout<<"rescaling matrix"<<std::endl;
    // for (Size basis_ray_id=0; basis_ray_id<parameters.nrays(); basis_ray_id++)//for second order feautrier, we technically need to separate the u and v bases, but their placement in the matrix follows the smae pattern
    // {
    //     for (Size basis_freq_id=0; basis_freq_id<n_freq_bases; basis_freq_id++)
    //     {
    //         for (Size basis_point_id=0; basis_point_id<n_point_bases[basis_freq_id]; basis_point_id++)
    //         {
    //             Size converted_point_id=local_basis_point_idx[basis_freq_id][basis_point_id];
    //             const Size mat_row_idx=get_mat_row_index(basis_ray_id, basis_freq_id, converted_point_id);
    //             const Size mat_col_idx=get_mat_col_index(basis_ray_id, basis_freq_id, basis_point_id);
    //             const Real centered_eval=collocation_mat.coeff(mat_row_idx,mat_col_idx);//evaluation of the basis function in the central point
    //             // basis_renorm_factor(mat_col_idx)=1.0/centered_eval;
    //             basis_renorm_factor.insert(mat_col_idx,mat_col_idx)=1.0/centered_eval;
    //         }
    //     }
    // }
    // std::cout<<"finished rescaling"<<std::endl;

    // //get inverse of diagonal as diagonal matrix
    // //FIXME: check whether we have zero on the diagonal; if yes, abort
    const auto diagonal_inverse=collocation_mat.diagonal().asDiagonal().inverse();
    std::cout<<"diag inverse rows: "<<diagonal_inverse.rows()<<std::endl;
    std::cout<<"diag inverse cols: "<<diagonal_inverse.cols()<<std::endl;
    std::cout<<"diagonal: "<<diagonal_inverse.diagonal()<<std::endl;
    rhs=diagonal_inverse*rhs;
    std::cout<<"rhs size: "<<rhs.size()<<std::endl;
    collocation_mat=diagonal_inverse*collocation_mat;

    std::cout<<"collocation rows: "<<collocation_mat.rows()<<std::endl;
    std::cout<<"collocation cols: "<<collocation_mat.cols()<<std::endl;
    // std::cout<<"nb basis renorm factors: "<<basis_renorm_factor.rows()<<", "<<basis_renorm_factor.cols()<<std::endl;
    // std::cout<<"basis renorm factor size: "<<basis_renorm_factor.size()<<std::endl;

    // collocation_mat=collocation_mat*(basis_renorm_factor.asDiagonal());
    // SparseMatrix<Real> temp_collocation_mat=collocation_mat*basis_renorm_factor;
    // collocation_mat=collocation_mat*basis_renorm_factor;
    // std::cout<<"test"<<std::endl;
    // collocation_mat=temp_collocation_mat;

    // std::cout<<"after rescaling mat"<<std::endl;

    //also printing out the rescaled version
    if (model.COMPUTING_SVD_AND_CONDITION_NUMBER)
    {
        Eigen::JacobiSVD<MatrixXr> svd(collocation_mat);
        std::cout<<"here"<<std::endl;
        Real cond = svd.singularValues()(0)
        / svd.singularValues()(svd.singularValues().size()-1);
        std::cout<<"condition number: "<<cond<<std::endl;
        // collocation_mat.data().squeeze();
        std::cout << MatrixXr(collocation_mat) << std::endl;
        std::cout<<"basis renorm factors: "<<VectorXr(basis_renorm_factor.diagonal())<<std::endl;
    }

}

///  Getter for the boundary conditions
///    @param[in] model  : reference to model object
///    @param[in] rayidx: ray index at which to evaluate the boundary condition
///    @param[in] freqidx: line frequency index at which to evaluate boundary condition
///    @param[in] pointidx: point index of the boundary point in the collocation grid
///    @returns incoming radiation intensity at the boundary
////////////////////////////////////////////////////////////////////////////
///  Copied from solver.tpp and modified a bit
inline Real Collocation :: boundary_intensity (Model& model, Size rayid, Size freqid, Size pointid)
{
    const Size point_in_grid=points_in_grid[pointid];
    Real freq=0;
    if (USING_COMOVING_FRAME)
    {
        freq=non_doppler_shifted_frequencies[freqid];
    }
    else
    {
        freq=frequencies[rayid][pointid][freqid];
    }

    const Size bdy_id = model.geometry.boundary.point2boundary[point_in_grid];

    switch (model.geometry.boundary.boundary_condition[bdy_id])
    {
        case Zero    : return 0.0;
        case Thermal : return planck (model.geometry.boundary.boundary_temperature[bdy_id], freq);
        default      : return planck (T_CMB, freq);
    }
}


///  Planck function
///    @param[in] temp : temperature of the corresponding black body
///    @param[in] freq : frequency at which to evaluate the function
///    @return Planck function evaluated at this frequency
///////////////////////////////////////////////////////////////////////////
///  Copied from solver.tpp
accel inline Real Collocation :: planck (Real temp, Real freq)
{
    return TWO_HH_OVER_CC_SQUARED * (freq*freq*freq) / expm1 (HH_OVER_KB*freq/temp);
}

// threaded_for (p, parameters.npoints(),
// {
//     for (Size l = 0; l < parameters.nlspecs(); l++)
//     {
//         for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
//         {
//             const Size lid = line_index (l, k);
//             emissivity (p, lid) = lineProducingSpecies[l].get_emissivity (p, k);
//                opacity (p, lid) = lineProducingSpecies[l].get_opacity    (p, k);
//         }
//     }
// })

inline void Collocation :: solve_collocation(Model& model)
{
    // std::cout<<"Eigen n threads: "<<Eigen::nbThreads()<<std::endl;
    // Eigen::BiCGSTAB<SparseMatrix<Real> > solver;//does not work on rectangular matrices
    //Either use some iterative method
    // Eigen::LeastSquaresConjugateGradient<SparseMatrix<Real>> solver;
    //Or use a direct solver
    // Eigen::SparseQR<SparseMatrix<Real>> solver;
    //TOOD: benchmark the two solvers versus each other!!!
    // solver.compute(collocation_mat);
    // basis_coefficients = solver.solve(rhs);
    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error()      << std::endl;

    // std::cout<<"is matrix compressed before makeCompressed: "<<collocation_mat.isCompressed()<<std::endl;
    // //SparseQR requires the matrix to be in compressed form. But set_from_triplets already puts it in compressed form

    Eigen :: SparseLU <SparseMatrix<Real>, COLAMDOrdering<int>> solver;

    std::cout<<"flag 1"<<std::endl;

    if (USING_LEAST_SQUARES)//then solve A^TA=A^Tb
    {
        auto col_mat_transpose=collocation_mat.transpose();
        //for solving least squares, we can just multiply both the rhs and the matrix by the transpose (leads to pseudoinverse)
        VectorXr rescaled_rhs=col_mat_transpose*rhs;
        SparseMatrix<Real> symm_collocation_mat=col_mat_transpose*collocation_mat;


        const auto diagonal_inverse=symm_collocation_mat.diagonal().asDiagonal().inverse();
        rescaled_rhs=diagonal_inverse*rescaled_rhs;
        symm_collocation_mat=diagonal_inverse*symm_collocation_mat;

        cout << "Analyzing system..."      << endl;
        solver.analyzePattern (symm_collocation_mat);

        cout << "Factorizing system..."    << endl;
        solver.factorize(symm_collocation_mat);

        if (solver.info() != Eigen::Success)
        {
            cout << "Factorization failed with error message:" << endl;
            cout << solver.lastErrorMessage()                  << endl;

            throw std::runtime_error ("Eigen solver ERROR.");
        }

        cout << "Solving equations..." << endl;

        basis_coefficients = solver.solve(rescaled_rhs);

        if (solver.info() != Eigen::Success)
        {
            cout << "Solving failed with error:" << endl;
            cout << solver.lastErrorMessage()    << endl;
            assert (false);
        }

        cout << "error: " << (symm_collocation_mat*basis_coefficients-rescaled_rhs).norm()/rhs.norm() << endl;
    }
    else
    {
        cout << "Analyzing system..."      << endl;
        solver.analyzePattern (collocation_mat);
        cout << "Factorizing system..."    << endl;

        solver.factorize (collocation_mat);

        if (solver.info() != Eigen::Success)
        {
            cout << "Factorization failed with error message:" << endl;
            cout << solver.lastErrorMessage()                  << endl;

            throw std::runtime_error ("Eigen solver ERROR.");
        }

        cout << "Solving equations..." << endl;

        basis_coefficients = solver.solve(rhs);

        if (solver.info() != Eigen::Success)
        {
            cout << "Solving failed with error:" << endl;
            cout << solver.lastErrorMessage()    << endl;
            assert (false);
        }

        cout << "error: " << (collocation_mat*basis_coefficients-rhs).norm()/rhs.norm() << endl;
    }


    // if (model.COMPUTING_SVD_AND_CONDITION_NUMBER)
    // {
    //     Eigen::JacobiSVD<MatrixXr> svd(symm_collocation_mat);
    //     std::cout<<"here"<<std::endl;
    //     Real cond = svd.singularValues()(0)
    //     / svd.singularValues()(svd.singularValues().size()-1);
    //     std::cout<<"condition number: "<<cond<<std::endl;
    //     // collocation_mat.data().squeeze();
    //     std::cout << MatrixXr(symm_collocation_mat) << std::endl;
    //     std::cout<<"basis renorm factors: "<<VectorXr(basis_renorm_factor.diagonal())<<std::endl;
    // }
    //
    // symm_collocation_mat.makeCompressed();


    // Err, it is unfortunately only reliable if we have some SPD property
    // Eigen :: LeastSquaresConjugateGradient<SparseMatrix<Real>> solver;
    // solver.compute(collocation_mat);
    // basis_coefficients=solver.solve(rhs);
    //
    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error()      << std::endl;

    // Eigen :: SparseQR <SparseMatrix<Real>, COLAMDOrdering<int>> solver;

    Size vectorsize=collocation_mat.cols();
    std::cout<<"computed coeffs"<<std::endl;
    for (Size i=0; i<vectorsize; i++)
    {
      std::cout<<"i: "<<i<<"value: "<<basis_coefficients(i)<<std::endl;
    }

    //denormalize the coefficients
    // basis_coefficients=basis_renorm_factor.asDiagonal()*basis_coefficients;
    // basis_coefficients=basis_renorm_factor*basis_coefficients;
    // std::cout<<"denormalized coeffs"<<std::endl;
    // for (Size i=0; i<vectorsize; i++)
    // {
    //   std::cout<<"i: "<<i<<"value: "<<basis_coefficients(i)<<std::endl;
    // }


}

// //TODO: calculate for all quadrature points the intensity; or exactly integrate the gaussian with another gaussian
inline void Collocation :: compute_J(Model& model)
{
  model.radiation.initialize_J();

  // //For testing purposes, initialize I too
  // for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
  // {
  //     for (Size pointidx=0; pointidx<parameters.npoints(); pointidx++)
  //     {
  //         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
  //         {
  //             model.radiation.I(rayidx, pointidx, lineidx)=0.0;
  //         }
  //     }
  // }
  //
  // for (Size rayidx=0; rayidx<parameters.hnrays(); rayidx++)
  // {
  //     for (Size pointidx=0; pointidx<parameters.npoints(); pointidx++)
  //     {
  //         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
  //         {
  //             model.radiation.u(rayidx, pointidx, lineidx)=0.0;
  //         }
  //     }
  // }

  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
      for (Size pointid=0; pointid<n_points_in_grid; pointid++)
      // threaded_for (p, parameters.npoints(),
      {
          const Size point_in_grid=points_in_grid[pointid];
          // const Size p=points_in_grid[idx];

          for (Size k = 0; k < lspec.linedata.nrad; k++)
          {
              // const Size1 freq_nrs = lspec.nr_line[p][k];
              const Size lineid = model.lines.line_index(l, k);
              const Size freqid_corresponding_to_line_center=center_line_freq_indices[lineid];
              // const Size1 freq_nrs = lspec.nr_line[point_in_grid][k];
              //const Size middle_linefreqidx=freq_nrs[parameters.nquads()/2];//for determining which frequencies interact with the line frequency; FIXME: remove this hack by explicitly calculating this
              //instead we might just use another slightly less hack: namely checking which ones interact with the center line frequency by checking which bases have a non-zero evalutation at a point at the center line frequency
              //FIXME: Instead compute explicitly which basis functions interact with the line center?

              // Initialize values
              lspec.Jlin[point_in_grid][k] = 0.0;

              // // Just sum for all directions over the non-zero basis functions
              // std::set<std::vector<Size>> nonzero_basis_triplets;
              // get_nonzero_basis_triplets(nonzero_basis_triplets, rayidx, lid, pointidx);


              std::cout<<"Started computing a single J"<<std::endl;
              //or just use a matrix-vector product
              Real average_I=0;
              for (Size rayid=0; rayid<parameters.nrays(); rayid++)
              {
                  Real lp_freq=0;
                  Real lp_freq_inv_width=0;
                  if (USING_COMOVING_FRAME)
                  {
                      // lp_freq=non_doppler_shifted_frequencies[lid];
                      // lp_freq_inv_width=non_doppler_shifted_inverse_widths[lid][p];
                      lp_freq=model.lines.line[lineid];
                      lp_freq_inv_width=model.lines.inverse_width(lineid,point_in_grid);
                  }
                  else
                  {
                      const Real doppler_shift_factor=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[point_in_grid],model.geometry.rays.direction[rayid]);
                      // std::cout<<"doppler shift factor: "<<doppler_shift_factor<<std::endl;
                      lp_freq=model.lines.line[lineid]*doppler_shift_factor;
                      lp_freq_inv_width=model.lines.inverse_width(lineid,point_in_grid)/doppler_shift_factor;
                  }
                  //testing out whether the inverse width is correct//seems to be correct
                  // lp_freq_inv_width=1.0/std::sqrt(2)*lp_freq_inv_width;
                  // Size hrayidx=rayidx;//rayidx for u
                  // if (hrayidx>=parameters.hnrays())
                  // {
                  //     hrayidx=model.geometry.rays.antipod[hrayidx];
                  // }

                  std::set<std::vector<Size>> nonzero_basis_triplets;

                  get_nonzero_basis_triplets(nonzero_basis_triplets, rayid, freqid_corresponding_to_line_center, pointid);
                  Real I=0;
                  Real temp_jlin=0;
                  for (std::vector<Size> triplet: nonzero_basis_triplets)
                  // for (Size z = 0; z < parameters.nquads(); z++)
                  {
                      const Size triplet_basis_index=get_mat_col_index(triplet[0], triplet[1], triplet[2]);
                      if (USING_COMOVING_FRAME)
                      {
                          // const Real doppler_shifted_frequency=lp_freq*calculate_doppler_shift_observer_frame(model.geometry.points.velocity[triplet[2]]-model.geometry.points.velocity[p]
                          //                                                                                       , model.geometry.rays.direction[rayidx]);
                          // Real doppler_shifted_frequency=lp_freq;//try it out non-doppler shifted

                          I+=basis_coefficients(triplet_basis_index)*basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          lspec.Jlin[point_in_grid][k]+=basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                          // I+=basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency);
                          // lspec.Jlin[p][k]+=basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], doppler_shifted_frequency, lp_freq_inv_width);
                      }
                      else
                      {
                          //TODO: figure out which values we actually want to keep...; note to self: I and u use the frequency index (nlines*nquads); so we are not able to compute it here
                          // model.radiation.I(rayidx, p, lid) +=basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          // model.radiation.u(hrayidx, p, lid)+=1.0/2*basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          //FIXME: doppler shift the line profile frequency such that it corresponds to the local frequency!!!!!!!!!

                          I+=basis_coefficients(triplet_basis_index)*basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                          lspec.Jlin[point_in_grid][k]+=basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                          // std::cout<<"The summation of this is I: "<<basis_coefficients(triplet_basis_index)*basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq)<<std::endl;
                          // std::cout<<"The summation of this should be J: "<<basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width)<<std::endl;

                          // lspec.Jlin[p][k]+=200000*basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p])*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                      }

                  }
                  std::cout<<"pointid: "<<pointid<<" I: "<<I<<std::endl;
                  std::cout<<"J: "<<lspec.Jlin[point_in_grid][k]<<std::endl;
                  average_I+=I/parameters.nrays();

              }

              std::cout<<"average I: "<<average_I<<std::endl;
              double diff = 0.0;

              // // Collect the approximated part
              // for (Size m = 0; m < lspec.lambda.get_size(p,k); m++)
              // {
              //     const Size I = lspec.index(lspec.lambda.get_nr(p,k,m), lspec.linedata.irad[k]);
              //
              //     diff += lspec.lambda.get_Ls(p,k,m) * lspec.population[I];
              // }

              lspec.Jeff[point_in_grid][k] = lspec.Jlin[point_in_grid][k] - HH_OVER_FOUR_PI * diff;
              lspec.Jdif[point_in_grid][k] = HH_OVER_FOUR_PI * diff;

              // if (lspec.Jeff[p][k]<0)
              // {
              //     lspec.Jeff[p][k]=0;
              // }
              // std::cout<<"Jeff: "<<lspec.Jeff[p][k]<<std::endl;
          }
      }//)
  }
  //print out j for first transition
  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
      for (Size pointid=0; pointid<n_points_in_grid; pointid++)
      {
          const Size point_in_grid=points_in_grid[pointid];
          std::cout<<"jeff: "<<lspec.Jeff[point_in_grid][0]<<std::endl;
      }
  }
}

// //TODO: calculate for all quadrature points the intensity; or exactly integrate the gaussian with another gaussian
inline void Collocation :: compute_J_2nd_feautrier(Model& model)
{
  model.radiation.initialize_J();

  // //For testing purposes, initialize I too
  // for (Size rayidx=0; rayidx<parameters.nrays(); rayidx++)
  // {
  //     for (Size pointidx=0; pointidx<parameters.npoints(); pointidx++)
  //     {
  //         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
  //         {
  //             model.radiation.I(rayidx, pointidx, lineidx)=0.0;
  //         }
  //     }
  // }
  //
  // for (Size rayidx=0; rayidx<parameters.hnrays(); rayidx++)
  // {
  //     for (Size pointidx=0; pointidx<parameters.npoints(); pointidx++)
  //     {
  //         for (Size lineidx=0; lineidx<parameters.nlines(); lineidx++)
  //         {
  //             model.radiation.u(rayidx, pointidx, lineidx)=0.0;
  //         }
  //     }
  // }

  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
      for (Size pointid=0; pointid<n_points_in_grid; pointid++)
      // threaded_for (p, parameters.npoints(),
      {
          const Size point_in_grid=points_in_grid[pointid];
          // const Size p=points_in_grid[idx];

          for (Size k = 0; k < lspec.linedata.nrad; k++)
          {
              const Size lineid = model.lines.line_index(lineid, k);
              const Size freqid_corresponding_to_line_center=center_line_freq_indices[lineid];
              // const Size1 freq_nrs = lspec.nr_line[p][k];
              // const Size middle_linefreqidx=freq_nrs[parameters.nquads()/2];//for determining which frequencies interact with the line frequency; FIXME: remove this hack by explicitly calculating this

              // Initialize values
              lspec.Jlin[point_in_grid][k] = 0.0;

              //Sum over ray indices for u
              for (Size rayid=0; rayid<parameters.hnrays(); rayid++)
              {

                  const Real doppler_shift_factor=calculate_doppler_shift_observer_frame(model.geometry.points.velocity[point_in_grid],model.geometry.rays.direction[rayid]);
                  // std::cout<<"doppler shift factor: "<<doppler_shift_factor<<std::endl;
                  Real lp_freq=model.lines.line[lineid]*doppler_shift_factor;
                  Real lp_freq_inv_width=model.lines.inverse_width(lineid,point_in_grid)/doppler_shift_factor;

                  Size hrayid=rayid;//rayidx for u and v
                  if (hrayid>=parameters.hnrays())
                  {
                      hrayid=model.geometry.rays.antipod[hrayid];
                  }

                  std::set<std::vector<Size>> nonzero_basis_triplets;

                  get_nonzero_basis_triplets(nonzero_basis_triplets, rayid, freqid_corresponding_to_line_center, pointid);
                  Real Ihatn=0;
                  Real Iminhatn=0;
                  Real temp_jlin=0;
                  for (std::vector<Size> triplet: nonzero_basis_triplets)
                  // for (Size z = 0; z < parameters.nquads(); z++)
                  {
                      const Size triplet_u_basis_index=get_mat_col_index_2nd_feautrier(triplet[0], triplet[1], triplet[2], false);
                      const Size triplet_v_basis_index=get_mat_col_index_2nd_feautrier(triplet[0], triplet[1], triplet[2], true);
                      // std::cout<<"triplet_u_basis_index: "<<triplet_u_basis_index<<std::endl;
                      // std::cout<<"u coeff: "<<basis_coefficients(triplet_u_basis_index)<<std::endl;
                      // std::cout<<"triplet_v_basis_index: "<<triplet_v_basis_index<<std::endl;
                      // std::cout<<"v coeff: "<<basis_coefficients(triplet_v_basis_index)<<std::endl;



                      //TODO: figure out which values we actually want to keep...; note to self: I and u use the frequency index (nlines*nquads); so we are not able to compute it here
                      // TODO: rewrite using second order feautrier vars
                      // model.radiation.I(rayidx, p, lid) +=basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                      // model.radiation.u(hrayidx, p, lid)+=1.0/2*basis_coefficients(triplet_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);

                      //same bases for u and v
                      Real basis_eval=basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                      Real sum_I_div_2=basis_coefficients(triplet_u_basis_index)*basis_eval;//basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                      Real dif_I_div_2=basis_coefficients(triplet_v_basis_index)*basis_eval;//basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                      Ihatn+=sum_I_div_2+dif_I_div_2;//basis_coefficients(triplet_u_basis_index)*basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction(triplet[0])*basis_freq(triplet[0], triplet[1], triplet[2], lp_freq);
                      Iminhatn+=sum_I_div_2-dif_I_div_2;

                      Real basis_dir_int_lp_int_eval=1/FOUR_PI*basis_point(triplet[0], triplet[1], triplet[2], point_locations[pointid], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                      lspec.Jlin[point_in_grid][k]+=2*basis_coefficients(triplet_u_basis_index)*basis_dir_int_lp_int_eval;//basis_coefficients(triplet_basis_index)/FOUR_PI*basis_point(triplet[2], point_locations[p], triplet[0], triplet[1], model.geometry)*basis_direction_int(model.geometry,triplet[0])*basis_freq_lp_int(triplet[0], triplet[1], triplet[2], lp_freq, lp_freq_inv_width);
                      //contribution from I(hatn)
                      //lspec.Jlin[p][k]+=(basis_coefficients(triplet_u_basis_index)+basis_coefficients(triplet_v_basis_index))*basis_dir_int_lp_int_eval;
                      //contribution from I(-hatn)
                      //lspec.Jlin[p][k]+=(basis_coefficients(triplet_u_basis_index)-basis_coefficients(triplet_v_basis_index))*basis_dir_int_lp_int_eval;

                  }
                  std::cout<<"pointid: "<<pointid<<" Ihatn: "<<Ihatn<<std::endl;
                  std::cout<<"pointid: "<<pointid<<" Iminhatn: "<<Iminhatn<<std::endl;

              }

              double diff = 0.0;

              // // Collect the approximated part
              // for (Size m = 0; m < lspec.lambda.get_size(p,k); m++)
              // {
              //     const Size I = lspec.index(lspec.lambda.get_nr(p,k,m), lspec.linedata.irad[k]);
              //
              //     diff += lspec.lambda.get_Ls(p,k,m) * lspec.population[I];
              // }

              lspec.Jeff[point_in_grid][k] = lspec.Jlin[point_in_grid][k] - HH_OVER_FOUR_PI * diff;
              lspec.Jdif[point_in_grid][k] = HH_OVER_FOUR_PI * diff;

              // if (lspec.Jeff[p][k]<0)
              // {
              //     lspec.Jeff[p][k]=0;
              // }
              // std::cout<<"Jeff: "<<lspec.Jeff[p][k]<<std::endl;
          }
      }//)
  }
  //print out j for first transition
  for (Size l = 0; l < parameters.nlspecs(); l++)
  {
      LineProducingSpecies& lspec=model.lines.lineProducingSpecies[l];
      for (Size pointid=0; pointid<n_points_in_grid; pointid++)
      {
          const Size point_in_grid=points_in_grid[pointid];
          std::cout<<"jeff: "<<lspec.Jeff[point_in_grid][0]<<std::endl;
      }
  }
}
