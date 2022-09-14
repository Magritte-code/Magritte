

/// Computes the cooling rates using the collisional approach
///   @param[in] model: the model for which to compute the cooling rates
inline void Cooling :: compute_cooling_collisional(Model& model)
{
    Temperature& temperature=model.thermodynamics.temperature;
    Double2 abundance=model.chemistry.species.abundance;

    //for every point, compute the cooling rate independently
    threaded_for(p, model.parameters->npoints(),
    {
        Real temp_cooling_rate=0.0;
        const Real tmp = temperature[p];

        for (LineProducingSpecies &lspec : model.lines.lineProducingSpecies)
        {
            //cooling is done by colliding with other species (in this formulation)
            //some simple setup, ensuring the collisional rates have been computed correctly
            for (CollisionPartner &colpar : lspec.linedata.colpar)
            {
                Real abn = abundance[p][colpar.num_col_partner];

                colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
                colpar.interpolate_collision_coefficients (tmp);

                for (Size k = 0; k < colpar.ncol; k++)
                {
                    const Real collision_rate_ij = colpar.Cd_intpld[k] * abn; //deexcitation->heating
                    const Real collision_rate_ji = colpar.Ce_intpld[k] * abn; //excitation  ->cooling

                    const Size I = lineProducingSpecies[l].index (p, colpar.icol[k]);
                    const Size J = lineProducingSpecies[l].index (p, colpar.jcol[k]);

                    const Real energy_diff=lineProducingSpecies[l].linedata.energy[colpar.icol[k]]-lineProducingSpecies[l].linedata.energy[colpar.jcol[k]];
                    // std::cout<<"energy_diff: "<<energy_diff<<std::endl;

                    temp_cooling_rate-=lineProducingSpecies[l].population[I]*collision_rate_ij*energy_diff;
                    temp_cooling_rate+=lineProducingSpecies[l].population[J]*collision_rate_ji*energy_diff;
                }
            }
        }
        cooling_rate[p]=temp_cooling_rate;
    });
}


/// Computes the cooling rates using the flux definition approach
///   @param[in] model: the model for which to compute the cooling rates
/// Warning: cannot compute LTE cooling rates
inline void Cooling :: compute_cooling_flux(Model& model)
{
    //FIXME: For now, we assume that the integrated mean intensity has been correctly computed
    // SO IMPLEMENT THIS (after discussing with frederik on how to handle post-processing on rays)
    //OR use the mean line intensity as a substitute
    // Matrix<Real>& J=model.radiation.J; Err, is not used in every solver..., also contains the mean (not line) intensity at each frequency
    //TODO: after refactor, add options to enable cooling computation for each solver
    Matrix<Real>& emissivity=model.lines.emissivity;
    Matrix<Real>& opacity=model.lines.opacity;
    //for every point, compute the cooling rate independently
    threaded_for(p, model.parameters->npoints(),
    {
        Real temp_cooling_rate=0.0;
        //err, we use the mean line intensity instead
        for (Size l = 0; l < model.parameters->nlspecs(); l++)
        {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Size lid = line_index (l, k);

                const Real line_emission = emissivity (p, lid);
                const Real line_absorption =  opacity (p, lid) * lineProducingSpecies[l].J_lin(p, k);

                temp_cooling_rate+=line_emission-line_absorption;
            }
        }

        cooling_rate[p]=temp_cooling_rate;
    });
}

/// Computes the cooling rates using the discretized flux definition approach
///   @param[in] model: the model for which to compute the cooling rates
/// Warning: cannot compute LTE cooling rates
inline void Cooling :: compute_cooling_grad_I(Model& model)
{
    //FIXME: For now the gradient has not yet been computed. This still has to be done!
    threaded_for(p, model.parameters->npoints(),
    {
        Real temp_cooling_rate=0.0;
        //err, we use the mean line intensity instead
        for (Size l = 0; l < model.parameters->nlspecs(); l++)
        {
            for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
            {
                const Size lid = line_index (l, k);

                const Real line_cooling =  opacity (p, lid) * line_grad_I(TODO);

                temp_cooling_rate+=line_cooling;
            }
        }

        cooling_rate[p]=temp_cooling_rate;
    });

}



// ///  OLD IMPLEMENTATION (use as reference for collisional approach)
// ///  Calculates the cooling rates
// ///  NOTE: call this only after computing the radiation field
// ///////////////////////////////////////////
// inline void Lines :: calculate_cooling_rates (const Double2 &abundance, const Vector<Real> &temperature)
// {
//     // Matrix<Real> partial_cooling_rates; //for each point and line species
//     // partial_cooling_rates.resize (parameters.npoints(),  parameters.nlspecs());
//     //FIXME: change the colpar.adjust_abundance... noncense to finally allow parallelisation ;-; same for interpolate_collision_coefficients
//     // threaded_for (p, parameters.npoints(),
//     for (Size p=0; p<parameters.npoints(); p++)
//     {
//         Real cooling_of_point=0;
//
//         const Real tmp = temperature[p];
//
//         // TODO: if certain species may not participate, exclude them here
//         for (Size l = 0; l < parameters.nlspecs(); l++)
//         {
//             Real cooling_of_lspec=0;
//
//             for (CollisionPartner &colpar : lineProducingSpecies[l].linedata.colpar)
//             {
//                 Real cooling_with_colpar=0;
//
//                 Real abn = abundance[p][colpar.num_col_partner];
//
//                 colpar.adjust_abundance_for_ortho_or_para (tmp, abn);
//                 colpar.interpolate_collision_coefficients (tmp);
//
//                 for (Size k = 0; k < colpar.ncol; k++)
//                 {
//                     const Real collision_rate_ij = colpar.Cd_intpld[k] * abn; //deexcitation->heating
//                     const Real collision_rate_ji = colpar.Ce_intpld[k] * abn; //excitation  ->cooling
//
//                     const Size I = lineProducingSpecies[l].index (p, colpar.icol[k]);
//                     const Size J = lineProducingSpecies[l].index (p, colpar.jcol[k]);
//
//                     const Real energy_diff=lineProducingSpecies[l].linedata.energy[colpar.icol[k]]-lineProducingSpecies[l].linedata.energy[colpar.jcol[k]];
//                     // std::cout<<"energy_diff: "<<energy_diff<<std::endl;
//
//                     cooling_with_colpar-=lineProducingSpecies[l].population[I]*collision_rate_ij*energy_diff;
//                     cooling_with_colpar+=lineProducingSpecies[l].population[J]*collision_rate_ji*energy_diff;
//
//                 }
//
//                 cooling_of_lspec+=cooling_with_colpar;
//             }
//
//             // for (Size k = 0; k < lineProducingSpecies[l].linedata.nrad; k++)
//             // {
//             //     const Size lid = line_index (l, k);
//             //
//             //     //TODO: replace this with deltaE(=hv)*Collisional rates()(*angle?)
//             //     //This calculates 'cooling' by computing what energy escapes by radiation (not by ckecking what energy is extracted from the gas (towards quantum states))
//             //     //This is can be wrong and should be replaced
//             //     cooling_of_lspec+=emissivity(p, lid)*line[lid]*FOUR_PI;//emissivity*freq(->energy)*4pi(->integrate out direction)
//             //     cooling_of_lspec-=opacity   (p, lid)*line[lid]*FOUR_PI*lineProducingSpecies[l].Jlin[p][k];//also times averaged intensity
//             // }
//
//             cooling_of_point+=cooling_of_lspec;
//             //TODO: figure out whether we need the individual cooling rates of the line producing species
//         }
//         // //TODO: add fancy reduction operation instead of this
//         cooling_rates[p]=cooling_of_point;
//     }//)
//
// }
