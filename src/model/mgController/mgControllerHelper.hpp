#pragma once


#include <memory>
#include "model/mgController/mgController.hpp"
#include "model/mgController/naiveMG/naiveMG.hpp"
#include "model/mgController/vCycle/vCycle.hpp"
#include "model/mgController/wCycle/wCycle.hpp"



//The Multigrid controller helper structure
//just for the practical implementation
// for getting where we currently are in the multigrid sequence
struct MgControllerHelper : virtual public MgController
{
  private:

    std::shared_ptr<MgController> implementation_instance_ptr=std::shared_ptr<MgController>(nullptr);

  public:

    //Inherited from mgController
    // enum class Actions {
    //   interpolate_levelpops,  //interpolating the levelpopulations
    //   interpolate_corrections,//interpolating the relative corrections
    //   restrict,               //restricting the levelpops and residuals
    //   stay,                   //remain on the same grid level; keep iterating
    //   finish,                 //the entire multigrid procedure is finished
    //   goto_coarsest           //signals that we start iterating at the coarsest grid//also possible after reset()
    // };

    //DO NOT USE THE DEFAULT IMPLEMENTATION!!!
    //FIXME: THROW ERROR WHEN TRYING TO USE DEFAULT IMPLEMENTATION
    MgControllerHelper()=default;

    //initializes the mgController
    //MgController(Size n_levels, Size n_pre_interpolation_steps);//TODO add much more
    // inline MgControllerHelper(MgController* implementation_instance_ptr);
    inline MgControllerHelper(std::shared_ptr<MgController> implementation_instance_ptr);

    //Because we cannot change constructor names, we use the default constructor just need to 'construct' using some regular functions
    // inline void UseNaiveMG(Size n_levels, Size finest_lvl);

    //returns the next action and updates what to do next
    inline Actions get_next_action() override;

    //returns the current level
    inline Size get_current_level() override;

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    inline void converged_on_current_grid() override;

    //Call when the solution has directly converged on the current grid.
    // In this case, we no longer need to use this and the coarser grid
    //inline void disable_current_and_coarser_grids();

    //resets the mgController
    //inline void reset();
};

#include "mgControllerHelper.tpp"
