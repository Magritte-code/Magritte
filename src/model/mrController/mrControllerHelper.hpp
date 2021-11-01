#pragma once


#include <memory>
#include "model/mrController/mrController.hpp"
#include "model/mrController/naiveMG/naiveMG.hpp"
#include "model/mrController/vCycle/vCycle.hpp"
#include "model/mrController/wCycle/wCycle.hpp"



//The multiresolution controller helper structure
//just for the practical implementation
//  for getting where we currently are in the multiresolution sequence
struct MrControllerHelper : virtual public MrController
{
  private:

    std::shared_ptr<MrController> implementation_instance_ptr=std::shared_ptr<MrController>(nullptr);

  public:

    //Inherited from mrController
    // enum class Actions {
    //   interpolate_levelpops,  //interpolating the levelpopulations
    //   interpolate_corrections,//interpolating the relative corrections
    //   restrict,               //restricting the levelpops and residuals
    //   stay,                   //remain on the same grid level; keep iterating
    //   finish,                 //the entire multiresolution procedure is finished
    //   goto_coarsest           //signals that we start iterating at the coarsest grid//also possible after reset()
    // };

    //DO NOT USE THE DEFAULT IMPLEMENTATION!!!
    //FIXME: THROW ERROR WHEN TRYING TO USE DEFAULT IMPLEMENTATION
    MrControllerHelper()=default;

    //initializes the mrController
    //mrController(Size n_levels, Size n_pre_interpolation_steps);//TODO add much more
    // inline mrControllerHelper(mrController* implementation_instance_ptr);
    inline MrControllerHelper(std::shared_ptr<MrController> implementation_instance_ptr);

    //returns the next action and updates what to do next
    inline Actions get_next_action() override;

    //returns the current level
    inline Size get_current_level() override;

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    inline void converged_on_current_grid() override;

    //resets the mrController TODO: maybe implment this (or just create a new one instead)
    //inline void reset();
};

#include "mrControllerHelper.tpp"
