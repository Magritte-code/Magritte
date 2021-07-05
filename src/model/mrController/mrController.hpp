#pragma once

//The abstract multiresolution controller structure
//  for getting where we currently are in the multiresolution sequence
//  For practical usage, use the mgControllerHelper class: that class is an implementation that implements all current multiresolution variants
struct MrController
{
  public:
    //All possible actions the multiresolution sequence is allowed to have
    enum class Actions {
      interpolate_levelpops,  //interpolating the levelpopulations
      interpolate_corrections,//interpolating the relative corrections
      restrict,               //restricting the levelpops and residuals
      stay,                   //remain on the same grid level; keep iterating
      finish,                 //the entire multiresolution procedure is finished
      goto_coarsest,          //signals that we start iterating at the coarsest grid//also possible after reset()
      do_nothing,             //signals that we do not need to do anything this iteration (for simplifying the state of some mgControllers)
      error                   //Somehow, we have an error
    };

    ///Returns the next action and updates what to do next
    virtual Actions get_next_action() = 0;

    ///Returns the current level
    virtual Size get_current_level() = 0;

    ///Call this when the solution is converged on the current grid.
    ///Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    virtual void converged_on_current_grid() = 0;

};
