#pragma once

//The Multigrid controller structure
// for getting where we currently are in the multigrid sequence
struct MgController
{
  public:

    //TODO think a bit more about this
    enum class Actions {
      interpolate_levelpops,  //interpolating the levelpopulations
      interpolate_corrections,//interpolating the relative corrections
      restrict,               //restricting the levelpops and residuals
      stay,                   //remain on the same grid level; keep iterating
      finish,                 //the entire multigrid procedure is finished
      goto_coarsest,          //signals that we start iterating at the coarsest grid//also possible after reset()
      error                   //Somehow, we have an error
    };

    //initializes the mgController
    //MgController(Size nb_levels, Size nb_pre_interpolation_steps);//TODO add much more

    //returns the next action and updates what to do next
    virtual Actions get_next_action() = 0;

    //returns the current level
    virtual Size get_current_level() = 0;

    //Call this when the solution is converged on the current grid.
    //Sets the state such that the next action will be something else than stay (skips the following 'stay's)
    virtual void converged_on_current_grid() = 0;

    //Call when the solution has directly converged on the current grid.
    // In this case, we no longer need to use this and the coarser grid
    //inline void disable_current_and_coarser_grids();

    //resets the mgController
    //inline void reset();
};
