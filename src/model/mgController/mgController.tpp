// Initializes the multigrid controller
MgController(Size nb_levels, Size nb_pre_interpolation_steps)
{
  //FIXME: check if nb_levels>=2 // otherwise we only have a single grid
  //FIXME: check if nb_pre_interpolation_steps>=1
  this.max_level=nb_levels-1;
  this.current_level=nb_levels-1;
  //this.min_level_visited=nb_levels-1;

  this.going_coarser=false;

  this.nb_pre_interpolation_steps=nb_pre_interpolation_steps;

}


//
//Currently implements a V-cycle with some initial steps on the coarsest grid
inline Actions MgController::next_action()
{

  if (is_next_action_set)
  {
    is_next_action_set=false;
    return next_action;
  }else{
    if (nb_pre_interpolation_steps==0)
    {//then change grid, depending on the direction
      if (current_level==0){
        going_coarser=true;
      } else if (current_level>=max_level) {
        going_coarser=false;
      }
      if (max_level==0){//there is only a single level left, so staying is the only option (except from finishing)
        return Actions::stay;
      }

      if (going_coarser){
        current_level++;
        return Actions::restrict;
      }
      else{
        current_level--;
        return Actions::interpolate_levelpops;
        //TODO:also interpolate corrections//add something to see whether it is the first time interpolating to that grid
      }

        //choose between interpolate_levelpops and interpolate_corrections

    }else{//otherwise just remain on the current grid
      nb_pre_interpolation_steps--;
      return Actions::stay;
    }
  }




}

//Call this when the solution is converged on the current grid.
//Sets the state such that the next action will be something else than stay (skips the following 'stay's)
inline void MgController::converged_on_current_grid()
{
  if (current_level==0)
  {
    next_action=finish;
    is_next_action_set=true;
    return;
  }
  if (going_coarser)
  {//then the iterations on this grids aren't doing much anymore, so disable the use of the current grid
    if (current_level>0)
    {
      next_action=interpolate_corrections;
      is_next_action_set=true;
      going_coarser=false;
      this.max_level=current_level-1;
      return;
    }
  }else{
    //converged on current grid, so no need for further iterations on it
    nb_pre_interpolation_steps=0;
  }
  return;
}
