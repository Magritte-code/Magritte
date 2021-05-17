#pragma once


// Initializes the multigrid controller
// Finest level is for debug purposes, restrict the finest level until we iterate
inline VCycle::VCycle(Size n_levels, Size finest_lvl, Size n_pre_interpolation_steps, Size max_n_iterations)
{
  //FIXME: check if n_levels>=2 // otherwise we only have a single grid
  //FIXME: check if n_pre_interpolation_steps>=1
  max_level=n_levels-1;
  current_level=n_levels-1;
  //this.min_level_visited=n_levels-1;
  this->finest_lvl=finest_lvl;

  this->max_n_iterations=max_n_iterations;

  this->going_coarser=true;
  this->n_pre_interpolation_steps=n_pre_interpolation_steps;

}


inline Size VCycle::get_current_level()
{
  return current_level;
}


//Currently implements a V-cycle with some initial steps on the coarsest grid
inline MgController::Actions VCycle::get_next_action()
{
  //If the max number of iterations is reached, stop//TODO: think about whether this part is still useful
  if ((current_n_iterations>=max_n_iterations)&&(!not_yet_iterated))
  {
    return Actions::finish;
  }

  //If converged, then use whatever recommended by converged_on_current_grid
  if (is_next_action_set)
  {
    is_next_action_set=false;
    return next_action;
  }else{
    if (n_pre_interpolation_steps==0)//if the initial pre-interpolation steps on the coarsest grid are already done
    {//then change grid, depending on the direction
      if (not_yet_iterated)
      {
        not_yet_iterated=false;
        return Actions::stay;
      }
      if (current_level<=finest_lvl){
        if (first_upward)
        {first_upward=false;}
        else//increasing the counter of number iterations every time the finest grid is reached (except the first time)
        {
          current_n_iterations++;
          if (current_n_iterations>=max_n_iterations)//If the max number of iterations is reached, stop
          {
            return Actions::finish;
          }
        }
        going_coarser=true;
      } else if (current_level>=max_level) {
        going_coarser=false;
      }
      if (max_level==finest_lvl){//there is only a single level left, so staying is the only option (except from finishing)
        return Actions::stay;
      }
      if (going_coarser){
        current_level++;
        not_yet_iterated=true;
        return Actions::restrict;
      }
      else{
        current_level--;
        not_yet_iterated=true;
        if (first_upward)
        {
          return Actions::interpolate_levelpops;
        }else{//if it is not the first time upward, interpolate the corrections instead
          return Actions::interpolate_corrections;
        }
      }
    }
    else
    {//still need to do some pre_interpolation steps
      not_yet_iterated=false;
      n_pre_interpolation_steps--;
      return Actions::stay;
    }
  }

}


//Call this when the solution is converged on the current grid.
//Sets the state such that the next action will be something else than stay (skips the following 'stay's)
inline void VCycle::converged_on_current_grid()
{
  n_pre_interpolation_steps=0;
  if (current_level==finest_lvl)
  {
    next_action=Actions::finish;
    is_next_action_set=true;
    return;
  }
  else if (current_level>finest_lvl)
  {
    if (first_upward)//note:convergence during the first time going upward should not be practically possible, but hey, it could happen...
    {
      next_action=Actions::interpolate_levelpops;
    }else{
      next_action=Actions::interpolate_corrections;
    }
    is_next_action_set=true;
    current_level--;
    max_level=current_level;
    return;
  }
  else
  {
    next_action=Actions::finish;
    //TODO also return error
    std::cout<<"Somehow, you are currently on a level finer than the finest level you allowed"<<std::endl;
    std::cout<<"Finishing the multigrid computation either way"<<std::endl;
    std::cout<<"Finest level: "<<finest_lvl<<"Current level: "<<current_level<<std::endl;
    return;
  }
}
