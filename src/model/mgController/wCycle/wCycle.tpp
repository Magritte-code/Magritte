#pragma once

#include<cmath>

// Initializes the multigrid controller
// Finest level is for debug purposes, restrict the finest level until we iterate
inline WCycle::WCycle(Size nb_levels, Size finest_lvl, Size nb_pre_interpolation_steps, Size max_nb_iterations)
{
  //FIXME: check if nb_levels>=2 // otherwise we only have a single grid
  //FIXME: check if nb_pre_interpolation_steps>=1

  max_level=nb_levels-1;
  initialize_w_cycle(max_level-finest_lvl);
  current_level=nb_levels-1;
  //this.min_level_visited=nb_levels-1;
  this->finest_lvl=finest_lvl;

  this->max_nb_iterations=max_nb_iterations;

  this->going_coarser=true;
  this->nb_pre_interpolation_steps=nb_pre_interpolation_steps;

}


//Returns the current level //TODO: DEPRECATE THIS
inline Size WCycle::get_current_level()
{
  return current_level;
}


//Currently implements a W-cycle with some initial steps on the coarsest grid
inline MgController::Actions WCycle::get_next_action()
{
  //When done enough iterations, just stop
  if ((current_nb_iterations>=max_nb_iterations)&&(!not_yet_iterated))
  {
    return Actions::finish;
  }

  //If converged, then use whatever recommended by converged_on_current_grid
  if (is_next_action_set)
  {
    is_next_action_set=false;
    return next_action;
  }else{//otherwise alternate between staying on the grid and going to the next finer/coarser grid as determined by the action_order
    if (nb_pre_interpolation_steps==0)
    {//
      if (not_yet_iterated)
      {
        not_yet_iterated=false;
        return Actions::stay;
      }
      if (first_upward)
      {
        if (current_level>finest_lvl)
        {
          current_level--;
          not_yet_iterated=true;
          return Actions::interpolate_levelpops;
        }else{//the first time upward has ended, so setting flag
          first_upward=false;
        }
      }

      //Note: only three actions should ever be in the action_order: restrict and interpolate_corrections
      //Otherwise we will need to rewrite this
      Actions current_action=action_order[curr_action_nb];
      curr_action_nb++;
      //curr_action_nb should wrap around // this also denotes the end of an iteration
      if (curr_action_nb>=action_order.size())
      {
        curr_action_nb=0;
        current_nb_iterations++;
      }
      not_yet_iterated=true;

      //Note: this if-statement is just for the internal level variable
      if (current_action==Actions::interpolate_corrections)
      {current_level--;}
      else if (current_action==Actions::restrict)//restriction
      {current_level++;}
      else {;}//do nothing
      return current_action;

    }
    else
    {//still need to do some pre_interpolation steps
      not_yet_iterated=false;
      nb_pre_interpolation_steps--;
      return Actions::stay;
    }
  }




}


//Call this when the solution is converged on the current grid.
//Sets the state such that the next action will be something else than stay (skips the following 'stay's)
inline void WCycle::converged_on_current_grid()
{
  nb_pre_interpolation_steps=0;
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
    Size temp_nb_level_diff=current_level-finest_lvl;

    //initialize new w-cycle
    initialize_w_cycle(temp_nb_level_diff);
    //and set position correct (just after all the restrictions in the beginning)
    curr_action_nb=temp_nb_level_diff;

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

/// Initializes the w-cycle
//Note: does not contain logic about staying on a grid level (which should happen every other iteration anyway)
//TODO: OPTIMIZE THIS A BIT (if people are asking after 20+ multigrid levels) // NEVER MIND: that means that a single cycle would be more than a million actions long...
inline void WCycle::initialize_w_cycle(Size nb_level_diff)
{ //FIXME: throw warning when nb_level_diff=0
  if (nb_level_diff==0)
  {
    vector<MgController::Actions> temp_action_order{MgController::Actions::do_nothing};
    action_order=temp_action_order;
    return;
  }
  Size curr_nb_levels=1;
  //start with mini V-cycle
  vector<MgController::Actions> temp_action_order{MgController::Actions::restrict, MgController::Actions::interpolate_corrections};
  //calculate order by copying and adding a single action to the front and the back

  //calculating the total size to reserve (always '*2+2') always copy and then add to front and back
  Size size_to_reserve=std::pow(2,nb_level_diff+1)-2;
  //explicit formula above // for (Size i=0; i<nb_level_diff<i++) {size_to_reserve=size_to_reserve*2+2;}
  temp_action_order.reserve(size_to_reserve);

  //example of what this does; nb_level_diff = 1) \/  2) \\/\//  3) \\\/\//\\/\///
  // in which '\' respresents the restriction and '/' respresents the correction interpolation
  while (curr_nb_levels<nb_level_diff)
  {
    //std::back_inserter invalidates the iterator when the size of the vector is exceeded->Undefined behaviour
    // by reserving enough space, we ignore this issue
    std::copy(temp_action_order.begin(), temp_action_order.end(), std::back_inserter(temp_action_order));
    //an additional restrict in front
    temp_action_order.insert(temp_action_order.begin(), MgController::Actions::restrict);
    //and an additional interpolate_corrections to the back
    temp_action_order.push_back(MgController::Actions::interpolate_corrections);
    curr_nb_levels++;
  }
  action_order=temp_action_order;
}
