#pragma once

#include<cmath>

///  Initializes the multiresolution controller
///    @param[in] n_levels: the number of level in the multi-resolution procedure; this includes the coarsest level
///    @param[in] finest_level: the finest level we consider for the multi-resolution procedure; might be useful
///  when only wanting the solution on a coarser level (for lesser computation time)
///    @param[in] max_n_iterations: the maximum number of iterations allowed in the scheme
inline WCycle::WCycle(Size n_levels, Size finest_lvl, Size n_pre_interpolation_steps, Size max_n_iterations)
{
    //FIXME: check if n_levels>=2 // otherwise we only have a single grid
    //FIXME: check if n_pre_interpolation_steps>=1
    max_level=n_levels-1;
    initialize_w_cycle(max_level-finest_lvl);
    current_level=n_levels-1;

    this->finest_lvl=finest_lvl;

    this->max_n_iterations=max_n_iterations;

    this->going_coarser=true;
    this->n_pre_interpolation_steps=n_pre_interpolation_steps;

}


///  Returns the current level
inline Size WCycle::get_current_level()
{
    return current_level;
}


///  Returns the next action and updates what to do next
inline MrController::Actions WCycle::get_next_action()
{
    //When done enough iterations, just stop
    if ((current_n_iterations>=max_n_iterations)&&(!not_yet_iterated))
    {
        return Actions::finish;
    }

    //If converged, then use whatever recommended by converged_on_current_grid
    if (is_next_action_set)
    {
        is_next_action_set=false;
        return next_action;
    }
    else
    {//otherwise alternate between staying on the grid and going to the next finer/coarser grid as determined by the action_order
        if (n_pre_interpolation_steps==0)
        {
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

            //Note: only two actions should ever be in the action_order: restrict and interpolate_corrections
            //Otherwise we will need to rewrite this
            Actions current_action=action_order[curr_action_nb];
            curr_action_nb++;
            //curr_action_nb should wrap around // this also denotes the end of an iteration
            if (curr_action_nb>=action_order.size())
            {
                curr_action_nb=0;
                current_n_iterations++;
            }
            not_yet_iterated=true;

            //Note: this if-statement is just for the internal level variable
            if (current_action==Actions::interpolate_corrections)
            {
                current_level--;
            }
            else if (current_action==Actions::restrict)//restriction
            {
                current_level++;
            }
            else {;}//? this shouldn't happen ? -> do nothing
            return current_action;

        }
        else
        {//still need to do some pre_interpolation steps
            not_yet_iterated=false;
            n_pre_interpolation_steps--;
            return Actions::stay;
        }
    }
}


///  Call this when the solution is converged on the current grid.
///  Sets the state such that the next action will be something else than stay (skips the following 'stay's)
inline void WCycle::converged_on_current_grid()
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
        }
        else
        {
            next_action=Actions::interpolate_corrections;
        }
        is_next_action_set=true;
        current_level--;
        Size temp_n_level_diff=current_level-finest_lvl;

        //initialize new w-cycle
        initialize_w_cycle(temp_n_level_diff);
        //and set position correct (just after all the restrictions in the beginning)
        curr_action_nb=temp_n_level_diff;

        return;
    }
    else
    {
        next_action=Actions::finish;
        //TODO also return error
        std::cout<<"Somehow, you are currently on a level finer than the finest level you allowed"<<std::endl;
        std::cout<<"Finishing the multiresolution computation either way"<<std::endl;
        std::cout<<"Finest level: "<<finest_lvl<<"Current level: "<<current_level<<std::endl;
        return;
    }
}

/// Initializes the w-cycle
//Note: does not contain logic about staying on a grid level (which should happen every other iteration anyway)
inline void WCycle::initialize_w_cycle(Size n_level_diff)
{ //FIXME: throw warning when n_level_diff=0
    if (n_level_diff==0)
    {
        vector<MrController::Actions> temp_action_order{MrController::Actions::do_nothing};
        action_order=temp_action_order;
        return;
    }
    Size curr_n_levels=1;
    //start with mini V-cycle
    vector<MrController::Actions> temp_action_order{MrController::Actions::restrict, MrController::Actions::interpolate_corrections};
    //calculate order by copying and adding a single action to the front and the back

    //calculating the total size to reserve (always '*2+2') always copy and then add to front and back
    Size size_to_reserve=std::pow(2,n_level_diff+1)-2;
    //explicit formula above // for (Size i=0; i<n_level_diff<i++) {size_to_reserve=size_to_reserve*2+2;}
    temp_action_order.reserve(size_to_reserve);

    //example of what this does; n_level_diff = 1) \/  2) \\/\//  3) \\\/\//\\/\///
    // in which '\' respresents the restriction and '/' respresents the correction interpolation
    while (curr_n_levels<n_level_diff)
    {
        //std::back_inserter invalidates the iterator when the size of the vector is exceeded->Undefined behaviour
        // by reserving enough space, we ignore this issue
        std::copy(temp_action_order.begin(), temp_action_order.end(), std::back_inserter(temp_action_order));
        //an additional restrict in front
        temp_action_order.insert(temp_action_order.begin(), MrController::Actions::restrict);
        //and an additional interpolate_corrections to the back
        temp_action_order.push_back(MrController::Actions::interpolate_corrections);
        curr_n_levels++;
    }
    action_order=temp_action_order;
}
