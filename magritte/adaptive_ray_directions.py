import healpy as hp
import numpy as np
import concurrent.futures

class AdaptiveRayDirectionHelper:
    def __init__(self, positions: np.ndarray, Ntop: int = 2, Nrefinements: int = 4, Ncomparisons: int = 4000) -> None:
        """Initializes the AdaptiveRayDirections class

        Args:
            positions (np.ndarray): 3D positions of the points (used for guessing the most important directions)
            Ntop (int, optional): Refinement parameter. Determines how much can be refined in each level. Should correspond to at least the number of interesting regions in the model. Defaults to 2.
            Nrefinements (int, optional): Refinement parameter. Determines the amount of refinement levels. Defaults to 4.
            Ncomparisons (int, optional): Randomly sample Ncomparison positions for determine the adaptive ray directions. Required computation time scales linearly. Defaults to 4000.
        """
        self.Nrefinements: int = Nrefinements
        self.Nside: int = 2**Nrefinements
        self.Npix: int = 12*self.Nside**2
        self.halfNpix: int = self.Npix//2
        self.Ntop: int = Ntop
        self.Ncomparisons: int = Ncomparisons
        self.N_adaptive_angles: int = self.get_number_adaptive_directions()
        self.half_N_adaptive_angles: int = self.N_adaptive_angles//2
        self.positions: np.ndarray[None, 3] = positions
        if Ncomparisons>=positions.shape[0]:
            self.refpositions = positions
        else:
            self.refpositions = self.positions[np.random.choice(self.positions.shape[0], self.Ncomparisons)]

    def get_number_adaptive_directions(self) -> int:
        """Returns the number of adaptive ray directions

        Returns:
            int: Number of adaptive ray directions
        """
        sum_angles: int = 12
        for i in range(self.Nrefinements):
            sum_angles += 6*min(self.Ntop, 6*4**(i))

        return sum_angles
    
    def calc_pixcount(self, posidx: int) -> np.ndarray:
        """Calculates the number of reference points in each pixel

        Args:
            posidx (int): Position index

        Returns:
            np.ndarray[self.Npix]: Number of counts per pixel
        """
        diffpos = self.refpositions - self.positions[posidx, :]
        pix = hp.vec2pix(self.Nside, diffpos[:,0], diffpos[:,1], diffpos[:,2])
        pixcount = np.bincount(pix, minlength = self.Npix)

        return pixcount
    
    def get_adaptive_direction_weight(self, pixcount):
        #upper bounds for array size; will need to be pruned later TODO: get exact size precomputed; all positions must have same size, due to the number of available directions being the same
        #Magritte uses symmetric ray directions; so I just added the oppisite bins at the start
        #To get the antipodal directions, just grab first half of healpix indices; then let healpix compute the corresponding -dir indices
        best_weights = np.zeros(self.N_adaptive_angles)
        best_directions = np.zeros((self.N_adaptive_angles, 3))
        
        #use nested ordering, for easier referinement operations
        full_nested_pixel_ordering = hp.reorder(pixcount, r2n=True)
        # Magritte uses antipodal rays, so compute them; unfortunately healpy api returns tuples for some reason... (-> vstack)
        antipodal_indices_nest = hp.vec2pix(self.Nside, *-np.vstack(hp.pix2vec(self.Nside, np.arange(self.halfNpix), nest=True)), nest=True)
        #algorithm mainly intended for improving ray discretization on outside of model, so use the maximum instead of the mean.
        nested_pixel_ordering = np.maximum(full_nested_pixel_ordering[:self.halfNpix], full_nested_pixel_ordering[antipodal_indices_nest])


        #TODO: think about which kind of average to use for the coarser level
        coarser_values = nested_pixel_ordering[::4] + nested_pixel_ordering[1::4] + nested_pixel_ordering[2::4] + nested_pixel_ordering[3::4]
        levels_to_pick = min(self.Ntop, self.halfNpix//4)#on the coarser grid (half the rays, 1/4 for coarsening)
        highest_coarser_indices = np.argpartition(coarser_values, -levels_to_pick)[-levels_to_pick:]
        already_picked_indices = 0
        
        if levels_to_pick>0:
            best_directions[:4*levels_to_pick, :] = np.vstack(hp.pix2vec(self.Nside, 4*np.repeat(highest_coarser_indices, 4)+np.tile(np.arange(4), levels_to_pick), nest=True)).T
        
            # All directions at a given level have equal weight
            best_weights[:4*levels_to_pick] = 1.0/self.Npix#*np.ones(4*levels_to_pick)

        already_picked_indices += 4*levels_to_pick

        #Limit the amount of iterations, as hp.pix2vec doesnt handle empty arrays.
        max_iter = self.Nrefinements - 1 - max(0, int(np.floor(np.log2(self.Ntop/6)/2)))

        for i in range(max_iter):
            # For consistency, the discretizations at a finer level must also have some corresponding discretization at a higher level
            already_planned_indices = np.unique(highest_coarser_indices//4)
            len_plan_indices = len(already_planned_indices)
            # Amount of discretizations to pick is limited by the number of already planned indices
            levels_to_pick = max(0, min(self.Ntop, self.halfNpix//(4**(i+2)))-len_plan_indices)
            proposed_indices = 4*np.repeat(already_planned_indices, 4)+np.tile(np.arange(4), len(already_planned_indices))
            proposed_indices = np.array(list(set(proposed_indices) - set(highest_coarser_indices)))#is on the current level
            len_prop_indices = len(proposed_indices)
            
            # Add the already planned directions (if nonzero amount)
            if len_prop_indices>0:
                best_directions[already_picked_indices:already_picked_indices+len_prop_indices, :] = np.vstack(hp.pix2vec(2**(self.Nrefinements-i-1), proposed_indices, nest=True)).T
                best_weights[already_picked_indices:already_picked_indices+len_prop_indices] = 4**(i+1)/self.Npix#*np.ones(len_prop_indices)
            
            already_picked_indices += len_prop_indices
            
            # And calculate statistics for the coarser level, in order to pick an additional levels_to_pick discretized directions
            coarser_values = coarser_values[::4] + coarser_values[1::4] + coarser_values[2::4] + coarser_values[3::4]
            # But make sure we never pick the indices we have already taken
            coarser_values[already_planned_indices] = -1
            # Hack: either use n (nonzero) last elements, or if zero, explicitly say no elements
            highest_coarser_indices = np.argpartition(coarser_values, -levels_to_pick)[-levels_to_pick or len(coarser_values):]
            
            # And now add these chosen directions
            best_directions[already_picked_indices:already_picked_indices+4*levels_to_pick, :] = np.vstack(hp.pix2vec(2**(self.Nrefinements-i-1), 4*np.repeat(highest_coarser_indices, 4)+np.tile(np.arange(4), levels_to_pick), nest=True)).T
            
            best_weights[already_picked_indices:already_picked_indices+4*levels_to_pick] = 4**(i+1)/self.Npix#*np.ones(4*levels_to_pick)
            already_picked_indices += 4*levels_to_pick
            
            #and add these two coarser level things together, to keep track of everything
            highest_coarser_indices = np.concatenate((highest_coarser_indices, already_planned_indices))
            
        remaining_indices = np.array(list(set(np.arange(6)) - set(highest_coarser_indices)))
        len_remaining = len(remaining_indices)
        
        # If clause gets triggered when topN < 6 (as then not all coarsest levels are already refined)
        if (len(remaining_indices)>0):
            best_directions[already_picked_indices:already_picked_indices+len_remaining, :] = np.vstack(hp.pix2vec(1, remaining_indices, nest=True)).T
            best_weights[already_picked_indices:already_picked_indices+len_remaining] = 1.0/12.0
        
        
        #Compute the antipodal directions, using the same weight
        best_directions[self.half_N_adaptive_angles:2*self.half_N_adaptive_angles, :] = - best_directions[:self.half_N_adaptive_angles]
        best_weights[self.half_N_adaptive_angles:2*self.half_N_adaptive_angles] = best_weights[:self.half_N_adaptive_angles]
        #print("final weight", best_weights, np.sum(best_weights), np.sum(best_weights>0.0))
        #print("final directions", best_directions)
        
        return best_directions, best_weights
    

    def get_adaptive_directions_block(self, start: int, end: int) -> 'tuple[np.ndarray, np.ndarray]':
        """Computes the adaptive ray directions and corresponding weights for a block of positions on a single thread

        Args:
            start (int): Start index of the block
            end (int): End index of the block

        Returns:
            tuple[np.ndarray, np.ndarray]: Tuple containing the adaptive ray directions and corresponding weights
        """
        best_directions = np.zeros((end-start, self.N_adaptive_angles, 3))
        best_weights = np.zeros((end-start, self.N_adaptive_angles))

        for i in range(start, end):
            if i%1000 == 0:
                print(i)
            pixcount = self.calc_pixcount(i)
            best_directions[i-start, :, :], best_weights[i-start, :] = self.get_adaptive_direction_weight(pixcount)

        return best_directions, best_weights
    

    def get_adaptive_directions(self, points_per_thread = 1000) -> 'tuple[np.ndarray, np.ndarray]':
        """Computes the adaptive ray directions and corresponding weights

        Returns:
            tuple[np.ndarray, np.ndarray]: Tuple containing the adaptive ray directions and corresponding weights
        """
        best_directions = np.zeros((self.positions.shape[0], self.N_adaptive_angles, 3))
        best_weights = np.zeros((self.positions.shape[0], self.N_adaptive_angles))

        # For each points_per_thread points, calculate the adaptive ray directions and weights
        index_starts = np.arange(0, self.positions.shape[0], points_per_thread)
        index_ends = np.append(index_starts[1:], self.positions.shape[0])

        print("Computing adaptive rays for points:")

        with concurrent.futures.ThreadPoolExecutor() as executor:
            adaptive_directions_weights = executor.map(self.get_adaptive_directions_block, index_starts, index_ends)
            for i, (directions, weights) in enumerate(adaptive_directions_weights):
                best_directions[index_starts[i]:index_ends[i], :, :] = directions
                best_weights[index_starts[i]:index_ends[i], :] = weights

        return best_directions, best_weights
    
    def get_antipod_indices(self) -> np.ndarray:
        """Returns the antipod indices, which are the same for each position

        Returns:
            np.ndarray: antipod indices; shape: (self.N_adaptive_angles,)
        """
        ind = np.arange(self.half_N_adaptive_angles)
        return np.concatenate((ind + self.half_N_adaptive_angles, ind))
