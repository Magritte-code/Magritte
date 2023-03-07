This series of notebooks accompagnies the Magritte 3 paper (in prep.), run using magritte 0.2.22 SHA: 8ab48f30dd6591d28e8c14f374d5a638c8d30cd9 .
Please note that it the notebooks have been cleaned and rerun later, so timings might be slightly off compared to the paper.

We expect the notebooks to be run in the following order:
- import_and_reduce_phantom
- compute_errors
- Plot_relative_differences

To show the speed of the remeshing algorithm, a Phantom model was remeshed and run using both GMSH and the recursive implementation.
This also includes a remesher based on Haar wavelets, but that implementation did not perform adequately, this did not make it into the paper.
As the recursive method is an improvement on the Haar wavelet remesher, we might deprecate the latter in the future.
