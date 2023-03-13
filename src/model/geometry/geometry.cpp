#include "geometry.hpp"
#include "tools/heapsort.hpp"

void Geometry :: read (const Io& io)
{
    points  .read (io);
    rays    .read (io);
    boundary.read (io);

    lengths.resize (parameters->hnrays(), parameters->npoints());

    //Setup only necessary for the comoving solvers, as we want a rough indicication of the size of the elements
    vector<Size> temp_sorted_position_indices;
    sorted_position_indices.resize(parameters->npoints());
    temp_sorted_position_indices.resize(parameters->npoints());
    vector<double> neg_mean_dist;
    neg_mean_dist.resize(parameters->npoints());
    //properly initialize sorted indices vector + mean distances
    //TODO use threaded for
    for (Size i=0; i<parameters->npoints(); i++)
  {
      temp_sorted_position_indices[i]=i;
        neg_mean_dist[i]=-get_mean_dist_to_point(i);
    }
    // I want to sort from largest mean dist to lowest, thus I sort based on - the distance
    heapsort(neg_mean_dist, temp_sorted_position_indices);
    // and copy the temporary vector to a Vector
    // std::cout<<"sorted position indices: "<<std::endl;
    for (Size i=0; i<parameters->npoints(); i++)
    {
       sorted_position_indices[i]=temp_sorted_position_indices[i];
       // std::cout<<temp_sorted_position_indices[i]<<std::endl;
    }
}


void Geometry :: write (const Io& io) const
{
    points  .write (io);
    rays    .write (io);
    boundary.write (io);
}
