#include <tuple>

//returns tuple of index of neighbor vector and how many neighbors exist
inline std::tuple<Size, Size> Points :: get_neighbors(Size pointid)
{
    return std::make_tuple(cum_n_neighbors[pointid], n_neighbors[pointid]);
}
