#ifndef TREE_D_HPP_INCLUDED
#define TREE_D_HPP_INCLUDED

#include <cstdint>

namespace qss::lattices::three_d
{

    struct sizes_t
    {
        using size_type = unsigned short;
        size_type x = 0;
        size_type y = 0;
        size_type z = 0;
    };
}

#endif