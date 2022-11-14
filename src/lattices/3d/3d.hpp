#ifndef TREE_D_HPP_INCLUDED
#define TREE_D_HPP_INCLUDED

#include <cstdint>

namespace qss::lattices::three_d
{
    template <typename size_t = std::uint8_t>
    struct sizes_t
    {
        size_t x = 0;
        size_t y = 0;
        size_t z = 0;
    };
}

#endif