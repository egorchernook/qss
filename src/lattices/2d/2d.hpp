#ifndef TWO_D_HPP_INCLUDED
#define TWO_D_HPP_INCLUDED

#include <cstdint>
#include <utility>

namespace qss::lattices::two_d
{
    template <typename size_t = std::uint8_t>
    struct sizes_t
    {
        size_t x = 0;
        size_t y = 0;
    };

}

#endif