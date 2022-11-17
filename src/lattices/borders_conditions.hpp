#ifndef BORDERS_CONDITIONS_HPP_INCLUDED
#define BORDERS_CONDITIONS_HPP_INCLUDED

#include <cstdint>
#include <type_traits>
#include <iterator>
#include <optional>

#include "2d/2d.hpp"
#include "2d/square.hpp"
#include "3d/3d.hpp"

namespace qss::lattices
{
    namespace two_d
    {
        template <typename coord_size_t = std::int16_t>
        std::optional<square_coords_t<coord_size_t>> use_periodic_conditions(square_coords_t<coord_size_t> coord,
                                                                             const sizes_t<> &sizes)
        {
            if (coord.x >= sizes.x)
            {
                coord.x -= sizes.x;
            }
            if (coord.y >= sizes.y)
            {
                coord.y -= sizes.y;
            }
            if (coord.x < 0)
            {
                coord.x += sizes.x;
            }
            if (coord.y < 0)
            {
                coord.y += sizes.y;
            }
            return coord;
        }
        template <typename coord_size_t = std::int16_t>
        std::optional<square_coords_t<coord_size_t>> use_sharp_conditions(square_coords_t<coord_size_t> coord,
                                                                          const sizes_t<> &sizes)
        {
            if (coord.x >= sizes.x || coord.x < 0 || coord.y >= sizes.y || coord.y < 0)
            {
                return {};
            }
            else
            {
                return coord;
            }
        }
    }
}

#endif