#ifndef BORDERS_CONDITIONS_HPP_INCLUDED
#define BORDERS_CONDITIONS_HPP_INCLUDED

#include <cstdint>
#include <type_traits>
#include <iterator>
#include <optional>

#include "2d/2d.hpp"
#include "2d/square.hpp"
#include "3d/3d.hpp"
#include "3d/fcc.hpp"

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

    namespace three_d
    {
        template <typename coord_size_t = std::int16_t>
        std::optional<fcc_coords_t<coord_size_t>> use_periodic_conditions(fcc_coords_t<coord_size_t> coord,
                                                                          const sizes_t<> &sizes)
        {
            sizes_t<coord_size_t> sublattice_size{};
            switch (coord.w)
            {
            case 0:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2 + sizes.x % 2),
                                   static_cast<coord_size_t>(sizes.y / 2 + sizes.y % 2),
                                   static_cast<coord_size_t>(sizes.z / 2 + sizes.z % 2)};
                break;
            case 1:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2),
                                   static_cast<coord_size_t>(sizes.y / 2),
                                   static_cast<coord_size_t>(sizes.z / 2 + sizes.z % 2)};
                break;
            case 2:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2 + sizes.x % 2),
                                   static_cast<coord_size_t>(sizes.y / 2),
                                   static_cast<coord_size_t>(sizes.z / 2)};
                break;
            case 3:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2),
                                   static_cast<coord_size_t>(sizes.y / 2 + sizes.y % 2),
                                   static_cast<coord_size_t>(sizes.z / 2)};
                break;
            default:
                throw std::out_of_range("coord.w out of range : " + std::to_string(coord.w));
                break;
            }

            if (coord.x >= sublattice_size.x)
            {
                coord.x -= sublattice_size.x;
            }
            if (coord.y >= sublattice_size.y)
            {
                coord.y -= sublattice_size.y;
            }
            if (coord.z >= sublattice_size.z)
            {
                coord.z -= sublattice_size.z;
            }
            if (coord.x < 0)
            {
                coord.x += sublattice_size.x;
            }
            if (coord.y < 0)
            {
                coord.y += sublattice_size.y;
            }
            if (coord.z < 0)
            {
                coord.z += sublattice_size.z;
            }
            return coord;
        }
        template <typename coord_size_t = std::int16_t>
        std::optional<fcc_coords_t<coord_size_t>> use_sharp_conditions(fcc_coords_t<coord_size_t> coord,
                                                                       const sizes_t<> &sizes)
        {
            sizes_t<coord_size_t> sublattice_size{};
            switch (coord.w)
            {
            case 0:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2 + sizes.x % 2),
                                   static_cast<coord_size_t>(sizes.y / 2 + sizes.y % 2),
                                   static_cast<coord_size_t>(sizes.z / 2 + sizes.z % 2)};
                break;
            case 1:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2),
                                   static_cast<coord_size_t>(sizes.y / 2),
                                   static_cast<coord_size_t>(sizes.z / 2 + sizes.z % 2)};
                break;
            case 2:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2 + sizes.x % 2),
                                   static_cast<coord_size_t>(sizes.y / 2),
                                   static_cast<coord_size_t>(sizes.z / 2)};
                break;
            case 3:
                sublattice_size = {static_cast<coord_size_t>(sizes.x / 2),
                                   static_cast<coord_size_t>(sizes.y / 2 + sizes.y % 2),
                                   static_cast<coord_size_t>(sizes.z / 2)};
                break;
            default:
                throw std::out_of_range("coord.w out of range : " + std::to_string(coord.w));
                break;
            }

            bool x_out_of_range = coord.x >= sublattice_size || coord.x < 0;
            bool y_out_of_range = coord.y >= sublattice_size || coord.y < 0;
            bool z_out_of_range = coord.z >= sublattice_size || coord.z < 0;
            if (x_out_of_range || y_out_of_range || z_out_of_range)
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