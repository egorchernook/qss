#ifndef BORDERS_CONDITIONS_HPP_INCLUDED
#define BORDERS_CONDITIONS_HPP_INCLUDED

#include <cstdint>
#include <type_traits>
#include <iterator>
#include <optional>
#include <concepts>

#include "2d/2d.hpp"
#include "2d/square.hpp"
#include "3d/3d.hpp"
#include "3d/fcc.hpp"

namespace qss::borders_conditions
{
    template <typename coord_size_t, typename size_type>
        requires std::is_convertible_v<coord_size_t, int> &&
                 std::is_convertible_v<size_type, int>
    struct linear_border_conditions
    {
        using coord_size_type = coord_size_t;
        using sizes_size_type = size_type;
    };
    template <typename coord_size_t, typename size_type>
    struct periodic : linear_border_conditions<coord_size_t, size_type>
    {
        [[nodiscard]] std::optional<coord_size_t>
        operator()(coord_size_t coord, const size_type &size) const noexcept
        {
            if (coord >= size)
            {
                return {coord - size};
            }
            if (coord < 0)
            {
                return {coord + size};
            }
            return coord;
        }
    };
    template <typename coord_size_t, typename size_type>
    struct sharp : linear_border_conditions<coord_size_t, size_type>
    {
        [[nodiscard]] std::optional<coord_size_t>
        operator()(coord_size_t coord, const size_type &size) const noexcept
        {
            if (coord < 0 || coord >= size)
            {
                return {};
            }
            return coord;
        }
    };

    template <typename x_conds_t, typename y_conds_t>
        requires std::is_same_v<typename x_conds_t::coord_size_type,
                                typename y_conds_t::coord_size_type> &&
                 std::is_same_v<typename x_conds_t::sizes_size_type,
                                typename y_conds_t::sizes_size_type>
    [[nodiscard]] std::optional<typename qss::lattices::two_d::square_coords_t>
    use_border_conditions(const typename qss::lattices::two_d::square_coords_t &coord,
                          const typename qss::lattices::two_d::sizes_t &sizes) noexcept
    {
        using coords_t_size_type = typename qss::lattices::two_d::square_coords_t::size_type;
        static const x_conds_t x_conds{};
        static const y_conds_t y_conds{};
        const auto x = x_conds(coord.x, sizes.x);
        const auto y = y_conds(coord.y, sizes.y);
        if (!x.has_value() + !y.has_value())
        {
            return {};
        }
        return typename qss::lattices::two_d::square_coords_t{static_cast<coords_t_size_type>(x.value()),
                                                              static_cast<coords_t_size_type>(y.value())};
    }

    template <typename x_conds_t, typename y_conds_t, typename z_conds_t>
        requires std::is_same_v<typename x_conds_t::coord_size_type,
                                typename y_conds_t::coord_size_type> &&
                 std::is_same_v<typename y_conds_t::coord_size_type,
                                typename z_conds_t::coord_size_type> &&
                 std::is_same_v<typename x_conds_t::sizes_size_type,
                                typename y_conds_t::sizes_size_type> &&
                 std::is_same_v<typename y_conds_t::sizes_size_type,
                                typename z_conds_t::sizes_size_type>
    [[nodiscard]] std::optional<typename qss::lattices::three_d::fcc_coords_t>
    use_border_conditions(const typename qss::lattices::three_d::fcc_coords_t &coord,
                          const typename qss::lattices::three_d::sizes_t &sizes)
    {
        using coords_t_size_type = typename qss::lattices::three_d::fcc_coords_t::size_type;
        using size_type = typename x_conds_t::sizes_size_type;
        typename qss::lattices::three_d::sizes_t sublattice_size{};
        switch (coord.w)
        {
        case 0:
            sublattice_size = {static_cast<size_type>(sizes.x / 2 + sizes.x % 2),
                               static_cast<size_type>(sizes.y / 2 + sizes.y % 2),
                               static_cast<size_type>(sizes.z / 2 + sizes.z % 2)};
            break;
        case 1:
            sublattice_size = {static_cast<size_type>(sizes.x / 2),
                               static_cast<size_type>(sizes.y / 2),
                               static_cast<size_type>(sizes.z / 2 + sizes.z % 2)};
            break;
        case 2:
            sublattice_size = {static_cast<size_type>(sizes.x / 2 + sizes.x % 2),
                               static_cast<size_type>(sizes.y / 2),
                               static_cast<size_type>(sizes.z / 2)};
            break;
        case 3:
            sublattice_size = {static_cast<size_type>(sizes.x / 2),
                               static_cast<size_type>(sizes.y / 2 + sizes.y % 2),
                               static_cast<size_type>(sizes.z / 2)};
            break;
        default:
            throw std::out_of_range("coord.w out of range : " + std::to_string(coord.w));
            break;
        }

        static const x_conds_t x_conds{};
        static const y_conds_t y_conds{};
        static const z_conds_t z_conds{};
        const auto x = x_conds(coord.x, sublattice_size.x);
        const auto y = y_conds(coord.y, sublattice_size.y);
        const auto z = z_conds(coord.z, sublattice_size.z);
        if (!x.has_value() + !y.has_value() + !z.has_value())
        {
            return {};
        }
        else
        {
            return typename qss::lattices::three_d::fcc_coords_t{static_cast<std::uint8_t>(coord.w),
                                                                 static_cast<coords_t_size_type>(x.value()),
                                                                 static_cast<coords_t_size_type>(y.value()),
                                                                 static_cast<coords_t_size_type>(z.value())};
        }
    }
}

#endif