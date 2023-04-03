#ifndef FILM_HPP_INCLUDED
#define FILM_HPP_INCLUDED

#include <concepts>
#include <type_traits>

#include "../lattices/3d/3d.hpp"
#include "../lattices/3d/fcc.hpp"
#include "../lattices/base_lattice.hpp"
#include "../lattices/borders_conditions.hpp"

namespace qss::inline nanostructures
{
template <typename lattice_t>
concept ThreeD_Lattice = std::same_as<typename lattice_t::sizes_t, qss::lattices::three_d::sizes_t>;

template <ThreeD_Lattice lattice_t> struct film final : public lattice_t
{
    using typename lattice_t::coords_t;
    using typename lattice_t::sizes_t;
    using xy_border_condition = typename qss::borders_conditions::periodic<typename lattice_t::coords_t::size_type,
                                                                           typename lattice_t::sizes_t::size_type>;
    using z_border_condition = typename qss::borders_conditions::sharp<typename lattice_t::coords_t::size_type,
                                                                       typename lattice_t::sizes_t::size_type>;
    double J;

    film(const lattice_t &lattice, double J_ = 1.0) : lattice_t{lattice}, J{J_} {};
    film(lattice_t &&lattice, double J_ = 1.0) noexcept : lattice_t{std::move(lattice)}, J{J_} {};

    typename lattice_t::sizes_t::size_type get_last_z(const typename qss::lattices::three_d::fcc_coords_t &coord) const
    {
        using cast_t = typename lattice_t::sizes_t::size_type;
        switch (coord.w)
        {
        case 0:
            return static_cast<cast_t>(this->sizes.z / 2 + this->sizes.z % 2 - 1);
            break;
        case 1:
            return static_cast<cast_t>(this->sizes.z / 2 + this->sizes.z % 2 - 1);
            break;
        case 2:
            return static_cast<cast_t>(this->sizes.z / 2 - 1);
            break;
        case 3:
            return static_cast<cast_t>(this->sizes.z / 2 - 1);
            break;
        default:
            throw std::out_of_range("coord.w out of range : " + std::to_string(coord.w));
            break;
        }
    }

    // returns sum of replased elems and their amount
    std::pair<typename lattice_t::value_t::magn_t, unsigned int> fill_plane(unsigned int z,
                                                                            const typename lattice_t::value_t &value)
    {
        const auto patterns = qss::lattices::three_d::get_plane_XY(z);
        auto amount = 0u;
        typename lattice_t::value_t::magn_t return_value{};
        for (const auto &pattern : patterns)
        {
            auto coord = pattern.get_coord();
            for (typename lattice_t::coords_t::size_type i{0}; i < this->sizes.x; ++i)
            {
                coord.x = i;
                for (typename lattice_t::coords_t::size_type j{0}; j < this->sizes.y; ++j)
                {
                    coord.x = j;
                    const auto exists =
                        qss::borders_conditions::use_border_conditions<z_border_condition, z_border_condition,
                                                                       z_border_condition>(coord, this->sizes);
                    if (exists)
                    {
                        return_value += this->get(coord);
                        amount++;
                        this->set(value, coord);
                    }
                }
            }
        }
        return {return_value, amount};
    }
};

template <typename spin_t, typename old_spin_t, template <typename = old_spin_t> class lattice_t>
requires ThreeD_Lattice<lattice_t<old_spin_t>>
[[nodiscard]] constexpr film<lattice_t<spin_t>> copy_structure(const film<lattice_t<old_spin_t>> &original) noexcept
{
    if constexpr (std::is_same_v<old_spin_t, spin_t>)
    {
        return original;
    }

    return film<lattice_t<spin_t>>{copy_structure<spin_t>(dynamic_cast<const lattice_t<old_spin_t> &>(original)),
                                   original.J};
}

template <ThreeD_Lattice lattice_t>
requires std::is_same_v<typename lattice_t::coords_t, qss::lattices::three_d::fcc_coords_t>
[[nodiscard]] std::optional<qss::lattices::three_d::fcc_coords_t> get_closest_neigbour_from_upper_film(
    const film<lattice_t> &other, const qss::lattices::three_d::fcc_coords_t &coord) noexcept
{
    auto result = coord;
    if (coord.z > 0 || coord.z < other.sublattices_sizes[coord.w].z - 1)
    {
        return {};
    }
    result.z = 0;
    return result;
}
template <ThreeD_Lattice lattice_t>
requires std::is_same_v<typename lattice_t::coords_t, qss::lattices::three_d::fcc_coords_t>
[[nodiscard]] std::optional<qss::lattices::three_d::fcc_coords_t> get_closest_neigbour_from_lower_film(
    const film<lattice_t> &other, const qss::lattices::three_d::fcc_coords_t &coord) noexcept
{
    auto result = coord;
    if (coord.z > 0 || coord.z < other.sublattices_sizes[coord.w].z - 1)
    {
        return {};
    }
    result.z = other.sublattices_sizes[coord.w].z - 1;
    return result;
}
} // namespace qss::inline nanostructures

#endif