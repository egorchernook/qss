#ifndef FILM_HPP_INCLUDED
#define FILM_HPP_INCLUDED

#include <concepts>
#include <type_traits>

#include "../lattices/base_lattice.hpp"
#include "../lattices/3d/3d.hpp"
#include "../lattices/3d/fcc.hpp"
#include "../lattices/borders_conditions.hpp"

namespace qss::nanostructures
{
    template <typename lattice_t>
    concept ThreeD_Lattice = std::same_as<typename lattice_t::sizes_t, qss::lattices::three_d::sizes_t>;

    template <ThreeD_Lattice lattice_t>
    struct film final : public lattice_t
    {
        using xy_border_condition = typename qss::border_conditions::periodic<
            typename lattice_t::coords_t::size_type,
            typename lattice_t::sizes_t::size_type>;
        using z_border_condition = typename qss::border_conditions::sharp<
            typename lattice_t::coords_t::size_type,
            typename lattice_t::sizes_t::size_type>;
        double J;

        film(const lattice_t &lattice, double J_ = 1.0) : lattice_t{lattice}, J{J_} {};
        film(lattice_t &&lattice, double J_ = 1.0) noexcept : lattice_t{std::move(lattice)}, J{J_} {};
    };

    template <typename old_spin_t,
              template <typename = old_spin_t> class lattice_t,
              typename spin_t>
    requires ThreeD_Lattice<lattice_t<old_spin_t>>
    [[nodiscard]] constexpr film<lattice_t<spin_t>>
    copy_structure(const film<lattice_t<old_spin_t>> &original) noexcept
    {
        if constexpr (std::is_same_v<old_spin_t, spin_t>)
        {
            return original;
        }

        return film<lattice_t<spin_t>>{
            qss::lattices::copy_structure<spin_t>(original),
            original.J};
    }

    template <ThreeD_Lattice lattice_t>
    requires std::is_same_v<typename lattice_t::coords_t, qss::lattices::three_d::fcc_coords_t>
    [[nodiscard]] std::optional<qss::lattices::three_d::fcc_coords_t>
    get_closest_neigbour_from_upper_film(const film<lattice_t> &other,
                                         const qss::lattices::three_d::fcc_coords_t &coord) noexcept
    {
        auto result = coord;
        if (coord.z > 0 || coord.z < other.sublattices_sizes[coord.w].z - 1)
        {
            return {};
        }
        result.z = other.sublattices_sizes[coord.w].z - 1;
        return result;
    }
    template <ThreeD_Lattice lattice_t>
    requires std::is_same_v<typename lattice_t::coords_t, qss::lattices::three_d::fcc_coords_t>
    [[nodiscard]] std::optional<qss::lattices::three_d::fcc_coords_t>
    get_closest_neigbour_from_lower_film(const film<lattice_t> &other,
                                         const qss::lattices::three_d::fcc_coords_t &coord) noexcept
    {
        auto result = coord;
        if (coord.z > 0 || coord.z < other.sublattices_sizes[coord.w].z - 1)
        {
            return {};
        }
        result.z = 0;
        return result;
    }
}

#endif