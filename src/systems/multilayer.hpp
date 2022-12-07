#ifndef MULTILAYER_HPP_INCLUDED
#define MULTILAYER_HPP_INCLUDED

#include <vector>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <concepts>
#include <algorithm>
#include <utility>
#include <type_traits>
#include <optional>
#include <algorithm>
#include <cstdint>

#include "../lattices/base_lattice.hpp"
#include "../lattices/3d/3d.hpp"
#include "../lattices/3d/fcc.hpp"
#include "../lattices/borders_conditions.hpp"
#include "../utility/functions.hpp"
#include "../random/random.hpp"
#include "../random/mersenne.hpp"
#include "film.hpp"
namespace qss::nanostructures
{

    template <ThreeD_Lattice lattice_t>
    struct multilayer_coords_t
    {
        using size_type = std::uint8_t;
        size_type idx;
        typename lattice_t::coords_t film_coord;
    };

    template <ThreeD_Lattice lattice_t>
    class multilayer final : private std::vector<film<lattice_t>>
    {
        using container_t = typename std::vector<film<lattice_t>>;
        std::vector<double> J_interlayers; // обменные интегралы взаимодействий плёнок.

        void check_sizes(const typename film<lattice_t>::sizes_t &value) const
        {
            if (value.x != this->at(0).sizes.x)
            {
                throw std::out_of_range("x size must be the same for all lattices : " + std::to_string(value.x) + " != " + std::to_string(this->at(0).sizes.x));
            }
            if (value.y != this->at(0).sizes.y)
            {
                throw std::out_of_range("y size must be the same for all lattices : " + std::to_string(value.y) + " != " + std::to_string(this->at(0).sizes.y));
            }
        }

    public:
        using film_t = film<lattice_t>;
        using coords_t = multilayer_coords_t<lattice_t>;
        [[nodiscard]] coords_t get(const coords_t &coord) const
        {
            if (coord.idx >= this->size())
            {
                throw std::out_of_range("coord.idx out of range : " + std::to_string(coord.idx) + " >= " + std::to_string(this->size()));
            }
            return this->at(coord.idx).get(coord.film_coord);
        }
        template <typename random_t = qss::random::mersenne::random_t<>>
        multilayer_coords_t<lattice_t> get_random_coord() const noexcept
        {
            using size_type = typename multilayer_coords_t<lattice_t>::size_type;
            static auto rand = random_t{qss::random::get_seed()};
            const size_type idx = static_cast<size_type>(rand(0, this->size()));
            const auto coord = this->at(idx).choose_random_node();
            return {idx, coord};
        }

        using container_t::begin;
        using container_t::cbegin;
        using container_t::cend;
        using container_t::crbegin;
        using container_t::crend;
        using container_t::end;
        using container_t::rbegin;
        using container_t::rend;
        using container_t::size;
        using container_t::operator[];
        using container_t::at;

        [[nodiscard]] constexpr multilayer(std::initializer_list<film_t> list,
                                           std::vector<double> &&J_interlayers_)
        {
            if (J_interlayers_.size() != list.size() - 1)
            {
                throw std::logic_error("number of interlayer exchange integrals must be number of films - 1 : " + std::to_string(J_interlayers.size()) + " != " + std::to_string(list.size() - 1));
            }
            this->reserve(list.size());
            for (auto &elem : list)
            {
                this->push_back(std::move(elem));
                check_sizes(this->back().sizes);
            }
            J_interlayers = std::move(J_interlayers_);
        }
        const std::vector<double> &get_J_interlayers() const noexcept
        {
            return J_interlayers;
        }

        [[nodiscard]] typename film_t::value_t::magn_t
        get_sum_of_closest_neighbours(const coords_t &central_) const noexcept
        {
            const auto idx = central_.idx;
            const auto central = central_.film_coord;
            const auto film = this->at(idx);
            bool has_upper = idx != size() - 1u;
            bool has_lower = idx != 0;

            using magn_t = typename film_t::value_t::magn_t;
            auto borders_conditions = [](const typename film_t::coords_t &coord_,
                                         const typename film_t::sizes_t &sizes)
            {
                return qss::border_conditions::use_border_conditions<
                    typename film_t::xy_border_condition,
                    typename film_t::xy_border_condition,
                    typename film_t::z_border_condition>(coord_, sizes);
            };
            auto neigs = get_closest_neighbours(central);
            const auto sizes = film.sizes;

            magn_t sum{};
            magn_t upper_neig_val{};
            magn_t lower_neig_val{};
            for (auto it = neigs.begin(); it != neigs.end(); ++it)
            {
                const auto coord = borders_conditions(*it, sizes);
                if (coord)
                {
                    sum += film.get(coord.value());
                }
                else
                {
                    if (has_upper)
                    {
                        const auto neig = get_closest_neigbour_from_upper_film(this->at(idx + 1u), central);
                        upper_neig_val += neig ? this->at(idx + 1u).get(neig.value()) : magn_t{};
                    }
                    if (has_lower)
                    {
                        const auto neig = get_closest_neigbour_from_lower_film(this->at(idx - 1u), central);
                        lower_neig_val += neig ? this->at(idx - 1u).get(neig.value()) : magn_t{};
                    }
                }
            }
            sum *= film.J;
            if (has_upper)
            {
                sum += J_interlayers[idx] * upper_neig_val;
            }
            if (has_lower)
            {
                sum += J_interlayers[idx - 1u] * lower_neig_val;
            }
            return sum;
        }
    };

    template <typename old_spin_t,
              template <typename = old_spin_t> class lattice_t,
              typename spin_t>
    requires ThreeD_Lattice<lattice_t<old_spin_t>>
    [[nodiscard]] constexpr multilayer<lattice_t<spin_t>>
    copy_structure(const multilayer<lattice_t<old_spin_t>> &original) noexcept
    {
        if constexpr (std::is_same_v<old_spin_t, spin_t>)
        {
            return original;
        }
        for (auto &film : original.get_films())
        {
        }
    }
}

#endif