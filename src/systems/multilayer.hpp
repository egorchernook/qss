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

#include "../lattices/3d/3d.hpp"
#include "../lattices/3d/fcc.hpp"
#include "../lattices/borders_conditions.hpp"
#include "../utility/functions.hpp"
#include "../random/random.hpp"
#include "../random/mersenne.hpp"

namespace qss::nanostructures
{
    template <typename lattice_t>
    concept ThreeD_Lattice = std::same_as<typename lattice_t::sizes_t, qss::lattices::three_d::sizes_t>;

    template <ThreeD_Lattice lattice_t>
    struct film final : public lattice_t
    {
        using xy_border_condition = typename qss::border_conditions::periodic<typename lattice_t::coords_t::size_type, typename lattice_t::sizes_t::size_type>;
        using z_border_condition = typename qss::border_conditions::sharp<typename lattice_t::coords_t::size_type, typename lattice_t::sizes_t::size_type>;
        double J;

        film(const lattice_t &lattice, double J_ = 1.0) : lattice_t{lattice}, J{J_} {};
        film(lattice_t &&lattice, double J_ = 1.0) noexcept : lattice_t{std::move(lattice)}, J{J_} {};
    };

    template <ThreeD_Lattice lattice_t>
    requires std::is_same_v<typename lattice_t::coords_t, qss::lattices::three_d::fcc_coords_t>
        std::optional<qss::lattices::three_d::fcc_coords_t>
        get_closest_neigbour_from_upper_film(const film<lattice_t> &other,
                                             const qss::lattices::three_d::fcc_coords_t &coord)
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
        std::optional<qss::lattices::three_d::fcc_coords_t>
        get_closest_neigbour_from_lower_film(const film<lattice_t> &other,
                                             const qss::lattices::three_d::fcc_coords_t &coord)
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
    class multilayer final
    {
        using container_t = typename std::vector<film<lattice_t>>;
        container_t films;
        std::vector<double> J_interlayers; // обменные интегралы взаимодействий плёнок.

        void check_sizes(const typename film<lattice_t>::sizes_t &value) const
        {
            if (value.x != films[0].sizes.x)
            {
                throw std::out_of_range("x size must be the same for all lattices : " + std::to_string(value.x) + " != " + std::to_string(films[0].sizes.x));
            }
            if (value.y != films[0].sizes.y)
            {
                throw std::out_of_range("y size must be the same for all lattices : " + std::to_string(value.y) + " != " + std::to_string(films[0].sizes.y));
            }
        }

    public:
        using film_t = film<lattice_t>;
        struct coords_t
        {
            typename film_t::coords_t film_coord;
            typename container_t::size_type idx;
        };

        constexpr multilayer(const film_t &film)
            : films{film},
              magns{calculate_magn(film)},
              energies{0.0} {};
        constexpr multilayer(film_t &&film)
            : films{std::move(film)},
              magns{calculate_magn(film)},
              energies{0.0} {};

        void add(const film_t &film, double J_interlayer)
        {
            check_sizes(film.sizes);
            films.push_back(film);
            magns.push_back(calculate_magn(film));
            energies.push_back(0.0);
            J_interlayers.push_back(J_interlayer);
        }
        void add(film_t &&film, double J_interlayer)
        {
            check_sizes(film.sizes);
            films.push_back(std::move(film));
            magns.push_back(calculate_magn(film));
            energies.push_back(0.0);
            J_interlayers.push_back(J_interlayer);
        }
        film_t pop()
        {
            film_t result = films.back();
            films.pop_back();
            magns.pop_back();
            energies.pop_back();
            J_interlayers.pop_back();
            return result;
        }

    private:
        typename film_t::value_t::magn_t
        get_sum_of_closest_neighbours(typename container_t::size_type idx,
                                      const typename film_t::coords_t &central) const
        {
            const auto film = films[idx];
            bool has_upper = idx != films.size() - 1u;
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
                        const auto neig = get_closest_neigbour_from_upper_film(films[idx + 1u], central);
                        upper_neig_val += neig ? films[idx + 1u].get(neig.value()) : magn_t{};
                    }
                    if (has_lower)
                    {
                        const auto neig = get_closest_neigbour_from_lower_film(films[idx - 1u], central);
                        lower_neig_val += neig ? films[idx - 1u].get(neig.value()) : magn_t{};
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

        std::vector<typename film_t::value_t::magn_t> magns;
        std::vector<double> energies;

    public:
        std::vector<typename film_t::value_t::magn_t> get_magns() const noexcept
        {
            return magns;
        }
        std::vector<double> get_energies() const noexcept
        {
            return energies;
        }
        double T = 0.0;

        using delta_h_t = auto(*)(const typename film_t::value_t::magn_t &,
                                  const typename film_t::value_t &,
                                  const typename film_t::value_t &) -> double;
        /*
         * использует алгоритм Метрополиса
         * необходимо установить температуру, перед использованием
         **/
        void evolve(delta_h_t delta_h = [](const typename film_t::value_t::magn_t &sum,
                                           const typename film_t::value_t &spin_old,
                                           const typename film_t::value_t &spin_new) -> double
                    {
                        return scalar_multiply(sum, spin_old - spin_new);
                    })
        {
            for (auto idx = 0u; idx < films.size(); ++idx)
            {
                auto delta_energy_f = [&delta_h, &idx, this](const lattice_t &lattice_,
                                                             const typename lattice_t::coords_t &central,
                                                             const typename lattice_t::value_t &new_spin)
                    -> double
                {
                    const auto sum = get_sum_of_closest_neighbours(idx, central);
                    return delta_h(sum, lattice_.get(central), new_spin);
                };

                auto [M, E] = qss::algorithms::metropolis::make_step(films[idx], delta_energy_f, T);
                magns[idx] += M / static_cast<double>(films[idx].get_amount_of_nodes());
                energies[idx] += -0.5 * E / static_cast<double>(films[idx].get_amount_of_nodes());
            }
        }
    };
}

#endif