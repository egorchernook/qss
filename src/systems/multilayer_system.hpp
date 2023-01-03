#ifndef MULTILAYER_SYSTEM_HPP_INCLUDED
#define MULTILAYER_SYSTEM_HPP_INCLUDED

#include <vector>

#include "multilayer.hpp"

namespace qss::inline nanostructures
{
    template <typename multilayer_t>
    struct multilayer_system
    {
        multilayer_t nanostructure;
        std::vector<typename multilayer_t::film_t::value_t::magn_t> magns;
        std::vector<double> energies;
        double T;

        [[nodiscard]] constexpr multilayer_system(multilayer_t &&structure) noexcept
            : nanostructure{std::move(structure)}
        {
            magns.reserve(nanostructure.size());
            energies.reserve(nanostructure.size());
            for (auto &film : nanostructure)
            {
                magns.push_back(calculate_magn(film));
                energies.push_back(0.0);
            }
        }
        [[nodiscard]] constexpr multilayer_system(const multilayer_t &structure) noexcept
            : nanostructure{structure}
        {
            magns.reserve(nanostructure.size());
            energies.reserve(nanostructure.size());
            for (auto &film : nanostructure)
            {
                magns.push_back(calculate_magn(film));
                energies.push_back(0.0);
            }
        }

        /*
        using delta_h_t = auto(*)(const typename multilayer_t::film_t::value_t::magn_t &,
                                  const typename multilayer_t::film_t::value_t &,
                                  const typename multilayer_t::film_t::value_t &) -> double;
        */
        /*
         * использует алгоритм Метрополиса
         * необходимо установить температуру, перед использованием
         **/
        template<typename delta_h_t>
        void evolve(delta_h_t delta_h = [](const typename multilayer_t::film_t::value_t::magn_t &sum,
                                           const typename multilayer_t::film_t::value_t &spin_old,
                                           const typename multilayer_t::film_t::value_t &spin_new) -> double
                    {
                        return scalar_multiply(sum, spin_old - spin_new);
                    }) noexcept
        {
            for (std::uint8_t idx = 0; idx < nanostructure.size(); ++idx)
            {
                auto delta_energy_f = [&delta_h,
                                       &idx,
                                       this](const typename multilayer_t::film_t &lattice_,
                                             const typename multilayer_t::film_t::coords_t &central,
                                             const typename multilayer_t::film_t::value_t &new_spin)
                    -> double
                {
                    const auto sum = nanostructure.get_sum_of_closest_neighbours({idx, central});
                    return delta_h(sum, lattice_.get(central), new_spin);
                };

                auto [M, E] = qss::algorithms::metropolis::make_step(nanostructure[idx], delta_energy_f, T);
                magns[idx] += M / static_cast<double>(nanostructure[idx].get_amount_of_nodes());
                energies[idx] += -0.5 * E / static_cast<double>(nanostructure[idx].get_amount_of_nodes());
            }
        }
    };

    template <typename spin_t,
              typename old_spin_t,
              template <typename = old_spin_t> class lattice_t>
    [[nodiscard]] constexpr multilayer_system<multilayer<lattice_t<spin_t>>>
    copy_structure(const multilayer_system<multilayer<lattice_t<old_spin_t>>> &original) noexcept
    {
        if constexpr (std::is_same_v<old_spin_t, spin_t>)
        {
            return original;
        }
        return {copy_structure<spin_t>(original.nanostructure)};
    }
}

#endif