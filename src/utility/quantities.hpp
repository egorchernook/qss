#ifndef QUANTITIES_HPP_INCLUDED
#define QUANTITIES_HPP_INCLUDED

#include <type_traits>

#include "functions.hpp"

namespace qss
{
    template <typename lattice_t,
              typename magn_t,
              std::enable_if_t<
                  std::is_convertible_v<typename lattice_t::value_t, magn_t>,
                  bool> = true> // TODO: добавить ограничений
    magn_t calculate_magn(const lattice_t &lattice)
    {
        magn_t magn{};
        for (const auto &elem : lattice)
        {
            magn += static_cast<magn_t>(elem);
        }
        return magn / static_cast<double>(lattice.get_amount_of_nodes());
    }

    template <typename lattice_t, typename magn_t, typename borders_conditions_t>
    double calculate_energy(const lattice_t &lattice, borders_conditions_t borders_conditions)
    {
        const auto coords_repr = lattice.as_coords();
        double energy = 0.0;
        for (auto elem : coords_repr)
        {
            const auto sum = qss::get_sum_of_closest_neighbours<magn_t>(
                lattice,
                elem,
                borders_conditions);
            energy += sum * static_cast<magn_t>(lattice.get(elem));
        }

        return 0.5 * energy / static_cast<double>(lattice.get_amount_of_nodes());
    }
}

#endif