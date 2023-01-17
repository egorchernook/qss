#ifndef QUANTITIES_HPP_INCLUDED
#define QUANTITIES_HPP_INCLUDED

#include <type_traits>
#include <numeric>
#include <algorithm>
#include <execution>

#include "functions.hpp"

namespace qss
{
    template <typename lattice_t> // TODO: добавить ограничений
    typename lattice_t::value_t::magn_t calculate_magn(const lattice_t &lattice)
    {
        using magn_t = typename lattice_t::value_t::magn_t;
        magn_t magn{};
        
        magn = std::reduce(std::execution::par_unseq, lattice.begin(), lattice.end(), magn);
        return magn / static_cast<double>(lattice.get_amount_of_nodes());
    }

    /*
    template <typename lattice_t, typename borders_conditions_t>
    double calculate_energy(const lattice_t &lattice, borders_conditions_t borders_conditions)
    {
        using magn_t = typename lattice_t::value_t::magn_t;
        const auto coords_repr = lattice.as_coords();
        double energy = 0.0;
        for (auto elem : coords_repr)
        {
            const auto sum = qss::get_sum_of_closest_neighbours(
                lattice,
                elem,
                borders_conditions);
            energy += sum * static_cast<magn_t>(lattice.get(elem));
        }

        return 0.5 * energy / static_cast<double>(lattice.get_amount_of_nodes());
    }
    */
}

#endif