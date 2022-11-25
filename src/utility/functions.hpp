#ifndef FUNCTIONS_HPP_INCLUDED
#define FUNCTIONS_HPP_INCLUDED

#include <type_traits>
#include <cassert>
#include <limits>

namespace qss
{
    template <typename lattice_t,
              typename borders_conditions_t> // TODO: ограничить typenames
    typename lattice_t::value_t::magn_t get_sum_of_closest_neighbours(const lattice_t &lattice,
                                         const typename lattice_t::coords_t &central,
                                         borders_conditions_t borders_conditions)
    {
        using magn_t = typename lattice_t::value_t::magn_t;
        auto neigs = get_closest_neighbours(central);
        const auto sizes = lattice.sizes;

        magn_t sum{};
        for (auto it = neigs.begin(); it != neigs.end(); ++it)
        {
            const auto coord = borders_conditions(*it, sizes);
            if (coord)
            {
                sum += lattice.get(coord.value());
            }
        }
        return sum;
    }
}

#endif