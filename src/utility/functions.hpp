#ifndef FUNCTIONS_HPP_INCLUDED
#define FUNCTIONS_HPP_INCLUDED

#include <type_traits>
#include <cassert>

namespace qss
{
    template <
        typename magn_t,
        typename lattice_t,
        typename borders_conditions_t,
        std::enable_if_t<
            std::is_convertible_v<typename lattice_t::value_t, magn_t>,
            bool> = true> // TODO: ограничить typenames
    magn_t get_sum_of_closest_neighbours(const lattice_t &lattice,
                                         const typename lattice_t::coords_t &central,
                                         borders_conditions_t borders_conditions)
    {
        auto neigs = lattice.get_closest_neighbours(central);
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