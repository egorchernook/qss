#ifndef METROPOLIS_HPP_INCLUDED
#define METROPOLIS_HPP_INCLUDED

#include <type_traits>

#include "../random/random.hpp"
#include "../random/mersenne.hpp"

namespace qss::algorithms::metropolis
{
    /*
     * Свободная процедура для прохождения одного шага Монте-Карло.
    **/
    template <typename lattice_t,
              typename delta_energy_f_t,
              Random random_t = qss::random::mersenne::random_t<>> // TODO: ограничить typename и auto
    void make_step(lattice_t &lattice,
                   delta_energy_f_t delta_energy_f,
                   double temperature)
    {
        static random_t rand{qss::random::get_seed()};
        for (auto _ = 0llu; _ < lattice.get_amount_of_nodes(); ++_)
        {
            const auto old_spin_coords = lattice.template choose_random_node<random_t>();
            const auto spin_new = lattice_t::value_t::template generate<random_t>();

            const double dE = delta_energy_f(lattice, old_spin_coords, spin_new); // E_old - E_new

            if (dE < 0.0 || rand() < std::exp(-dE / temperature))
            {
                lattice.set(spin_new, old_spin_coords);
            }
        }
    }
}

#endif