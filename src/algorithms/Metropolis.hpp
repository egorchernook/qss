#ifndef METROPOLIS_HPP_INCLUDED
#define METROPOLIS_HPP_INCLUDED

#include "../random/mersenne.hpp"
#include "../random/random.hpp"

#include <type_traits>

namespace qss {
inline namespace algorithms {
namespace metropolis {
/*
 * Свободная процедура для прохождения одного шага Монте-Карло.
 * возвращает std::pair{ изменение намагниченности (ненормированная), изменение энергии }
 **/
template<
    typename lattice_t,
    typename delta_energy_f_t,
    Random random_t = qss::random::mersenne::random_t<>> // TODO: ограничить typename и auto
std::pair<typename lattice_t::value_t::magn_t, double>
make_step(lattice_t& lattice, delta_energy_f_t delta_energy_f, double temperature)
{
    static random_t rand{qss::random::get_seed()};
    double delta_energy = 0.0;
    typename lattice_t::value_t::magn_t delta_magn{};
    for (auto _ = 0llu; _ < lattice.get_amount_of_nodes(); ++_) {
        const auto old_spin_coords = lattice.template choose_random_node<random_t>();
        const auto spin_new = lattice_t::value_t::template generate<random_t>();

        const double dE = delta_energy_f(lattice, old_spin_coords, spin_new); // E_old - E_new
        const auto old_spin = lattice.get(old_spin_coords);
        if (dE < 0.0 || rand() < std::exp(-dE / temperature)) {
            lattice.set(spin_new, old_spin_coords);
            delta_energy += dE;
            delta_magn += spin_new - old_spin;
        }
    }
    return std::pair{delta_magn, delta_energy};
}
} // namespace metropolis
} // namespace algorithms
} // namespace qss

#endif