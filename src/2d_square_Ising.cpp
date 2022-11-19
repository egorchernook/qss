#include <fstream>
#include <iostream>
#include <numeric>

#include "algorithms/Metropolis.hpp"
#include "models/ising.hpp"
#include "lattices/2d/square.hpp"
#include "lattices/2d/2d.hpp"
#include "lattices/borders_conditions.hpp"
#include "utility/quantities.hpp"
#include "utility/functions.hpp"

std::vector<double> get_temperatures(const double T_begin = 1.5,
                                     const double T_end = 4.0,
                                     const double delta_T = 0.25)
{
    std::vector<double> result{};
    const auto T_amount = static_cast<int>((T_end - T_begin) / delta_T);
    for (int i = 0; i < T_amount; ++i)
    {
        result.push_back(T_begin + i * delta_T);
    }
    return result;
}

int main()
{
    using spin_t = qss::models::ising::spin;
    using magn_t = qss::models::ising::magn;
    using lattice_t = qss::lattices::two_d::square<spin_t>;
    using sizes_t = qss::lattices::two_d::sizes_t<>;

    constexpr static sizes_t sizes{64, 64};
    lattice_t lattice{spin_t{1}, sizes};

    constexpr static std::uint32_t mcs_amount = 2'000;
    const std::vector temperatures = get_temperatures();

    auto delta_energy_f =
        [](const lattice_t &lattice_,
           const lattice_t::coords_t &central,
           const spin_t &new_spin)
        -> double
    {
        const auto sum =
            qss::get_sum_of_closest_neighbours<magn_t>(
                lattice_,
                central,
                qss::lattices::two_d::use_periodic_conditions<>);

        return sum * (lattice_.get(central) - new_spin);
    };

    std::ofstream output{"m.txt"};
    for (auto T : temperatures)
    {
        std::cout << "T = " << T << std::endl;
        for (std::size_t mcs = 0; mcs <= mcs_amount; ++mcs)
        {
            if (mcs % 100 == 0)
            {
                output << mcs << "\t"
                       << T << "\t"
                       << std::abs(qss::calculate_magn<lattice_t, magn_t>(lattice))
                       << "\n";
            }
            qss::algorithms::metropolis::make_step(lattice, delta_energy_f, T);
        }
    }
    output.flush();
    output.close();

    return 0;
}