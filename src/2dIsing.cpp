#include "algorithms/Metropolis.hpp"
#include "models/ising.hpp"
#include "lattices/2d/square.hpp"
#include "lattices/2d/2d.hpp"

#include <fstream>
#include <iostream>
#include <numeric>

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

    auto calculate_magn = [](const lattice_t &lattice) -> magn_t
    {
        magn_t magn{};
        for (const auto &elem : lattice)
        {
            magn += static_cast<magn_t>(elem);
        }
        return std::abs(magn / lattice.get_amount_of_nodes());
    };
    auto delta_energy_f = [](const lattice_t &lattice,
                             const lattice_t::coords_t &central,
                             const spin_t &new_spin)
        -> double
    {
        auto neig = lattice.get_closest_neighbours(central);
        // периодические граничные условия
        for (auto &coord : neig)
        {
            if (coord.x == lattice.sizes.x)
            {
                coord.x = 0;
            }
            if (coord.y == lattice.sizes.y)
            {
                coord.y = 0;
            }
            if (coord.x == -1)
            {
                coord.x = lattice.sizes.x - 1;
            }
            if (coord.y == -1)
            {
                coord.y = lattice.sizes.y - 1;
            }
        }
        // --------------------------------
        magn_t sum{};
        for (auto &elem : neig)
        {
            sum += lattice.get(elem);
        }
        return sum * (lattice.get(central) - new_spin);
    };

    std::ofstream output{"m.txt"};
    for (auto T : temperatures)
    {
        std::cout << "T = " << T << std::endl;
        for (std::size_t mcs = 0; mcs <= mcs_amount; ++mcs)
        {
            if (mcs % 100 == 0)
            {
                output << mcs << "\t" << T << "\t" << calculate_magn(lattice) << "\n";
            }
            qss::algorithms::metropolis::make_step(lattice, delta_energy_f, T);
        }
    }
    //output << T << "\t" << calculate_magn(lattice) << "\n";
    output.flush();
    output.close();

    return 0;
}