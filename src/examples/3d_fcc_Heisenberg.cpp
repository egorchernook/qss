#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "../algorithms/Metropolis.hpp"
#include "../models/heisenberg.hpp"
#include "../lattices/3d/fcc.hpp"
#include "../lattices/3d/3d.hpp"
#include "../lattices/borders_conditions.hpp"
#include "../utility/quantities.hpp"
#include "../utility/functions.hpp"

std::vector<double> get_temperatures(const double T_begin = 1.0,
                                     const double T_end = 5.0,
                                     const double delta_T = 0.25)
{
    std::vector<double> result{};
    const auto T_amount = static_cast<int>((T_end - T_begin) / delta_T) + 1;
    for (int i = 0; i < T_amount; ++i)
    {
        result.push_back(T_begin + i * delta_T);
    }
    return result;
}

int main()
{
    using spin_t = qss::models::heisenberg::spin;
    using lattice_t = qss::lattices::three_d::fcc<spin_t>;
    using sizes_t = qss::lattices::three_d::sizes_t;
    using periodic = qss::borders_conditions::periodic<typename lattice_t::coords_t::size_type,
                                                       typename sizes_t::size_type>;
    using sharp = qss::borders_conditions::sharp<typename lattice_t::coords_t::size_type,
                                                 typename sizes_t::size_type>;

    constexpr static sizes_t sizes{16, 16, 3};
    lattice_t lattice{spin_t{1.0, 0.0, 0.0}, sizes};

    constexpr static std::uint32_t mcs_amount = 5'000;
    const std::vector temperatures = get_temperatures();

    auto delta_energy_f =
        [](const lattice_t &lattice_,
           const lattice_t::coords_t &central,
           const spin_t &new_spin)
        -> double
    {
        const auto sum =
            qss::get_sum_of_closest_neighbours(
                lattice_,
                central,
                qss::borders_conditions::use_border_conditions<periodic, periodic, sharp>);

        return scalar_multiply(sum, (lattice_.get(central) - new_spin));
    };

    std::ofstream output{"m.txt"};
    for (auto T : temperatures)
    {
        std::cout << "T = " << T << std::endl;
        for (std::size_t mcs = 0; mcs <= mcs_amount; ++mcs)
        {
            if (mcs % 100 == 0)
            {
                const auto magn = qss::calculate_magn(lattice);
                const auto absl = abs(magn);
                output << mcs << "\t"
                       << T << "\t"
                       << absl << "\t"
                       << magn << "\t"
                       << "\n";
            }
            qss::algorithms::metropolis::make_step(lattice, delta_energy_f, T);
        }
    }
    output.flush();
    output.close();

    return 0;
}