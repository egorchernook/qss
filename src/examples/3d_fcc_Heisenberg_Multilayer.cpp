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
#include "../systems/multilayer.hpp"
#include "../systems/multilayer_system.hpp"

int main()
{
    using spin_t = qss::heisenberg::spin;
    using lattice_t = qss::lattices::three_d::fcc<spin_t>;
    using sizes_t = qss::lattices::three_d::sizes_t;
    using qss::film;
    using qss::multilayer;
    using qss::multilayer_system;

    constexpr static sizes_t sizes{64, 64, 3};
    multilayer_system<multilayer<lattice_t>> system{
        multilayer{{film<lattice_t>{lattice_t{spin_t{1.0, 0.0, 0.0}, sizes}, 1.0},
                    film<lattice_t>{lattice_t{spin_t{-1.0, 0.0, 0.0}, sizes}, 1.0}},
                   {-0.1}}};

    constexpr static std::uint32_t mcs_amount = 2'000;
    constexpr static double Delta = 0.665;
    system.T = 0.5;
    std::ofstream out_magn{"m.txt"};
    for (std::size_t mcs = 0; mcs <= mcs_amount; ++mcs)
    {
        if (mcs % 10 == 0)
        {
            std::cout << mcs << "\n";
        }
        system.evolve([](const typename spin_t::magn_t &sum,
                            const spin_t &spin_old,
                            const spin_t &spin_new) -> double
                         {
                                auto diff = spin_old - spin_new;
                                diff.z *= (1.0 - Delta);
                                return scalar_multiply(sum, diff); });
        const auto magn1 = system.magns[0];
        const auto magn2 = system.magns[1];

        out_magn << mcs << "\t"
                 << abs(magn1) - abs(magn2) << "\t"
                 << abs(magn1) << "\t"
                 << magn1 << "\t"
                 << abs(magn2) << "\t"
                 << magn2 << std::endl;
    }
    out_magn.flush();
    out_magn.close();
    return 0;
}