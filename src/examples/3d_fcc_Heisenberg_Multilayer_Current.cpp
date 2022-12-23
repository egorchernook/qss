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
#include "../algorithms/spin_transport.hpp"
#include "../models/electron_dencity.hpp"

int main()
{
    using spin_t = qss::models::heisenberg::spin;
    using lattice_t = qss::lattices::three_d::fcc<spin_t>;
    using sizes_t = qss::lattices::three_d::sizes_t;
    using qss::nanostructures::film;
    using qss::nanostructures::multilayer;
    using qss::nanostructures::multilayer_system;
    using ed_t = qss::models::electron_dencity;
    using electron_dencity_t = qss::lattices::three_d::fcc<ed_t>;

    constexpr static sizes_t sizes{64, 64, 3};
    multilayer_system<multilayer<lattice_t>> system{
        multilayer{{film<lattice_t>{lattice_t{spin_t{1.0, 0.0, 0.0}, sizes}, 1.0},
                    film<lattice_t>{lattice_t{spin_t{-1.0, 0.0, 0.0}, sizes}, 1.0}},
                   {-0.1}}};
    multilayer n_up{{film<electron_dencity_t>{electron_dencity_t{ed_t{0.0}, sizes}, 1.0},
                     film<electron_dencity_t>{electron_dencity_t{ed_t{0.0}, sizes}, 1.0}},
                    {-0.1}};
    multilayer n_down{{film<electron_dencity_t>{electron_dencity_t{ed_t{0.0}, sizes}, 1.0},
                       film<electron_dencity_t>{electron_dencity_t{ed_t{0.0}, sizes}, 1.0}},
                      {-0.1}};

    auto sys = qss::algorithms::spin_transport::prepare_proxy_structure<'x'>(system, n_up, n_down);

    constexpr static std::uint32_t mcs_amount = 10'000;
    constexpr static double Delta = 0.665;
    system.T = 0.5;
    sys.T = system.T;
    std::ofstream out_magn{"m.txt"};
    std::ofstream out_j{"j.txt"};
    for (std::size_t mcs = 0; mcs <= mcs_amount; ++mcs)
    {
        if (mcs % 10 == 0)
        {
            std::cout << mcs << "\n";
        }
        double j_up_all = 0.0;
        double j_down_all = 0.0;
        if (mcs >= 2000)
        {
            const auto temp_magn1 = abs(system.magns[0]);
            const auto temp_magn2 = -abs(system.magns[1]);
            if (mcs == 2000)
            {
                n_up[0].fill(ed_t{0.5 * (1.0 + temp_magn1)});
                n_up[1].fill(ed_t{0.5 * (1.0 + temp_magn1)});
                n_down[0].fill(ed_t{0.5 * (1.0 - temp_magn2)});
                n_down[1].fill(ed_t{0.5 * (1.0 - temp_magn2)});
                system.T = 1.0;
                sys.T = system.T;
            }
            n_up[0].fill_plane(0, ed_t{0.5 * (1.0 + temp_magn1)});
            n_down[0].fill_plane(0, ed_t{0.5 * (1.0 - temp_magn2)});
            const auto [j_up, j_down] = qss::algorithms::spin_transport::perform(sys);
            j_up_all += j_up / (sizes.x * sizes.y);
            j_down_all += j_down / (sizes.x * sizes.y);
            out_j << mcs << "\t"
                  << j_up_all << "\t"
                  << j_down_all << std::endl;
        }

        system.evolve([](const typename spin_t::magn_t &sum,
                         const spin_t &spin_old,
                         const spin_t &spin_new) -> double
                      {
                                auto diff = spin_old - spin_new;
                                diff.y *= 0.8;
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
    out_j.flush();
    out_magn.close();
    out_j.close();

    return 0;
}