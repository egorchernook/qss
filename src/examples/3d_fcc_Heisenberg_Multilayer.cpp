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

int main()
{
    using spin_t = qss::models::heisenberg::spin;
    using lattice_t = qss::lattices::three_d::fcc<spin_t>;
    using sizes_t = qss::lattices::three_d::sizes_t;
    using qss::nanostructures::film;
    using qss::nanostructures::multilayer;

    constexpr static sizes_t sizes{16, 16, 3};
    const film<lattice_t> fst_film{lattice_t{spin_t{1.0, 0.0, 0.0}, sizes}, 1.0};
    const film<lattice_t> snd_film{lattice_t{spin_t{-1.0, 0.0, 0.0}, sizes}, 1.0};

    multilayer structure{fst_film};
    structure.add(snd_film, -0.1);

    constexpr static std::uint32_t mcs_amount = 5'000;
    
    structure.T = 0.5;
    for (std::size_t mcs = 0; mcs <= mcs_amount; ++mcs)
    {
        structure.evolve();
    }

    return 0;
}