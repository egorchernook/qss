#ifndef SYSTEM_HPP_INCLUDED
#define SYSTEM_HPP_INCLUDED

#include <tuple>
#include <array>

namespace qss
{
    template <typename delta_energy_f, typename ... Tp>
    struct system {
        std::tuple<Tp...> data;
        delta_energy_f delta_energy_function;
        std::array<double, 2 * sizeof...(Tp) - 1> exchange_integrals;
    };
}

#endif