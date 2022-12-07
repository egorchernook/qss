#ifndef SPIN_TRANSPORT_HPP_INCLUDED
#define SPIN_TRANSPORT_HPP_INCLUDED

#include "../models/electron_dencity.hpp"
#include "../models/heisenberg.hpp"

namespace qss::algorithms::spin_transport
{
    template <typename multilayer_t, typename electron_dencity_multilayer_t>
    void make(const multilayer_t &system,
              const electron_dencity_multilayer_t &n_up,
              const electron_dencity_multilayer_t &n_down) noexcept
    {
        const auto coord = system.template get_random_coord<>();
    }
}

#endif