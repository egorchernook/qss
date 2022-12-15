#ifndef SPIN_TRANSPORT_HPP_INCLUDED
#define SPIN_TRANSPORT_HPP_INCLUDED

#include <concepts>
#include <functional>
#include <stdexcept>

#include "../models/heisenberg.hpp"
#include "../systems/multilayer_system.hpp"
#include "../systems/multilayer.hpp"
#include "../systems/film.hpp"
#include "../random/random.hpp"
#include "../random/mersenne.hpp"
#include "../models/electron_dencity.hpp"

namespace qss::algorithms::spin_transport
{
    struct proxy_spin
    {
        using magn_t = double;
        std::reference_wrapper<double> value;
        std::reference_wrapper<typename qss::models::electron_dencity> up;
        std::reference_wrapper<typename qss::models::electron_dencity> down;

        constexpr proxy_spin(double &value_,
                             typename qss::models::electron_dencity &up_,
                             typename qss::models::electron_dencity &down_) noexcept
            : value{std::ref(value_)}, up{std::ref(up_)}, down{std::ref(down_)} {}

        operator magn_t() const
        {
            return value * (typename qss::models::electron_dencity::magn_t(up.get()) - typename qss::models::electron_dencity::magn_t(down.get()));
        }
    };

    inline typename proxy_spin::magn_t operator+(const proxy_spin &lhs, const proxy_spin &rhs) noexcept
    {
        typename proxy_spin::magn_t result{};
        result = typename proxy_spin::magn_t(lhs) + typename proxy_spin::magn_t(rhs);
        return result;
    }
    inline typename proxy_spin::magn_t operator-(const proxy_spin &lhs, const proxy_spin &rhs) noexcept
    {
        typename proxy_spin::magn_t result{};
        result = typename proxy_spin::magn_t(lhs) - typename proxy_spin::magn_t(rhs);
        return result;
    }
    inline double scalar_multiply(const proxy_spin &lhs, const proxy_spin &rhs) noexcept
    {
        return typename proxy_spin::magn_t(lhs) * typename proxy_spin::magn_t(rhs);
    }
    inline double scalar_multiply(const typename proxy_spin::magn_t &lhs, const typename proxy_spin::magn_t &rhs) noexcept
    {
        return lhs * rhs;
    }
    inline typename proxy_spin::magn_t operator*(const double &lhs, const proxy_spin &rhs) noexcept
    {
        return lhs * typename proxy_spin::magn_t(rhs);
    }
    inline typename proxy_spin::magn_t operator*(const proxy_spin &lhs, const double &rhs) noexcept
    {
        return (rhs * lhs);
    }
    inline typename proxy_spin::magn_t operator/(const proxy_spin &lhs, const double &rhs) noexcept
    {
        return typename proxy_spin::magn_t(lhs) / rhs;
    }

    template <template <typename> class lattice_t, typename spin_t>
        requires qss::nanostructures::ThreeD_Lattice<lattice_t<spin_t>>
    using nanostructure_type =
        qss::nanostructures::multilayer_system<
            qss::nanostructures::multilayer<
                lattice_t<spin_t>>>;

    template <char spin_component_name,
              typename old_spin_t,
              template <typename = old_spin_t> class film_t>
    [[nodiscard]] nanostructure_type<film_t, proxy_spin>
    prepare_proxy_structure(nanostructure_type<film_t, old_spin_t> &system,
                            typename qss::nanostructures::multilayer<
                                film_t<
                                    qss::models::electron_dencity>> &n_up,
                            typename qss::nanostructures::multilayer<
                                film_t<
                                    qss::models::electron_dencity>> &n_down)
    {
        auto result = copy_structure<proxy_spin>(system);
        {
            auto result_iter = result.nanostructure.begin();
            auto system__iter = system.nanostructure.begin();
            auto system_up_iter = n_up.begin();
            auto system_down_iter = n_down.begin();
            for (;
                 result_iter != result.nanostructure.end();
                 ++result_iter, ++system__iter, ++system_up_iter, ++system_down_iter)
            {
                auto iter = result_iter->begin();
                auto sys_iter = system__iter->begin();
                auto up_iter = system_up_iter->begin();
                auto down_iter = system_down_iter->begin();
                for (;
                     iter != result_iter->end();
                     ++iter, ++sys_iter, ++up_iter, ++down_iter)
                {
                    if (spin_component_name == 'x')
                    {
                        *iter = proxy_spin{sys_iter->x, *up_iter, *down_iter};
                    }
                    else if (spin_component_name == 'y')
                    {
                        *iter = proxy_spin{sys_iter->y, *up_iter, *down_iter};
                    }
                    else if (spin_component_name == 'z')
                    {
                        *iter = proxy_spin{sys_iter->z, *up_iter, *down_iter};
                    }
                    else
                    {
                        throw std::logic_error("spin_component_name should be 'x' or 'y' or 'z' but was : " + std::to_string(spin_component_name));
                    }
                }
            }
        }
        return result;
    }

    struct result_t
    {
        double up;
        double down;
    };

    template <template <typename> class film_t,
              typename random_t = qss::random::mersenne::random_t<>>
    result_t perform(const nanostructure_type<film_t, proxy_spin> &system) noexcept
    {
        result_t result{};
        const auto layers = system.nanostructure;
        const auto amount = layers[layers.template get_random_coord<random_t>().idx].get_amount_of_nodes();
        for (auto i = 0u; i < amount; ++i)
        {
            const auto coord = layers.template get_random_coord<random_t>();
            const auto E1 = scalar_multiply(layers.get_sum_of_closest_neighbours(coord),
                                            layers.get(coord));
            double E2{};
            if (coord.film_coord.z == layers[coord.idx].get_last_z(coord.film_coord))
            {
                E2 = 0.0;
            }
            else
            {
                auto next_coord = coord.film_coord;
                next_coord.z += 1;
                E2 = scalar_multiply(layers.get_sum_of_closest_neighbours(coord),
                                     layers.get({coord.idx, next_coord}));
            }
            const auto delta_E = E2 - E1;
            if (delta_E < 0.0)
            {
                result.up += layers.get(coord).up.get();
                result.down += layers.get(coord).down.get();
            }
        }
        return result;
    }
}

#endif