#ifndef SPIN_TRANSPORT_HPP_INCLUDED
#define SPIN_TRANSPORT_HPP_INCLUDED

// #include <concepts>
#include <numeric>
#include <stdexcept>

#include "../models/electron_dencity.hpp"
#include "../models/heisenberg.hpp"
#include "../random/mersenne.hpp"
#include "../random/random.hpp"
#include "../systems/film.hpp"
#include "../systems/multilayer.hpp"
#include "../systems/multilayer_system.hpp"

namespace qss::inline algorithms::spin_transport
{
    // используется ТОЛЬКО для создания proxy структуры, для расчёта тока
    class proxy_spin
    {
        double *value = nullptr;
        typename qss::models::electron_dencity *up = nullptr;
        typename qss::models::electron_dencity *down = nullptr;

    public:
        proxy_spin() noexcept = default;
        proxy_spin(const proxy_spin &) noexcept = default;
        proxy_spin(proxy_spin &&) noexcept = default;
        ~proxy_spin() noexcept = default;
        proxy_spin &operator=(const proxy_spin &other) noexcept
        {
            value = other.value;
            up = other.up;
            down = other.down;
            return *this;
        }

        using magn_t = double;

        double get_val() const noexcept
        {
            return *value;
        }
        typename qss::models::electron_dencity get_up() const noexcept
        {
            return *up;
        }
        typename qss::models::electron_dencity get_down() const noexcept
        {
            return *down;
        }
        void set_up(typename qss::models::electron_dencity value_) noexcept
        {
            *up = value_;
        }
        void set_down(typename qss::models::electron_dencity value_) noexcept
        {
            *down = value_;
        }

        constexpr proxy_spin(double &value_, typename qss::models::electron_dencity &up_,
                             typename qss::models::electron_dencity &down_) noexcept
            : value{&value_}, up{&up_}, down{&down_}
        {
        }

        operator magn_t() const
        {
            if (value == nullptr || up == nullptr || down == nullptr)
            {
                return 0.0;
            }
            if (*value > 0.0)
            {
                return typename qss::models::electron_dencity::magn_t(*up) -
                       typename qss::models::electron_dencity::magn_t(*down);
            }
            else
            {
                return -typename qss::models::electron_dencity::magn_t(*up) +
                       typename qss::models::electron_dencity::magn_t(*down);
            }
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
        // requires qss::nanostructures::ThreeD_Lattice<lattice_t<spin_t>>
    using nanostructure_type = qss::nanostructures::multilayer_system<qss::nanostructures::multilayer<lattice_t<spin_t>>>;

    /*  proxy структура ТОЛЬКО для использования в функции spin_transport
     *  содержит указатели на нужные значения спина, и электронных плотностей
     *  spin_component_name принимает значения x, y или z,
     *  чтобы указать какую составляющую спина использовать далее
     **/
    template <typename old_spin_t, template <typename = old_spin_t> class lattice_t>
    [[nodiscard]] nanostructure_type<lattice_t, proxy_spin> prepare_proxy_structure(
        nanostructure_type<lattice_t, old_spin_t> &system,
        typename qss::nanostructures::multilayer<lattice_t<qss::models::electron_dencity>> &n_up,
        typename qss::nanostructures::multilayer<lattice_t<qss::models::electron_dencity>> &n_down,
        char spin_component_name = 'x')
    {
        auto result = copy_structure<proxy_spin>(system);
        {
            auto result_iter = result.nanostructure.begin();
            auto system__iter = system.nanostructure.begin();
            auto system_up_iter = n_up.begin();
            auto system_down_iter = n_down.begin();
            for (; result_iter != result.nanostructure.end();
                 ++result_iter, ++system__iter, ++system_up_iter, ++system_down_iter)
            {
                auto iter = result_iter->begin();
                auto sys_iter = system__iter->begin();
                auto up_iter = system_up_iter->begin();
                auto down_iter = system_down_iter->begin();
                for (; iter != result_iter->end(); ++iter, ++sys_iter, ++up_iter, ++down_iter)
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
                        throw std::logic_error("spin_component_name should be 'x' or 'y' or 'z' but was : " +
                                               std::to_string(spin_component_name));
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

    template <template <typename> class film_t, typename random_t = qss::random::mersenne::random_t<>>
    result_t perform(nanostructure_type<film_t, proxy_spin> &system) noexcept
    {
        result_t result{};
        auto layers = system.nanostructure;
        const auto amount =
            std::accumulate(layers.begin(), layers.end(), layers.begin()->get_amount_of_nodes(),
                            []([[maybe_unused]] auto first, auto second)
                            { return second.get_amount_of_nodes(); });
        for (auto _ = 0u; _ < amount; ++_)
        {
            const auto coord = layers.template get_random_coord<random_t>();
            double E1 = layers.get_sum_of_closest_neighbours(coord) -
                        layers[static_cast<std::size_t>(coord.film_coord.z)].J * layers.get(coord);

            auto next_coord = coord;
            auto next_film_coord = coord.film_coord;
            auto next_coord_J = layers[coord.idx].J;
            if (next_film_coord.z == layers[coord.idx].get_last_z(coord.film_coord))
            {
                next_film_coord.z = 0;
                next_coord = {static_cast<typename decltype(system.nanostructure)::coords_t::size_type>(coord.idx + 1),
                              next_film_coord};
                next_coord_J = layers.get_J_interlayers()[coord.idx];
            }
            else
            {
                next_film_coord.z += 1;
                next_coord = {coord.idx, next_film_coord};
            }
            double E2{};
            if (coord.film_coord.z == layers[coord.idx].get_last_z(coord.film_coord))
            {
                E2 = 0.0;
            }
            else
            {
                E2 = layers.get_sum_of_closest_neighbours(next_coord) -
                     next_coord_J * layers.get(next_coord);
                E2 -= next_coord_J * layers.get(coord);
                E1 -= next_coord_J * layers.get(next_coord);
            }
            const auto delta_E = E2 - E1;
            static random_t rand{qss::random::get_seed()};
            if (delta_E < 0.0 || rand() < std::exp(-delta_E / system.T))
            {
                auto chosen = layers.get(coord);
                result.up += layers.get(coord).get_up();
                result.down += layers.get(coord).get_down();
                if (next_coord.idx < layers.size())
                {
                    const auto next = layers.get(next_coord);
                    chosen.set_up(typename qss::models::electron_dencity{chosen.get_up() + next.get_up()});
                    chosen.set_down(typename qss::models::electron_dencity{chosen.get_down() + next.get_down()});
                    layers.set(next_coord, chosen);
                }
                chosen.set_up(typename qss::models::electron_dencity{0.0});
                chosen.set_down(typename qss::models::electron_dencity{0.0});
                layers.set(coord, chosen);
            }
        }
        return result;
    }
} // namespace qss::inline algorithms::spin_transport

#endif