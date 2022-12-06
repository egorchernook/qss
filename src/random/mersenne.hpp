#ifndef MERSENNE_HPP_INCLUDED
#define MERSENNE_HPP_INCLUDED

#include <cmath>
#include <random>
#include <type_traits>
#include <chrono>

namespace qss::random::mersenne
{

#ifndef M_PI
#define M_PI 3.1'415'926'535'8979'323'846
#endif

    template <typename genrand_t = std::mt19937> // TODO: добавить ограничение на typename
    class random_t
    {
        static constexpr const typename genrand_t::result_type seed_default = 0;
        genrand_t m_genrand{};

    public:
        random_t() noexcept : m_genrand(random_t<genrand_t>::seed_default) {}

        [[nodiscard]] random_t(const random_t<genrand_t> &) noexcept = default;
        [[nodiscard]] random_t(random_t<genrand_t> &&) noexcept = default;

        [[nodiscard]] random_t(const typename genrand_t::result_type seed) noexcept : m_genrand(seed) {}

        ~random_t() noexcept = default;

        random_t<genrand_t> &operator=(const random_t<genrand_t> &) noexcept = default;
        random_t<genrand_t> &operator=(random_t<genrand_t> &&) noexcept = default;

        void cooldown(const typename genrand_t::result_type seed) noexcept
        {
            this->m_genrand
                .seed(static_cast<typename genrand_t::result_type>(seed));
        }
        //возвращает double в полуинтервале [0;1)
        double operator()() noexcept
        {
            constexpr double temp_denominator =
                1.0 / (static_cast<double>(genrand_t::max() - genrand_t::min()) + 1.0);

            return temp_denominator * static_cast<double>(m_genrand() - genrand_t::min());
        }
        //возвращает double в полуинтервале [_begin;_end)
        double operator()(const double _begin, const double _end) noexcept
        {
            return _begin + (_end - _begin) * (*this)();
        }
        //возвращает int в полуинтервале [_begin;_end)
        int operator()(const int _begin, const int _end) noexcept
        {
            return _begin + static_cast<int>((_end - _begin) * (*this)());
        }
        //возвращает double в полуинтервале [0; 2pi)
        double get_angle_2pi() noexcept
        {
            constexpr double angle = 2 * M_PI;

            constexpr double temp_denominator =
                1.0 / (static_cast<double>(genrand_t::max() - genrand_t::min()) + 1.0);

            return angle * temp_denominator *
                   static_cast<double>(m_genrand() - genrand_t::min());
        }
        //возвращает double в полуинтервале [0; pi)
        double get_angle_pi() noexcept
        {
            constexpr double angle = M_PI;

            constexpr double temp_denominator =
                1.0 / (static_cast<double>(genrand_t::max() - genrand_t::min()) + 1.0);

            return angle * temp_denominator *
                   static_cast<double>(m_genrand() - genrand_t::min());
        }
    };
}

#endif