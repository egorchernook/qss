#ifndef RANDOM_HPP_INCLUDED
#define RANDOM_HPP_INCLUDED

// #include <concepts>
#include <chrono>
#include <cmath>
#include <stdexcept>

#define Random typename // to migrate to c++17

namespace qss
{
    // template <typename T>
    // concept Random = requires (T rand) {
    //     { rand() } -> std::convertible_to<double>;
    //     { rand(0.0, 1.0) } -> std::convertible_to<double>;
    //     { rand(0, 2) } -> std::convertible_to<int>;
    //     { rand.get_angle_2pi() } -> std::convertible_to<double>;
    //     { rand.get_angle_pi() } -> std::convertible_to<double>;
    // };

    namespace random 
    {
        inline std::size_t get_seed(const std::size_t init = 0, const std::size_t top_limiter = 2'004'991) noexcept
        {
            static std::size_t counter = 0;

            using std::chrono::nanoseconds;
            using clock_type = std::chrono::system_clock;
            using time_type = std::chrono::time_point<clock_type, nanoseconds>;

            const auto time_start = time_type{nanoseconds{init}};
            const auto time_end = static_cast<time_type>(clock_type::now());

            const auto time_interval =
                static_cast<std::size_t>(std::chrono::duration_cast<nanoseconds>(time_end - time_start).count());

            const auto seed = (time_interval + counter) % top_limiter;

            counter += init + 1;

            if (counter > top_limiter)
            {
                counter -= top_limiter;
            }

            return seed;
        }
    }
}
#endif