#ifndef SPIN_HPP_INCLUDED
#define SPIN_HPP_INCLUDED

#include <concepts>
#include <type_traits>
namespace qss
{
    template <typename T>
    concept Spin = requires(T a, T b) {
        { T::generate() } -> std::convertible_to<T>;
    };
}
#endif