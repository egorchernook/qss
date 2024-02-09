#ifndef ISING_HPP_INCLUDED
#define ISING_HPP_INCLUDED

#include "../random/mersenne.hpp"
#include "../random/random.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <sstream>

namespace qss {
inline namespace models {
namespace ising {
using magn = double;
struct spin {
    using magn_t = magn;
    std::int8_t value = 1;

    template<Random random_t = qss::random::mersenne::random_t<>>
    static spin generate() noexcept
    {
        static random_t rand{qss::random::get_seed()};
        auto number = rand(0, 2);
        if (number == 0) {
            number = -1;
        }
        return spin{static_cast<std::int8_t>(number)};
    }

    operator magn() const noexcept
    {
        return static_cast<magn>(value);
    }
};

inline magn operator+(const spin& lhs, const spin& rhs) noexcept
{
    magn result = lhs.value;
    result += rhs.value;
    return result;
}
inline magn operator-(const spin& lhs, const spin& rhs) noexcept
{
    magn result = lhs.value;
    result -= rhs.value;
    return result;
}
inline double scalar_multiply(const spin& lhs, const spin& rhs) noexcept
{
    magn result = lhs.value;
    result *= rhs.value;
    return result;
}
inline magn operator*(const int& lhs, const spin& rhs) noexcept
{
    magn result = rhs.value;
    result *= lhs;
    return result;
}
inline magn operator*(const double& lhs, const spin& rhs) noexcept
{
    magn result = rhs.value;
    result *= lhs;
    return result;
}
inline magn operator*(const spin& lhs, const double& rhs) noexcept
{
    return (rhs * lhs);
}

inline bool is_almost_equals(const magn& fst, const magn& snd, double epsilon = 1e-4) noexcept
{
    const auto diff = fst - snd;
    return diff < epsilon;
}
inline std::string to_string(const magn& data) noexcept
{
    std::ostringstream stream{};
    stream << "( " << data << " )";
    return stream.str();
}
inline std::ostream& operator<<(std::ostream& out, const spin& data) noexcept
{
    out << "( " << data.value << " )";
    return out;
}
inline std::istream& operator>>(std::istream& in, spin& data) noexcept
{
    in >> data.value;
    return in;
}
inline double abs(const magn& val) noexcept
{
    return std::abs(val);
}
} // namespace ising
} // namespace models
} // namespace qss
#endif