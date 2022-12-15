#ifndef ISING_HPP_INCLUDED
#define ISING_HPP_INCLUDED

#include "../random/random.hpp"
#include "../random/mersenne.hpp"

#include <iostream>
#include <cstdint>
#include <cmath>

namespace qss::models::ising
{
    using magn = double;
    struct spin
    {
        using magn_t = magn;
        std::int8_t value = 1;

        template <Random random_t = qss::random::mersenne::random_t<>>
        static spin generate() noexcept
        {
            static random_t rand{qss::random::get_seed()};
            auto number = rand(0, 2);
            if (number == 0)
            {
                number = -1;
            }
            return spin{static_cast<std::int8_t>(number)};
        }

        operator magn() const noexcept
        {
            return static_cast<magn>(value);
        }
    };

    inline magn operator+(const spin &lhs, const spin &rhs) noexcept
    {
        magn result = lhs.value;
        result += rhs.value;
        return result;
    }
    inline magn operator-(const spin &lhs, const spin &rhs) noexcept
    {
        magn result = lhs.value;
        result -= rhs.value;
        return result;
    }
    inline double scalar_multiply(const spin &lhs, const spin &rhs) noexcept
    {
        magn result = lhs.value;
        result *= rhs.value;
        return result;
    }
    inline magn operator*(const int &lhs, const spin &rhs) noexcept
    {
        magn result = rhs.value;
        result *= lhs;
        return result;
    }
    inline magn operator*(const double &lhs, const spin &rhs) noexcept
    {
        magn result = rhs.value;
        result *= lhs;
        return result;
    }
    inline magn operator*(const spin &lhs, const double &rhs) noexcept
    {
        return (rhs * lhs);
    }

    inline std::ostream &operator<<(std::ostream &out, const spin &data) noexcept
    {
        out << "( " << data.value << " )";
        return out;
    }
    inline std::istream &operator>>(std::istream &in, spin &data) noexcept
    {
        in >> data.value;
        return in;
    }
    inline double abs(const magn &val) noexcept
    {
        return std::abs(val);
    }
}

#endif