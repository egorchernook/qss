#ifndef ISING_HPP_INCLUDED
#define ISING_HPP_INCLUDED

#include "../random.hpp"
#include "../random/mersenne.hpp"

#include <iostream>
#include <cstdint>

namespace qss::models::ising
{
     using magn = double;

    struct spin
    {
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

        operator magn() const {
            return static_cast<magn>(value);
        }
    };

    magn operator+(const spin &lhs, const spin &rhs)
    {
        magn result = lhs.value;
        result += rhs.value;
        return result;
    }
    magn operator-(const spin &lhs, const spin &rhs)
    {
        magn result = lhs.value;
        result -= rhs.value;
        return result;
    }
    magn scalar_multiply(const spin &lhs, const spin &rhs)
    {
        magn result = lhs.value;
        result *= rhs.value;
        return result;
    }
    magn operator*(const int &lhs, const spin &rhs)
    {
        magn result = rhs.value;
        result *= lhs;
        return result;
    }
    magn operator*(const spin &lhs, const int &rhs)
    {
        return (rhs * lhs);
    }
    
    std::ostream &operator<<(std::ostream &out, const spin &data)
    {
        out << "( " << data.value << " )";
        return out;
    }
    std::istream &operator>>(std::istream &in, spin &data)
    {
        in >> data.value;
        return in;
    }
}

#endif