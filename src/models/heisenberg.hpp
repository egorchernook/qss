#ifndef HEISENBERG_HPP_INLCUDE
#define HEISENBERG_HPP_INLCUDE

#include <iostream>

#include "../random.hpp"
#include "../random/mersenne.hpp"

namespace qss::models::heisenberg
{
    struct spin
    {
        double x = 1.0;
        double y = 0.0;
        double z = 0.0;

        spin &operator+=(const spin &rhs)
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        spin &operator-=(const spin &rhs)
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }

        static spin create(const double &phi = 0.0, const double &eta = 0.0)
        {
            const double sinus_eta = std::sin(eta / 2);

            spin result{};
            result.x = sinus_eta * std::cos(phi / 2);
            result.y = sinus_eta * std::sin(phi / 2);
            result.z = std::cos(eta / 2);
            
            return result;
        }

        template<Random random_t = qss::random::mersenne::random_t<>>
        static spin generate()
        {
            random_t rand{qss::random::get_seed()};
            return create(rand.get_angle_2pi(), rand.get_angle_pi());
        }
    };

    spin operator+(const spin &lhs, const spin &rhs)
    {
        spin result{};
        result.x = lhs.x + rhs.x;
        result.y = lhs.y + rhs.y;
        result.z = lhs.z + rhs.z;
        return result;
    }

    spin operator-(const spin &lhs, const spin &rhs)
    {
        spin result{};
        result.x = lhs.x - rhs.x;
        result.y = lhs.y - rhs.y;
        result.z = lhs.z - rhs.z;
        return result;
    }

    double scalar_multiply(const spin &lhs, const spin &rhs)
    {
        double sum = 0.0;
        sum += lhs.x * rhs.x;
        sum += lhs.y * rhs.y;
        sum += lhs.z * rhs.z;
        return sum;
    }

    spin operator*(const double &lhs, const spin &rhs)
    {
        spin result;
        result.x = lhs * rhs.x;
        result.y = lhs * rhs.y;
        result.z = lhs * rhs.z;
        return (result);
    }
    spin operator*(const spin &lhs, const double &rhs)
    {
        return (rhs * lhs);
    }

    std::ostream &operator<<(std::ostream &out, const spin &data)
    {
        out << "( " << data.x << " , " << data.y << " , " << data.z << " )";
        return out;
    }
    std::istream &operator>>(std::istream &in, spin &data)
    {
        in >> data.x >> data.y >> data.z;
        return in;
    }
}
#endif