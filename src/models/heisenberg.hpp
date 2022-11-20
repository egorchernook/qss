#ifndef HEISENBERG_HPP_INLCUDE
#define HEISENBERG_HPP_INLCUDE

#include <iostream>

#include "../random/random.hpp"
#include "../random/mersenne.hpp"

namespace qss::models::heisenberg
{
    struct magn
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        magn &operator+=(const magn &rhs)
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }

        magn &operator-=(const magn &rhs)
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }
    };
    struct spin
    {
        double x = 1.0;
        double y = 0.0;
        double z = 0.0;

        template <Random random_t = qss::random::mersenne::random_t<>>
        static spin generate() noexcept
        {
            static random_t rand{qss::random::get_seed()};
            const double phi = rand.get_angle_2pi();
            const double eta = rand.get_angle_pi();
            const double sinus_eta = std::sin(eta / 2);

            spin result{};
            result.x = sinus_eta * std::cos(phi / 2);
            result.y = sinus_eta * std::sin(phi / 2);
            result.z = std::cos(eta / 2);

            return result;
        }

        operator magn() const
        {
            return magn{x, y, z};
        }
    };

    magn operator+(const spin &lhs, const spin &rhs)
    {
        magn result{};
        result.x = lhs.x + rhs.x;
        result.y = lhs.y + rhs.y;
        result.z = lhs.z + rhs.z;
        return result;
    }
    magn operator+(const magn &lhs, const magn &rhs)
    {
        magn result{};
        result.x = lhs.x + rhs.x;
        result.y = lhs.y + rhs.y;
        result.z = lhs.z + rhs.z;
        return result;
    }
    magn operator-(const spin &lhs, const spin &rhs)
    {
        magn result{};
        result.x = lhs.x - rhs.x;
        result.y = lhs.y - rhs.y;
        result.z = lhs.z - rhs.z;
        return result;
    }
    magn operator-(const magn &lhs, const magn &rhs)
    {
        magn result{};
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
    double scalar_multiply(const magn &lhs, const magn &rhs)
    {
        double sum = 0.0;
        sum += lhs.x * rhs.x;
        sum += lhs.y * rhs.y;
        sum += lhs.z * rhs.z;
        return sum;
    }
    magn operator*(const double &lhs, const spin &rhs)
    {
        magn result;
        result.x = lhs * rhs.x;
        result.y = lhs * rhs.y;
        result.z = lhs * rhs.z;
        return (result);
    }
    magn operator*(const double &lhs, const magn &rhs)
    {
        magn result;
        result.x = lhs * rhs.x;
        result.y = lhs * rhs.y;
        result.z = lhs * rhs.z;
        return (result);
    }
    magn operator*(const spin &lhs, const double &rhs)
    {
        return (rhs * lhs);
    }
    magn operator*(const magn &lhs, const double &rhs)
    {
        return (rhs * lhs);
    }
    magn operator/(const spin &lhs, const double &rhs){
        magn result{};
        result.x = lhs.x / rhs;
        result.y = lhs.y / rhs;
        result.z = lhs.z / rhs;
        return result;
    }
    magn operator/(const magn &lhs, const double &rhs){
        magn result{};
        result.x = lhs.x / rhs;
        result.y = lhs.y / rhs;
        result.z = lhs.z / rhs;
        return result;
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
    std::ostream &operator<<(std::ostream &out, const magn &data)
    {
        out << data.x << "\t" << data.y << "\t" << data.z;
        return out;
    }
    std::istream &operator>>(std::istream &in, magn &data)
    {
        in >> data.x >> data.y >> data.z;
        return in;
    }
    double abs(const magn &val)
    {
        return std::sqrt( val.x * val.x + val.y * val.y + val.z * val.z);
    }
}
#endif