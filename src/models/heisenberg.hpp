#ifndef HEISENBERG_HPP_INLCUDE
#define HEISENBERG_HPP_INLCUDE

#include <iostream>

#include "../random/random.hpp"
#include "../random/mersenne.hpp"

namespace qss::inline models::heisenberg
{
    struct magn
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;

        magn &operator+=(const magn &rhs) noexcept
        {
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
            return *this;
        }
        magn &operator-=(const magn &rhs) noexcept
        {
            x -= rhs.x;
            y -= rhs.y;
            z -= rhs.z;
            return *this;
        }
        magn &operator*=(const double &rhs) noexcept
        {
            x *= rhs;
            y *= rhs;
            z *= rhs;
            return *this;
        }
        magn &operator/=(const double &rhs) noexcept
        {
            x /= rhs;
            y /= rhs;
            z /= rhs;
            return *this;
        }
    };
    struct spin
    {
        using magn_t = magn;
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

        operator magn_t() const noexcept
        {
            return magn_t{x, y, z};
        }
    };

    inline magn operator+(const spin &lhs, const spin &rhs) noexcept
    {
        magn result{};
        result.x = lhs.x + rhs.x;
        result.y = lhs.y + rhs.y;
        result.z = lhs.z + rhs.z;
        return result;
    }
    inline magn operator+(const magn &lhs, const magn &rhs) noexcept
    {
        magn result{};
        result.x = lhs.x + rhs.x;
        result.y = lhs.y + rhs.y;
        result.z = lhs.z + rhs.z;
        return result;
    }
    inline magn operator-(const spin &lhs, const spin &rhs) noexcept
    {
        magn result{};
        result.x = lhs.x - rhs.x;
        result.y = lhs.y - rhs.y;
        result.z = lhs.z - rhs.z;
        return result;
    }
    inline magn operator-(const magn &lhs, const magn &rhs) noexcept
    {
        magn result{};
        result.x = lhs.x - rhs.x;
        result.y = lhs.y - rhs.y;
        result.z = lhs.z - rhs.z;
        return result;
    }
    inline double scalar_multiply(const spin &lhs, const spin &rhs) noexcept
    {
        double sum = 0.0;
        sum += lhs.x * rhs.x;
        sum += lhs.y * rhs.y;
        sum += lhs.z * rhs.z;
        return sum;
    }
    inline double scalar_multiply(const magn &lhs, const magn &rhs) noexcept
    {
        double sum = 0.0;
        sum += lhs.x * rhs.x;
        sum += lhs.y * rhs.y;
        sum += lhs.z * rhs.z;
        return sum;
    }
    inline magn operator*(const double &lhs, const spin &rhs) noexcept
    {
        magn result;
        result.x = lhs * rhs.x;
        result.y = lhs * rhs.y;
        result.z = lhs * rhs.z;
        return (result);
    }
    inline magn operator*(const double &lhs, const magn &rhs) noexcept
    {
        magn result;
        result.x = lhs * rhs.x;
        result.y = lhs * rhs.y;
        result.z = lhs * rhs.z;
        return (result);
    }
    inline magn operator*(const spin &lhs, const double &rhs) noexcept
    {
        return (rhs * lhs);
    }
    inline magn operator*(const magn &lhs, const double &rhs) noexcept
    {
        return (rhs * lhs);
    }
    inline magn operator/(const spin &lhs, const double &rhs) noexcept
    {
        magn result{};
        result.x = lhs.x / rhs;
        result.y = lhs.y / rhs;
        result.z = lhs.z / rhs;
        return result;
    }
    inline magn operator/(const magn &lhs, const double &rhs) noexcept
    {
        magn result{};
        result.x = lhs.x / rhs;
        result.y = lhs.y / rhs;
        result.z = lhs.z / rhs;
        return result;
    }
    inline std::ostream &operator<<(std::ostream &out, const spin &data) noexcept
    {
        out << "( " << data.x << " , " << data.y << " , " << data.z << " )";
        return out;
    }
    inline std::istream &operator>>(std::istream &in, spin &data) noexcept
    {
        in >> data.x >> data.y >> data.z;
        return in;
    }
    inline std::ostream &operator<<(std::ostream &out, const magn &data) noexcept
    {
        out << data.x << "\t" << data.y << "\t" << data.z;
        return out;
    }
    inline std::istream &operator>>(std::istream &in, magn &data) noexcept
    {
        in >> data.x >> data.y >> data.z;
        return in;
    }
    inline double abs(const magn &val) noexcept
    {
        return std::sqrt(val.x * val.x + val.y * val.y + val.z * val.z);
    }
}
#endif