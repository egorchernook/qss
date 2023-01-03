#ifndef ELECTRON_DENCITY_HPP_INCLUDED
#define ELECTRON_DENCITY_HPP_INCLUDED

#include <iostream>

namespace qss::inline models
{
    struct electron_dencity
    {
        using magn_t = double;
        double value;

        operator magn_t() const noexcept
        {
            return static_cast<magn_t>(value);
        }
    };

    inline typename electron_dencity::magn_t operator+(const electron_dencity &lhs,
                                                       const electron_dencity &rhs) noexcept
    {
        typename electron_dencity::magn_t result = lhs.value;
        result += rhs.value;
        return result;
    }
    inline typename electron_dencity::magn_t operator-(const electron_dencity &lhs,
                                                       const electron_dencity &rhs) noexcept
    {
        typename electron_dencity::magn_t result = lhs.value;
        result -= rhs.value;
        return result;
    }
    inline double scalar_multiply(const electron_dencity &lhs,
                                  const electron_dencity &rhs) noexcept
    {
        typename electron_dencity::magn_t result = lhs.value;
        result *= rhs.value;
        return result;
    }
    inline typename electron_dencity::magn_t operator*(const int &lhs,
                                                       const electron_dencity &rhs) noexcept
    {
        typename electron_dencity::magn_t result = rhs.value;
        result *= lhs;
        return result;
    }
    inline typename electron_dencity::magn_t operator*(const double &lhs,
                                                       const electron_dencity &rhs) noexcept
    {
        typename electron_dencity::magn_t result = rhs.value;
        result *= lhs;
        return result;
    }
    inline typename electron_dencity::magn_t operator*(const electron_dencity &lhs,
                                                       const double &rhs) noexcept
    {
        return (rhs * lhs);
    }

    inline std::ostream &operator<<(std::ostream &out,
                                    const electron_dencity &data) noexcept
    {
        out << "( " << data.value << " )";
        return out;
    }
    inline std::istream &operator>>(std::istream &in,
                                    electron_dencity &data) noexcept
    {
        in >> data.value;
        return in;
    }
    inline double abs(const typename electron_dencity::magn_t &val) noexcept
    {
        return std::abs(val);
    }

}

#endif