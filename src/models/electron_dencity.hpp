#ifndef ELECTRON_DENCITY_HPP_INCLUDED
#define ELECTRON_DENCITY_HPP_INCLUDED

namespace qss::models
{
    struct electron_dencity
    {
        using magn_t = double;
        double value;

        operator magn_t() const
        {
            return static_cast<magn_t>(value);
        }
    };

}

#endif