#ifndef MULTILAYER_HPP_INCLUDED
#define MULTILAYER_HPP_INCLUDED

#include <vector>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <concepts>
#include <algorithm>
#include <utility>

#include "../lattices/3d/3d.hpp"
#include "../lattices/borders_conditions.hpp"
#include "../utility/functions.hpp"

namespace qss::nanostructures
{
    template <typename lattice_t>
    concept ThreeD_Lattice = std::same_as<typename lattice_t::sizes_t, qss::lattices::three_d::sizes_t>;

    template <ThreeD_Lattice lattice_t>
    struct film : public lattice_t
    {
        using xy_border_condition = typename qss::border_conditions::periodic;
        using z_border_condition = typename qss::border_conditions::sharp;
        double J;
    };

    template <ThreeD_Lattice lattice_t>
    typename film<lattice_t>::value_t::magn_t
    get_sum_of_closest_neighbours(const film<lattice_t> &film,
                                  const typename film<lattice_t>::coords_t &central)
    {
        using film_t = typename film<lattice_t>;
        return J * get_sum_of_closest_neighbours(film, central,
                                                 qss::border_conditions::use_border_conditions<
                                                     film_t::xy_border_condition,
                                                     film_t::xy_border_condition,
                                                     film_t::z_border_condition>);
    }

    template <ThreeD_Lattice lattice_t>
    class multilayer : private std::vector<film<lattice_t>>
    {
        using container_t = typename std::vector<film<lattice_t>>;
        container_t films;
        std::vector<double> J_interlayers; // обменные интегралы взаимодействий плёнок.

        void check_sizes(const sizes_t &value) const
        {
            if (value.x != sizes.x)
            {
                throw std::logic_error("x size must be the same for all lattices : " + std::to_string(value.x) + " != " + std::to_string(sizes.x));
            }
            if (value.y != sizes.y)
            {
                throw std::logic_error("y size must be the same for all lattices : " + std::to_string(value.y) + " != " + std::to_string(sizes.y));
            }
        }

    public:
        using film_t = typename film<lattice_t>;

        constexpr multilayer(const film_t &film) : films{film} {};
        constexpr multilayer(film_t &&film) : films{std::move(film)} {};

        void add(const film_t &film, double J_interlayer)
        {
            check_sizes(film.sizes);
            films.push_back(film);
            J_interlayers.push_back(J_interlayer);
        }
        void add(film_t &&film, double J_interlayer)
        {
            check_sizes(film.sizes);
            films.push_back(std::move(film));
            J_interlayers.push_back(J_interlayer);
        }
        film_t pop()
        {
            film_t result = films.back();
            films.pop_back();
            J_interlayers.pop_back();
            return result;
        }
        const film_t &get(container_t::size_type idx)
        {
            return films.at(idx);
        }

        
    };
}

#endif