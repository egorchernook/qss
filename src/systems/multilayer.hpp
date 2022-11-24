#ifndef MULTILAYER_HPP_INCLUDED
#define MULTILAYER_HPP_INCLUDED

#include <vector>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <concepts>
#include <algorithm>

#include "../lattices/3d/3d.hpp"
#include "../lattices/borders_conditions.hpp"

namespace qss::nanostructures
{
    template <typename lattice_t>
    concept ThreeD_Lattice = std::same_as<typename lattice_t::sizes_t, qss::lattices::three_d::sizes_t>;

    template <ThreeD_Lattice lattice_t>
    struct multilayer : public lattice_t
    {
        using base_t = lattice_t;
        using xy_border_condition = typename qss::border_conditions::periodic;

    private:
        std::vector<typename base_t::sizes_t::size_type> films_thickness;
        std::vector<double> J_films;       // обменные интегралы плёнок.
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
        typename base_t::base_t::container_t::size_type get_amount_of_node_in_film(const typename base_t::sizes_t::size_type &number) const
        {
            return films_thickness[number] * sizes.x * sizes.y;
        }

    public:
        constexpr multilayer(lattice_t &&first_lattice, double first_film_J)
            : base_t{std::move(first_lattice)},
              sizes{first_lattice.sizes},
              films_thickness{first_lattice.sizes.z},
              J_films{first_film_J},
              J_interlayers{} {}

        void add(lattice_t &&film, double J_film, double J_interlayer)
        {
            check_sizes(film.sizes);
            this->reserve(this->size() + film.get_amount_of_nodes());
            for(auto elem&& : film){
                this->push_back(std::move(elem));
            }
            films_thickness.push_back(std::move(lattice.sizes.z));
            J_films.push_back(J_film);
            J_interlayers.push_back(J_interlayer);
        }
    };
}

#endif