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
    concept ThreeD_Lattice = requires std::same_as<typename lattice_t::sizes_t, qss::lattices::three_d::sizes_t>;

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

    public:
        multilayer(std::initializer_list<lattice_t> list,
                   std::vector<double> &&J_films_,
                   std::vector<double> &&J_interlayers_)
        {
            if (J_films_.size() != list.size())
            {
                throw std::logic_error("you must give exchange integrals for all film : " + std::to_string(J_films_.size()) + " != " + std::to_string(list.size()));
            }
            if (J_interlayer_.size() != list.size() - 1)
            {
                throw std::logic_error("you must give exchange integrals for all interfilms : " + std::to_string(J_interlayers_.size()) + " != " + std::to_string(list.size() - 1));
            }
            J_films = std::move(J_films_);
            J_interlayers = std::move(J_interlayers_);
            for (auto &&elem : list)
            {
                check_sizes(elem.sizes);
                sizes.z += elem_sizes.z;
                films_thickness.push_back(elem_sizes.z);
                this->push_back(std::move(elem));
            }
        }
    }
}

#endif