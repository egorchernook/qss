#ifndef SQUARE_HPP_INCLUDED
#define SQUARE_HPP_INCLUDED

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <optional>

#include "2d.hpp"
#include "../base_lattice.hpp"
#include "../../random/random.hpp"
#include "../../random/mersenne.hpp"

namespace qss::lattices::two_d
{
    template <typename coord_size_t = std::int16_t>
    struct square_coords_t
    {
        coord_size_t x = 0;
        coord_size_t y = 0;
    };

    template <typename node_t,
              typename size_t = std::uint8_t,
              typename coord_size_t = std::int16_t> // TODO: добавить require для типа node_t
    struct square : public base_lattice_t<node_t, square_coords_t<coord_size_t>>
    {
        using base_t = base_lattice_t<node_t, square_coords_t<coord_size_t>>;
        using typename base_t::value_t;
        using typename base_t::coords_t;
        const two_d::sizes_t<size_t> sizes;
        using sizes_t = two_d::sizes_t<size_t>;
        
        constexpr square(const value_t &initial_spin, const size_t &size_x, const size_t &size_y)
            : base_t{static_cast<typename base_t::size_type>(size_x * size_y), initial_spin}, sizes{size_x, size_y}
        {
            this->shrink_to_fit();
        }
        constexpr square(const value_t &initial_spin, const sizes_t &sizes_)
            : square{initial_spin, sizes_.x, sizes_.y} {}
        // работает когда есть default параметры конструктора node_t
        constexpr explicit square(const sizes_t &sizes_)
            : square{value_t{}, sizes_.x, sizes_.y} {}

        value_t get(const coords_t &coords) const
        {
            if (coords.x < 0 || coords.x >= sizes.x)
            {
                throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
            }
            if (coords.y < 0 || coords.y >= sizes.y)
            {
                throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
            }
            const auto idx = sizes.x * coords.y + coords.x;
            return this->at(idx);
        }
        void set(const value_t &value, const coords_t &coords)
        {
            if (coords.x < 0 || coords.x >= sizes.x)
            {
                throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
            }
            if (coords.y < 0 || coords.y >= sizes.y)
            {
                throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
            }
            const auto idx = sizes.x * coords.y + coords.x;
            this->at(idx) = value;
        }

        std::vector<coords_t> get_closest_neighbours(const coords_t &coords) const
        {
            return {coords_t{static_cast<coord_size_t>(coords.x - 1), coords.y},
                    coords_t{coords.x, static_cast<coord_size_t>(coords.y - 1)},
                    coords_t{static_cast<coord_size_t>(coords.x + 1), coords.y},
                    coords_t{coords.x, static_cast<coord_size_t>(coords.y + 1)}};
        }

        template <typename random_t = qss::random::mersenne::random_t<>>
        coords_t choose_random_node() const
        {
            static auto rand = random_t{qss::random::get_seed()};
            return coords_t{
                static_cast<coord_size_t>(rand(0, sizes.x)),
                static_cast<coord_size_t>(rand(0, sizes.y))};
        }
    };
}

#endif