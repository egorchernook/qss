#ifndef SQUARE_HPP_INCLUDED
#define SQUARE_HPP_INCLUDED

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include "2d.hpp"
#include "../base_lattice.hpp"
#include "../../random/random.hpp"
#include "../../random/mersenne.hpp"

namespace qss::lattices::two_d
{
    struct square_coords_t
    {
        using size_type = int;
        size_type x = 0;
        size_type y = 0;
    };

    std::vector<square_coords_t> get_closest_neighbours(const square_coords_t &coords)
    {
        std::vector<square_coords_t> result{};
        result.push_back({coords.x - 1, coords.y});
        result.push_back({coords.x, coords.y - 1});
        result.push_back({coords.x + 1, coords.y});
        result.push_back({coords.x, coords.y + 1});
        return result;
    }

    template <typename node_t> // TODO: добавить require для типа node_t
    struct square : public base_lattice_t<node_t, square_coords_t>
    {
        using base_t = base_lattice_t<node_t, square_coords_t>;
        using typename base_t::coords_t;
        using typename base_t::value_t;
        const two_d::sizes_t sizes;
        using sizes_t = two_d::sizes_t;

    private:
        void bounds_check(const coords_t &coords) const
        {
            if (coords.x < 0 || coords.x >= sizes.x)
            {
                throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
            }
            if (coords.y < 0 || coords.y >= sizes.y)
            {
                throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
            }
        }

    public:
        constexpr square(const value_t &initial_spin,
                         const typename sizes_t::size_type &size_x,
                         const typename sizes_t::size_type &size_y)
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
            bounds_check(coords);
            const auto idx = static_cast<typename base_t::size_type>(sizes.x * coords.y + coords.x);
            return this->at(idx);
        }
        void set(const value_t &value, const coords_t &coords)
        {
            bounds_check(coords);
            const auto idx = static_cast<typename base_t::size_type>(sizes.x * coords.y + coords.x);
            this->at(idx) = value;
        }

        template <typename random_t = qss::random::mersenne::random_t<>>
        coords_t choose_random_node() const
        {
            static auto rand = random_t{qss::random::get_seed()};
            return coords_t{
                static_cast<typename coords_t::size_type>(rand(0, sizes.x)),
                static_cast<typename coords_t::size_type>(rand(0, sizes.y))};
        }
    };
}

#endif