#ifndef FCC_HPP_INCLUDED
#define FCC_HPP_INCLUDED

#include <vector>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <cassert>

#include "3d.hpp"
#include "../base_lattice.hpp"
#include "../../random/mersenne.hpp"
#include "../../random/random.hpp"

namespace qss::lattices::three_d
{
    /*
     * реализует координаты на ГЦК решётке
     * {w} = 0 -- базовая решётка
     * {w} = 1 -- смещённая в Oxy решётка
     * {w} = 2 -- смещённая в Oyz решётка
     * {w} = 3 -- смещённая в Oxz решётка
     **/
    struct fcc_coords_t
    {
        using size_type = int;
        std::uint8_t w = 0; // номер простой подрешётки
        size_type x = 0;    //
        size_type y = 0;    // координаты в простой подрешётке
        size_type z = 0;    //
    };

    inline std::vector<fcc_coords_t> get_closest_neighbours(const fcc_coords_t &coords)
    {
        std::vector<fcc_coords_t> result{};
        switch (coords.w)
        {
        case 0:
            return {
                {1, coords.x, coords.y, coords.z},
                {1, coords.x, coords.y - 1, coords.z},
                {1, coords.x - 1, coords.y, coords.z},
                {1, coords.x - 1, coords.y - 1, coords.z},
                {2, coords.x, coords.y, coords.z},
                {2, coords.x, coords.y - 1, coords.z},
                {2, coords.x, coords.y, coords.z - 1},
                {2, coords.x, coords.y - 1, coords.z - 1},
                {3, coords.x, coords.y, coords.z},
                {3, coords.x - 1, coords.y, coords.z},
                {3, coords.x, coords.y, coords.z - 1},
                {3, coords.x - 1, coords.y, coords.z - 1}};
            break;
        case 1:
            return {
                {0, coords.x, coords.y, coords.z},
                {0, coords.x + 1, coords.y, coords.z},
                {0, coords.x, coords.y + 1, coords.z},
                {0, coords.x + 1, coords.y + 1, coords.z},
                {2, coords.x, coords.y, coords.z},
                {2, coords.x + 1, coords.y, coords.z},
                {2, coords.x, coords.y, coords.z - 1},
                {2, coords.x + 1, coords.y, coords.z - 1},
                {3, coords.x, coords.y, coords.z},
                {3, coords.x, coords.y + 1, coords.z},
                {3, coords.x, coords.y, coords.z - 1},
                {3, coords.x, coords.y + 1, coords.z - 1}};
            break;
        case 2:
            return {
                {0, coords.x, coords.y, coords.z},
                {0, coords.x, coords.y + 1, coords.z},
                {0, coords.x, coords.y, coords.z + 1},
                {0, coords.x, coords.y + 1, coords.z + 1},
                {1, coords.x, coords.y, coords.z},
                {1, coords.x - 1, coords.y, coords.z},
                {1, coords.x, coords.y, coords.z + 1},
                {1, coords.x - 1, coords.y, coords.z + 1},
                {3, coords.x, coords.y, coords.z},
                {3, coords.x - 1, coords.y, coords.z},
                {3, coords.x, coords.y + 1, coords.z},
                {3, coords.x - 1, coords.y + 1, coords.z}};
            break;
        case 3:
            return {
                {0, coords.x, coords.y, coords.z},
                {0, coords.x + 1, coords.y, coords.z},
                {0, coords.x, coords.y, coords.z + 1},
                {0, coords.x + 1, coords.y, coords.z + 1},
                {1, coords.x, coords.y, coords.z},
                {1, coords.x, coords.y - 1, coords.z},
                {1, coords.x, coords.y, coords.z + 1},
                {1, coords.x, coords.y - 1, coords.z + 1},
                {2, coords.x, coords.y, coords.z},
                {2, coords.x + 1, coords.y, coords.z},
                {2, coords.x, coords.y - 1, coords.z},
                {2, coords.x + 1, coords.y - 1, coords.z}};
            break;
        default:
            throw std::out_of_range("coords.w out of range : " + std::to_string(coords.w));
            break;
        }
        return result;
    }

    /*
     * реализует границентрированную решётку.
     * её можно представить как 4 простые кубические
     * со смещением относительно друг друга
     * шаблонный параметр {node_t} -- тип хранимого узла (обычно просто спин нужной модели)
     **/
    template <typename node_t> // TODO: добавить require для типа node_t
    struct face_centric_cubic : public base_lattice_t<node_t, fcc_coords_t>
    {
        using base_t = base_lattice_t<node_t, fcc_coords_t>;
        using typename base_t::coords_t;
        using typename base_t::value_t;
        const three_d::sizes_t sizes;
        using sizes_t = three_d::sizes_t;

    private:
        void bounds_check(const coords_t &coords) const
        {
            if (coords.w >= 4)
            {
                throw std::out_of_range("coords.w out of range : " + std::to_string(coords.w));
            }
            if (coords.x < 0 || coords.x >= sublattices_sizes[coords.w].x)
            {
                throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
            }
            if (coords.y < 0 || coords.y >= sublattices_sizes[coords.w].y)
            {
                throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
            }
            if (coords.z < 0 || coords.z >= sublattices_sizes[coords.w].z)
            {
                throw std::out_of_range("coords.z out of range : " + std::to_string(coords.z));
            }
        }
        typename base_t::size_type get_amount_of_sublattice_nodes(const sizes_t &sublattice_size) const
        {
            return sublattice_size.x * sublattice_size.y * sublattice_size.z;
        }
        typename base_t::size_type calc_idx(const sizes_t &sublattice_size, const coords_t &coords) const
        {
            return static_cast<typename base_t::size_type>(sublattice_size.x * sublattice_size.y * coords.z + sublattice_size.x * coords.y + coords.x);
        }
        typename base_t::size_type calc_shift(const std::uint8_t &w) const
        {
            typename base_t::size_type result{};
            for (auto i = 0u; i < w; ++i)
            {
                result += get_amount_of_sublattice_nodes(sublattices_sizes[w]);
            }
            return result;
        };

    public:
        const std::array<sizes_t, 4> sublattices_sizes; // размеры подрешёток

        constexpr face_centric_cubic(const value_t &initial_spin,
                                     const typename sizes_t::size_type &size_x,
                                     const typename sizes_t::size_type &size_y,
                                     const typename sizes_t::size_type &size_z)
            : base_t(static_cast<typename base_t::size_type>(
                         size_x * size_y * size_z / 2 + (size_x * size_y * size_z) % 2),
                     initial_spin),
              sizes{size_x, size_y, size_z},
              sublattices_sizes{
                  sizes_t{static_cast<typename sizes_t::size_type>(size_x / 2 + size_x % 2),
                          static_cast<typename sizes_t::size_type>(size_y / 2 + size_y % 2),
                          static_cast<typename sizes_t::size_type>(size_z / 2 + size_z % 2)},
                  sizes_t{static_cast<typename sizes_t::size_type>(size_x / 2),
                          static_cast<typename sizes_t::size_type>(size_y / 2),
                          static_cast<typename sizes_t::size_type>(size_z / 2 + size_z % 2)},
                  sizes_t{static_cast<typename sizes_t::size_type>(size_x / 2 + size_x % 2),
                          static_cast<typename sizes_t::size_type>(size_y / 2),
                          static_cast<typename sizes_t::size_type>(size_z / 2)},
                  sizes_t{static_cast<typename sizes_t::size_type>(size_x / 2),
                          static_cast<typename sizes_t::size_type>(size_y / 2 + size_y % 2),
                          static_cast<typename sizes_t::size_type>(size_z / 2)}}
        {
        }
        constexpr face_centric_cubic(const value_t &initial_spin, const sizes_t &sizes_)
            : face_centric_cubic{initial_spin, sizes_.x, sizes_.y, sizes_.z} {}
        // работает когда есть default параметры конструктора node_t
        constexpr face_centric_cubic(const sizes_t &sizes_)
            : face_centric_cubic{value_t{}, sizes_.x, sizes_.y, sizes_.z} {}

        value_t get(const coords_t &coords) const
        {
            bounds_check(coords);
            const typename base_t::size_type idx = calc_shift(coords.w) + calc_idx(sublattices_sizes[coords.w], coords);
            assert(idx <= this->size());
            return this->at(idx);
        }
        void set(const value_t &value, const coords_t &coords)
        {
            bounds_check(coords);
            const typename base_t::size_type idx = calc_shift(coords.w) + calc_idx(sublattices_sizes[coords.w], coords);
            assert(idx <= this->size());
            this->at(idx) = value;
        }

        template <typename random_t = qss::random::mersenne::random_t<>>
        coords_t choose_random_node() const
        {
            using coord_size_t = typename coords_t::size_type;
            static auto rand = random_t{qss::random::get_seed()};
            const auto w = static_cast<std::uint8_t>(rand(0, 4));
            return coords_t{w, static_cast<coord_size_t>(rand(0, sublattices_sizes[w].x)),
                            static_cast<coord_size_t>(rand(0, sublattices_sizes[w].y)),
                            static_cast<coord_size_t>(rand(0, sublattices_sizes[w].z))};
        }
    };

    template <typename node_t>
    using fcc = face_centric_cubic<node_t>;
}

#endif