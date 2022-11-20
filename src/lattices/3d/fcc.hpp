#ifndef FCC_HPP_INCLUDED
#define FCC_HPP_INCLUDED

#include <vector>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <cassert>
#include <array>

#include "3d.hpp"
#include "../base_lattice.hpp"
#include "../../random/mersenne.hpp"
#include "../../random/random.hpp"

namespace qss::lattices::three_d
{
    /* реализует координаты на ГЦК решётке
     * {w} = 0 -- базовая решётка
     * {w} = 1 -- смещённая в Oxy решётка
     * {w} = 2 -- смещённая в Oyz решётка
     * {w} = 3 -- смещённая в Oxz решётка
     */
    template <typename coord_size_t = std::int16_t>
    struct fcc_coords_t
    {
        std::uint8_t w = 0; // номер простой подрешётки
        coord_size_t x = 0; //
        coord_size_t y = 0; // координаты в простой подрешётке
        coord_size_t z = 0; //
    };

    /* реализует границентрированную решётку.
     * её можно представить как 4 простые кубические
     * со смещением относительно друг друга
     * шаблонный параметр {node_t} -- тип хранимого узла (обычно просто спин нужной модели)
     */
    template <typename node_t,
              typename size_t = std::uint8_t,
              typename coord_size_t = std::int16_t> // TODO: добавить require для типа node_t
    struct face_centric_cubic : public base_lattice_t<node_t, fcc_coords_t<coord_size_t>>
    {
        using base_t = base_lattice_t<node_t, fcc_coords_t<coord_size_t>>;
        using typename base_t::coords_t;
        using typename base_t::value_t;
        const three_d::sizes_t<size_t> sizes;
        using sizes_t = three_d::sizes_t<size_t>;

        const std::array<sizes_t, 4> sublattices_sizes; // размеры подрешёток

        constexpr face_centric_cubic(const value_t &initial_spin,
                                     size_t size_x,
                                     size_t size_y,
                                     size_t size_z)
            : base_t{static_cast<typename base_t::size_type>(
                         size_x * size_y * size_z / 2 + (size_x * size_y * size_z) % 2),
                     initial_spin},
              sizes{size_x, size_y, size_z},
              sublattices_sizes{
                  sizes_t{static_cast<size_t>(size_x / 2 + size_x % 2),
                          static_cast<size_t>(size_y / 2 + size_y % 2),
                          static_cast<size_t>(size_z / 2 + size_z % 2)},
                  sizes_t{static_cast<size_t>(size_x / 2),
                          static_cast<size_t>(size_y / 2),
                          static_cast<size_t>(size_z / 2 + size_z % 2)},
                  sizes_t{static_cast<size_t>(size_x / 2 + size_x % 2),
                          static_cast<size_t>(size_y / 2),
                          static_cast<size_t>(size_z / 2)},
                  sizes_t{static_cast<size_t>(size_x / 2),
                          static_cast<size_t>(size_y / 2 + size_y % 2),
                          static_cast<size_t>(size_z / 2)}}
        {
            this->shrink_to_fit();
        }
        constexpr face_centric_cubic(const value_t &initial_spin, const sizes_t &sizes_)
            : face_centric_cubic{initial_spin, sizes_.x, sizes_.y, sizes_.z} {}
        // работает когда есть default параметры конструктора node_t
        constexpr face_centric_cubic(const sizes_t &sizes_)
            : face_centric_cubic{value_t{}, sizes_.x, sizes_.y, sizes_.z} {}

        value_t get(const coords_t &coords) const
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

            auto get_amount_of_sublattice_nodes = [](const sizes_t &sublattice_size)
                -> typename base_t::size_type
            {
                return sublattice_size.x * sublattice_size.y * sublattice_size.z;
            };
            auto calc_idx = [&coords](const sizes_t &sublattice_size) -> typename base_t::size_type
            {
                return static_cast<typename base_t::size_type>(sublattice_size.x * sublattice_size.y * coords.z + sublattice_size.x * coords.y + coords.x);
            };
            const auto shift = [*this, &get_amount_of_sublattice_nodes](const std::uint8_t &w)
            {
                typename base_t::size_type result{};
                for (auto i = 0u; i < w; ++i)
                {
                    result += get_amount_of_sublattice_nodes(sublattices_sizes[w]);
                }
                return result;
            }(coords.w);
            const typename base_t::size_type idx = shift + calc_idx(sublattices_sizes[coords.w]);
            assert(idx <= this->size());
            return this->at(idx);
        }
        void set(const value_t &value, const coords_t &coords)
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

            auto get_amount_of_sublattice_nodes = [](const sizes_t &sublattice_size)
                -> typename base_t::size_type
            {
                return sublattice_size.x * sublattice_size.y * sublattice_size.z;
            };
            auto calc_idx = [&coords](const sizes_t &sublattice_size) -> typename base_t::size_type
            {
                return static_cast<typename base_t::size_type>(sublattice_size.x * sublattice_size.y * coords.z + sublattice_size.x * coords.y + coords.x);
            };
            const auto shift = [*this, &get_amount_of_sublattice_nodes](const std::uint8_t &w)
            {
                typename base_t::size_type result{};
                for (auto i = 0u; i < w; ++i)
                {
                    result += get_amount_of_sublattice_nodes(sublattices_sizes[w]);
                }
                return result;
            }(coords.w);
            const typename base_t::size_type idx = shift + calc_idx(sublattices_sizes[coords.w]);
            assert(idx <= this->size());
            this->at(idx) = value;
        }

        /* возвращает массив соседних вершин с координатами
         * для получения суммы учитываете ваши граничные условия
         */
        std::vector<coords_t> get_closest_neighbours(const coords_t &coords) const
        {
            std::vector<coords_t> result{};
            switch (coords.w)
            {
            case 0:
                result.push_back({1, coords.x, coords.y, coords.z});
                result.push_back({1, coords.x, static_cast<coord_size_t>(coords.y - 1), coords.z});
                result.push_back({1, static_cast<coord_size_t>(coords.x - 1), coords.y, coords.z});
                result.push_back({1, static_cast<coord_size_t>(coords.x - 1), static_cast<coord_size_t>(coords.y - 1), coords.z});
                result.push_back({2, coords.x, coords.y, coords.z});
                result.push_back({2, coords.x, static_cast<coord_size_t>(coords.y - 1), coords.z});
                result.push_back({2, coords.x, coords.y, static_cast<coord_size_t>(coords.z - 1)});
                result.push_back({2, coords.x, static_cast<coord_size_t>(coords.y - 1), static_cast<coord_size_t>(coords.z - 1)});
                result.push_back({3, coords.x, coords.y, coords.z});
                result.push_back({3, static_cast<coord_size_t>(coords.x - 1), coords.y, coords.z});
                result.push_back({3, coords.x, coords.y, static_cast<coord_size_t>(coords.z - 1)});
                result.push_back({3, static_cast<coord_size_t>(coords.x - 1), coords.y, static_cast<coord_size_t>(coords.z - 1)});
                break;
            case 1:
                result.push_back({0, coords.x, coords.y, coords.z});
                result.push_back({0, static_cast<coord_size_t>(coords.x + 1), coords.y, coords.z});
                result.push_back({0, coords.x, static_cast<coord_size_t>(coords.y + 1), coords.z});
                result.push_back({0, static_cast<coord_size_t>(coords.x + 1), static_cast<coord_size_t>(coords.y + 1), coords.z});
                result.push_back({2, coords.x, coords.y, coords.z});
                result.push_back({2, static_cast<coord_size_t>(coords.x + 1), coords.y, coords.z});
                result.push_back({2, coords.x, coords.y, static_cast<coord_size_t>(coords.z - 1)});
                result.push_back({2, static_cast<coord_size_t>(coords.x + 1), coords.y, static_cast<coord_size_t>(coords.z - 1)});
                result.push_back({3, coords.x, coords.y, coords.z});
                result.push_back({3, coords.x, static_cast<coord_size_t>(coords.y + 1), coords.z});
                result.push_back({3, coords.x, coords.y, static_cast<coord_size_t>(coords.z - 1)});
                result.push_back({3, coords.x, static_cast<coord_size_t>(coords.y + 1), static_cast<coord_size_t>(coords.z - 1)});
                break;
            case 2:
                result.push_back({0, coords.x, coords.y, coords.z});
                result.push_back({0, coords.x, static_cast<coord_size_t>(coords.y + 1), coords.z});
                result.push_back({0, coords.x, coords.y, static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({0, coords.x, static_cast<coord_size_t>(coords.y + 1), static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({1, coords.x, coords.y, coords.z});
                result.push_back({1, static_cast<coord_size_t>(coords.x - 1), coords.y, coords.z});
                result.push_back({1, coords.x, coords.y, static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({1, static_cast<coord_size_t>(coords.x - 1), coords.y, static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({3, coords.x, coords.y, coords.z});
                result.push_back({3, static_cast<coord_size_t>(coords.x - 1), coords.y, coords.z});
                result.push_back({3, coords.x, static_cast<coord_size_t>(coords.y + 1), coords.z});
                result.push_back({3, static_cast<coord_size_t>(coords.x - 1), static_cast<coord_size_t>(coords.y + 1), coords.z});
                break;
            case 3:
                result.push_back({0, coords.x, coords.y, coords.z});
                result.push_back({0, static_cast<coord_size_t>(coords.x + 1), coords.y, coords.z});
                result.push_back({0, coords.x, coords.y, static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({0, static_cast<coord_size_t>(coords.x + 1), coords.y, static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({1, coords.x, coords.y, coords.z});
                result.push_back({1, coords.x, static_cast<coord_size_t>(coords.y - 1), coords.z});
                result.push_back({1, coords.x, coords.y, static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({1, coords.x, static_cast<coord_size_t>(coords.y - 1), static_cast<coord_size_t>(coords.z + 1)});
                result.push_back({2, coords.x, coords.y, coords.z});
                result.push_back({2, static_cast<coord_size_t>(coords.x + 1), coords.y, coords.z});
                result.push_back({2, coords.x, static_cast<coord_size_t>(coords.y - 1), coords.z});
                result.push_back({2, static_cast<coord_size_t>(coords.x + 1), static_cast<coord_size_t>(coords.y - 1), coords.z});
                break;
            default:
                throw std::out_of_range("coords.w out of range : " + std::to_string(coords.w));
                break;
            }
            return result;
        }

        template <typename random_t = qss::random::mersenne::random_t<>>
        coords_t choose_random_node() const
        {
            static auto rand = random_t{qss::random::get_seed()};
            const auto w = static_cast<std::uint8_t>(rand(0, 4));
            return coords_t{w, static_cast<coord_size_t>(rand(0, sublattices_sizes[w].x)),
                            static_cast<coord_size_t>(rand(0, sublattices_sizes[w].y)),
                            static_cast<coord_size_t>(rand(0, sublattices_sizes[w].z))};
        }
    };

    template <typename node_t, typename size_t = std::uint8_t>
    using fcc = face_centric_cubic<node_t, size_t>;
}

#endif