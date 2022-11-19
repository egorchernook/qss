#ifndef FCC_HPP_INCLUDED
#define FCC_HPP_INCLUDED

#include <vector>
#include <stdexcept>
#include <utility>
#include <algorithm>

#include "3d.hpp"
#include "../../random/mersenne.hpp"
#include "../../random.hpp"

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
    struct face_centric_cubic : private std::vector<node_t>
    {
        const three_d::sizes_t<size_t> sizes;
        using sizes_t = three_d::sizes_t<size_t>;

        const sizes_t base_sublattice_size; // у базовой подрешётки и подрешётке, смещённой в Oxy, одинаковые размеры
        // у подрешёток, смещённых в Oyz и Oxz, одинаковые размеры и их легче считать

        using value_t = node_t;
        using container_t = std::vector<node_t>;
        using container_t::begin;
        using container_t::cbegin;
        using container_t::cend;
        using container_t::crbegin;
        using container_t::crend;
        using container_t::end;
        using container_t::rbegin;
        using container_t::rend;
        using coords_t = fcc_coords_t<coord_size_t>;

        constexpr face_centric_cubic(const node_t &initial_spin,
                                     size_t size_x,
                                     size_t size_y,
                                     size_t size_z)
            : container_t{static_cast<container_t::size_type>(
                              size_x * size_y * size_z / 2 + (size_x * size_y * size_z) % 2),
                          initial_spin},
              sizes{size_x, size_y, size_z},
              base_sublattice_size{static_cast<size_t>(size_x / 2 + size_x % 2),
                                   static_cast<size_t>(size_y / 2 + size_y % 2),
                                   static_cast<size_t>(size_z / 2 + size_z % 2)}
        {
            this->shrink_to_fit();
        }
        constexpr face_centric_cubic(const node_t &initial_spin, const sizes_t &sizes_)
            : face_centric_cubic{initial_spin, sizes_.x, sizes_.y, sizes_.z} {}
        // работает когда есть default параметры конструктора node_t
        constexpr face_centric_cubic(const sizes_t &sizes_)
            : face_centric_cubic{node_t{}, sizes_.x, sizes_.y, sizes_.z} {}

        std::size_t get_amount_of_nodes() const
        {
            return this->size();
        }

        node_t get(const coords_t &coords) const
        {
            switch (coords.w)
            {
            case 0:
                [[fallthrough]];
            case 1:
                if (coords.x < 0 || coords.x >= base_sublattice_size.x)
                {
                    throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
                }
                if (coords.y < 0 || coords.y >= base_sublattice_size.y)
                {
                    throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
                }
                if (coords.z < 0 || coords.z >= base_sublattice_size.z)
                {
                    throw std::out_of_range("coords.z out of range : " + std::to_string(coords.z));
                }
                break;
            case 2:
                [[fallthrough]];
            case 3:
                if (coords.x < 0 || coords.x >= sizes.x / 2)
                {
                    throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
                }
                if (coords.y < 0 || coords.y >= sizes.y / 2)
                {
                    throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
                }
                if (coords.z < 0 || coords.z >= sizes.z / 2)
                {
                    throw std::out_of_range("coords.z out of range : " + std::to_string(coords.z));
                }
                break;
            default:
                throw std::out_of_range("coords.w out of range : " + std::to_string(coords.w));
                break;
            }
            const auto amount_of_base_lattice_nodes = base_sublattice_size.x * base_sublattice_size.y * base_sublattice_size.z;
            auto idx = base_sublattice_size.x * base_sublattice_size.y * coords.x + base_sublattice_size.x * coords.y + coords.z;
            switch (coords.w)
            {
            case 0:
                break;
            case 1:
                idx += amount_of_base_lattice_nodes;
                break;
            case 2:
                idx += amount_of_base_lattice_nodes * 2;
                break;
            case 3:
                idx += amount_of_base_lattice_nodes * 2 + sizes.x * sizes.y * sizes.z / 8;
                break;
            default:
                break;
            }
            return this->at(idx);
        }
        void set(const node_t &value, const coords_t &coords)
        {
            switch (coords.w)
            {
            case 0:
                [[fallthrough]];
            case 1:
                if (coords.x < 0 || coords.x >= base_sublattice_size.x)
                {
                    throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
                }
                if (coords.y < 0 || coords.y >= base_sublattice_size.y)
                {
                    throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
                }
                if (coords.z < 0 || coords.z >= base_sublattice_size.z)
                {
                    throw std::out_of_range("coords.z out of range : " + std::to_string(coords.z));
                }
                break;
            case 2:
                [[fallthrough]];
            case 3:
                if (coords.x < 0 || coords.x >= sizes.x / 2)
                {
                    throw std::out_of_range("coords.x out of range : " + std::to_string(coords.x));
                }
                if (coords.y < 0 || coords.y >= sizes.y / 2)
                {
                    throw std::out_of_range("coords.y out of range : " + std::to_string(coords.y));
                }
                if (coords.z < 0 || coords.z >= sizes.z / 2)
                {
                    throw std::out_of_range("coords.z out of range : " + std::to_string(coords.z));
                }
                break;
            default:
                throw std::out_of_range("coords.w out of range : " + std::to_string(coords.w));
                break;
            }
            const auto amount_of_base_lattice_nodes = base_sublattice_size.x * base_sublattice_size.y * base_sublattice_size.z;
            auto idx = base_sublattice_size.x * base_sublattice_size.y * coords.x + base_sublattice_size.x * coords.y + coords.z;
            switch (coords.w)
            {
            case 0:
                break;
            case 1:
                idx += amount_of_base_lattice_nodes;
                break;
            case 2:
                idx += amount_of_base_lattice_nodes * 2;
                break;
            case 3:
                idx += amount_of_base_lattice_nodes * 2 + sizes.x * sizes.y * sizes.z / 8;
                break;
            default:
                break;
            }
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
            switch (w)
            {
            case 0:
                [[fallthrough]];
            case 1:
                return coords_t{w, static_cast<coord_size_t>(rand(0, base_sublattice_size.x)),
                                static_cast<coord_size_t>(rand(0, base_sublattice_size.y)),
                                static_cast<coord_size_t>(rand(0, base_sublattice_size.z))};
                break;
            case 2:
                [[fallthrough]];
            case 3:
                return coords_t{w, static_cast<coord_size_t>(rand(0, sizes.x / 2)),
                                static_cast<coord_size_t>(rand(0, sizes.y / 2)),
                                static_cast<coord_size_t>(rand(0, sizes.z / 2))};
                break;
            default:
                throw std::out_of_range("coords.w out of range : " + std::to_string(w));
                break;
            }
        }
    };

    template <typename node_t, typename size_t = std::uint8_t>
    using fcc = face_centric_cubic<node_t, size_t>;
}

#endif