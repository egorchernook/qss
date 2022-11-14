#ifndef FCC_HPP_INCLUDED
#define FCC_HPP_INCLUDED

#include <vector>
#include <stdexcept>
#include <utility>
#include <array>
#include <algorithm>

#include "3d.hpp"
#include "../../random/mersenne.hpp"
#include "../../spin.hpp"

namespace qss::lattices::three_d
{
    /* реализует границентрированную решётку.
     * её можно представить как 4 простые кубические
     * со смещением относительно друг друга
     * шаблонный параметр {node_t} -- тип хранимого узла (обычно просто спин нужной модели)
     */
    template <typename node_t, typename size_t = std::uint8_t> // TODO: добавить require для типа node_t
    struct face_centric_cubic : private std::vector<node_t>
    {
        const three_d::sizes_t<size_t> sizes;
        using sizes_t = three_d::sizes_t<size_t>;

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

        /* реализует координаты на ГЦК решётке
         * {w} = 0 -- базовая решётка
         * {w} = 1 -- смещённая в Oxy решётка
         * {w} = 2 -- смещённая в Oyz решётка
         * {w} = 3 -- смещённая в Oxz решётка
         */
        struct coords_t
        {
            std::uint8_t w = 0; // номер простой подрешётки
            std::int16_t x = 0;       //
            std::int16_t y = 0;       // координаты в простой подрешётке
            std::int16_t z = 0;       //
        };

        face_centric_cubic(const node_t &initial_spin, size_t size_x, size_t size_y, size_t size_z)
        {
            sizes = {size_x, size_y, size_z};
            auto temp = size_x * size_y * size_z;
            std::size_t data_length =  temp / 2 + temp % 2;

            this->reserve(data_length);
            std::fill(begin(), end(), initial_spin);
            this->shrink_to_fit();
        }

        face_centric_cubic(const node_t &initial_spin,const  sizes_t &sizes_) 
        : face_centric_cubic_t{initial_spin, sizes_.x, sizes_.y, sizes_.z} {}

        // работает когда есть default параметры конструктора node_t
        face_centric_cubic(const sizes_t &sizes_) 
        : face_centric_cubic_t{node_t{}, sizes_.x, sizes_.y, sizes_.z} {}

        auto get_amount_of_nodes() const // TODO: ограничить auto
        {
            return sizes.x * sizes.y * sizes.z;
        }

        node_t get(const coords_t &coords) const // TODO: добавить проверку на coords
        {
            const auto idx = sizes.w * sizes.x * sizes.y * coords.w 
                           + sizes.x * sizes.y * coords.x 
                           + sizes.y * coords.y 
                           + coords.z;
            return at(idx);
        }
        void set(const node_t &value, const coords_t &coords) // TODO: добавить проверку на coords
        {
            const auto idx = sizes.w * sizes.x * sizes.y * coords.w 
                           + sizes.x * sizes.y * coords.x 
                           + sizes.y * coords.y 
                           + coords.z;
            at(idx) = value;
        }

        /* возвращает массив соседних вершин с координатами
         * для получения суммы учитываете ваши граничные условия
         */
        std::array<coords_t, 12> get_closest_neighbours(const coords_t &coords) const
        {
            const auto central = coords.w * sizes.x * sizes.y * coords.w 
                               + sizes.x * sizes.y * coords.x 
                               + sizes.y * coords.y 
                               + coords.z;

            std::array<coords_t, 12> result{};
            switch (coords.w)
            {
            case 0:
                result[0] = {1, coords.x, coords.y, coords.z};
                result[1] = {1, coords.x, coords.y - 1, coords.z};
                result[2] = {1, coords.x - 1, coords.y, coords.z};
                result[3] = {1, coords.x - 1, coords.y - 1, coords.z};
                result[4] = {2, coords.x, coords.y, coords.z};
                result[5] = {2, coords.x, coords.y - 1, coords.z};
                result[6] = {2, coords.x, coords.y, coords.z - 1};
                result[7] = {2, coords.x, coords.y - 1, coords.z - 1};
                result[8] = {3, coords.x, coords.y, coords.z};
                result[9] = {3, coords.x - 1, coords.y, coords.z};
                result[10] = {3, coords.x, coords.y, coords.z - 1};
                result[11] = {3, coords.x - 1, coords.y, coords.z - 1};
                break;
            case 1:
                result[0] = {0, coords.x, coords.y, coords.z};
                result[1] = {0, coords.x + 1, coords.y, coords.z};
                result[2] = {0, coords.x, coords.y + 1, coords.z};
                result[3] = {0, coords.x + 1, coords.y + 1, coords.z};
                result[4] = {2, coords.x, coords.y, coords.z};
                result[5] = {2, coords.x + 1, coords.y, coords.z};
                result[6] = {2, coords.x, coords.y, coords.z - 1};
                result[7] = {2, coords.x + 1, coords.y, coords.z - 1};
                result[8] = {3, coords.x, coords.y, coords.z};
                result[9] = {3, coords.x, coords.y + 1, coords.z};
                result[10] = {3, coords.x, coords.y, coords.z - 1};
                result[11] = {3, coords.x, coords.y + 1, coords.z - 1};
                break;
            case 2:
                result[0] = {0, coords.x, coords.y, coords.z};
                result[1] = {0, coords.x, coords.y + 1, coords.z};
                result[2] = {0, coords.x, coords.y, coords.z + 1};
                result[3] = {0, coords.x, coords.y + 1, coords.z + 1};
                result[4] = {1, coords.x, coords.y, coords.z};
                result[5] = {1, coords.x - 1, coords.y, coords.z};
                result[6] = {1, coords.x, coords.y, coords.z + 1};
                result[7] = {1, coords.x - 1, coords.y, coords.z + 1};
                result[8] = {3, coords.x, coords.y, coords.z};
                result[9] = {3, coords.x - 1, coords.y, coords.z};
                result[10] = {3, coords.x, coords.y + 1, coords.z};
                result[11] = {3, coords.x - 1, coords.y + 1, coords.z};
                break;
            case 3:
                result[0] = {0, coords.x, coords.y, coords.z};
                result[1] = {0, coords.x + 1, coords.y, coords.z};
                result[2] = {0, coords.x, coords.y, coords.z + 1};
                result[3] = {0, coords.x + 1, coords.y, coords.z + 1};
                result[4] = {1, coords.x, coords.y, coords.z};
                result[5] = {1, coords.x, coords.y - 1, coords.z};
                result[6] = {1, coords.x, coords.y, coords.z + 1};
                result[7] = {1, coords.x, coords.y - 1, coords.z + 1};
                result[8] = {2, coords.x, coords.y, coords.z};
                result[9] = {2, coords.x + 1, coords.y, coords.z};
                result[10] = {2, coords.x, coords.y - 1, coords.z};
                result[11] = {2, coords.x + 1, coords.y - 1, coords.z};
                break;
            default:
                throw std::out_of_range("w must be in range [0-4]");
                break;
            }
            return result;
        }

        template <typename random_t = qss::random::mersenne>
        coords_t choose_random_node()
        {
            auto rand = random_t{};

            return coords_t{
                rand(0, 4),
                rand(0, sizes.x / 4),
                rand(0, sizes.x / 4),
                rand(0, sizes.x / 4)};
        }
    };

    template <typename node_t, typename size_t = std::uint8_t>
    using fcc = face_centric_cubic<node_t, size_t>;
}

#endif