#ifndef BASE_LATTICE_HPP_INCLUDED
#define BASE_LATTICE_HPP_INCLUDED

#include <concepts>
#include <vector>

namespace qss::lattices
{
    template <typename node_t, typename coordinates_t>
    struct base_lattice_t : protected std::vector<node_t>
    {
        using container_t = std::vector<node_t>;
        using container_t::container_t;
        using container_t::begin;
        using container_t::cbegin;
        using container_t::rbegin;
        using container_t::crbegin;
        using container_t::end;
        using container_t::cend;
        using container_t::rend;
        using container_t::crend;

        using value_t = node_t;
        using coords_t = coordinates_t;
        constexpr typename container_t::size_type get_amount_of_nodes() const
        {
            return this->size();
        }

        virtual ~base_lattice_t() {};
        virtual value_t get(const coords_t &coord) const = 0;
        virtual void set(const node_t &value, const coords_t &coords) = 0;
        virtual std::vector<coords_t> get_closest_neighbours(const coords_t &coords) const = 0;
    };

    template <typename T>
    concept Lattice = std::is_base_of_v<
        base_lattice_t<typename T::value_t, typename T::coords_t>,
        T> && requires(T obj)
    {
        {
            obj.choose_random_node()
            } -> std::convertible_to<typename T::coords_t>;
    };

}

#endif