#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <lattice/graph.hpp>

template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
  os << "[ ";
  for (const T &item : vec)
    os << item << ", ";
  os << "]";
  return os;
}

namespace lattice
{

    struct bond_t
    {
        std::vector<std::size_t> bonds;
        int bond_type;

        friend std::ostream& operator<<(std::ostream& os, const bond_t& b) {
            os << "Bond Type: " << b.bond_type << " Bonds: [ ";
            for (auto& bond : b.bonds) {
                os << bond << ", ";
            }
            os << "]";
            return os;
        }
    };

    std::tuple<std::vector<bond_t>, std::vector<std::vector<double>>> create_graph(const std::string &file, const std::string &lattice_name, const std::string &cell_name, const std::vector<size_t> &size, const boundary_t &boundary);

} // namespace lattice

