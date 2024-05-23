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
            os << "Bond Type: " << b.bond_type << ", Bonds: {";
            for (size_t i = 0; i < b.bonds.size(); ++i) {
                os << b.bonds[i];
                if (i < b.bonds.size() - 1) os << ", ";
            }
            os << "}";
            return os;
        }
    };

    struct loop_t
    {
        std::vector<std::size_t> loops;
        int loop_type;

        friend std::ostream& operator<<(std::ostream& os, const loop_t& l) {
            os << "Loop Type: " << l.loop_type << ", Loops: {";
            for (size_t i = 0; i < l.loops.size(); ++i) {
                os << l.loops[i];
                if (i < l.loops.size() - 1) os << ", ";
            }
            os << "}";
            return os;
        }
    };

    std::tuple<std::vector<bond_t>, std::vector<loop_t>, std::vector<std::vector<double>>> create_graph(const std::string &file, const std::string &lattice_name, const std::string &cell_name, const std::vector<size_t> &size, const boundary_t &boundary);

} // namespace lattice

