#include "../include/graph.hpp"
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <lattice/basis.hpp>
#include <lattice/graph_xml.hpp>

namespace lattice
{

    std::tuple<std::vector<bond_t>, std::vector<std::vector<double>>> create_graph(const std::string &file, const std::string &lattice_name, const std::string &cell_name, const std::vector<size_t> &size, const boundary_t &boundary)
    {
        std::ifstream is(file);
        if (!is.is_open())
        {
            throw std::runtime_error("Failed to open lattice configuration file: " + file);
        }

        boost::property_tree::ptree pt;
        boost::property_tree::read_xml(is, pt);

        lattice::basis bs;
        read_xml(pt, lattice_name, bs);

        lattice::unitcell cell;
        read_xml(pt, cell_name, cell);

        if (size.size() != cell.dimension())
        {
            throw std::invalid_argument("Size vector does not match the lattice dimension.");
        }

        lattice::span_t extent;
        switch (size.size())
        {
        case 1:
            extent = lattice::extent(size[0]);
            break;
        case 2:
            extent = lattice::extent(size[0], size[1]);
            break;
        case 3:
            extent = lattice::extent(size[0], size[1], size[2]);
            break;
        case 4:
            extent = lattice::extent(size[0], size[1], size[2], size[3]);
            break;
        case 5:
            extent = lattice::extent(size[0], size[1], size[2], size[3], size[4]);
            break;
        default:
            throw std::invalid_argument("Unsupported number of dimensions for lattice extent.");
        }

        lattice::graph lat(bs, cell, extent, boundary);

        std::vector<bond_t> bonds;
        std::vector<std::vector<double>> coordinates;

        for (size_t i = 0; i < lat.num_bonds(); i++)
        {
            bonds.push_back({{lat.source(i), lat.target(i)}, lat.bond_type(i)});
        }

        for (size_t i = 0; i < lat.num_sites(); i++)
        {
            const Eigen::VectorXd &coord = lat.coordinate(i);
            std::vector<double> coord_vec(coord.data(), coord.data() + coord.size());
            coordinates.push_back(coord_vec);
        }

        return std::make_tuple(bonds, coordinates);
    }

} // namespace lattice