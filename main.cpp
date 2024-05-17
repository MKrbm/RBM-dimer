#include <iostream>
#include <vector>
#include <test.hpp>
#include <lattice/basis.hpp>
#include <lattice/graph_xml.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

int main() {
    std::ifstream is("../lattices.xml");
    if (!is) {
        std::cerr << "Error opening file lattices.xml" << std::endl;
        return 1;
    }

    boost::property_tree::ptree pt;
    boost::property_tree::read_xml(is, pt);

    std::string basis_name = "chain lattice";
    std::string cell_name = "simple1d";

    lattice::basis bs;
    read_xml(pt, basis_name, bs);

    lattice::unitcell cell;
    read_xml(pt, cell_name, cell);

    std::cout << cell.dimension() << std::endl;

    return 0;
}
